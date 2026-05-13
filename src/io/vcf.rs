//! Shared VCF/BCF output-format and indexing helpers.

use rust_htslib::htslib;
use std::ffi::CString;

/// Output container/encoding inferred from an output path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputKind {
    PlainVcf,
    BgzfVcf,
    Bcf,
}

impl OutputKind {
    /// Stable label used in summaries and docs.
    pub fn as_str(self) -> &'static str {
        match self {
            OutputKind::PlainVcf => "vcf",
            OutputKind::BgzfVcf => "vcf.gz",
            OutputKind::Bcf => "bcf",
        }
    }

    /// Whether htslib can build a random-access variant index for this output.
    pub fn is_indexable(self) -> bool {
        matches!(self, OutputKind::BgzfVcf | OutputKind::Bcf)
    }
}

/// Variant index type to build after writing an output file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputIndexKind {
    Csi,
    Tbi,
}

impl OutputIndexKind {
    /// Stable label used in summaries and docs.
    pub fn as_str(self) -> &'static str {
        match self {
            OutputIndexKind::Csi => "csi",
            OutputIndexKind::Tbi => "tbi",
        }
    }

    fn min_shift(self) -> i32 {
        match self {
            OutputIndexKind::Csi => 14,
            OutputIndexKind::Tbi => 0,
        }
    }
}

/// User-facing indexing policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputIndexPolicy {
    /// Index indexable file outputs by default.
    Auto,
    /// Do not index outputs.
    Off,
    /// Force a specific index kind.
    On(OutputIndexKind),
}

/// Infer VCF/BCF output kind from an optional output path.
pub fn infer_output_kind(output_path: Option<&str>) -> OutputKind {
    match output_path {
        None | Some("-") => OutputKind::PlainVcf,
        Some(path) => {
            let lower = path.to_ascii_lowercase();
            if lower.ends_with(".bcf") {
                OutputKind::Bcf
            } else if lower.ends_with(".vcf.gz") || lower.ends_with(".vcf.bgz") {
                OutputKind::BgzfVcf
            } else {
                OutputKind::PlainVcf
            }
        }
    }
}

fn unchecked_resolved_output_index(
    output_path: Option<&str>,
    policy: OutputIndexPolicy,
) -> Option<OutputIndexKind> {
    match policy {
        OutputIndexPolicy::Off => None,
        OutputIndexPolicy::On(kind) => Some(kind),
        OutputIndexPolicy::Auto => match (output_path, infer_output_kind(output_path)) {
            (Some(path), kind) if path != "-" && kind.is_indexable() => Some(OutputIndexKind::Csi),
            _ => None,
        },
    }
}

/// Validate an explicit user index request.
///
/// Automatic indexing silently does nothing for stdout/plain VCF; explicit
/// requests keep failing loudly so mistakes are visible.
pub fn validate_index_request(
    output_path: Option<&str>,
    policy: OutputIndexPolicy,
) -> Result<(), String> {
    let OutputIndexPolicy::On(index_kind) = policy else {
        return Ok(());
    };
    let Some(path) = output_path.filter(|p| *p != "-") else {
        return Err(
            "--write-index requires -o/--output with .vcf.gz, .vcf.bgz, or .bcf output".to_string(),
        );
    };
    match (infer_output_kind(output_path), index_kind) {
        (OutputKind::PlainVcf, _) => Err(format!(
            "--write-index cannot index plain VCF output '{path}'; use .vcf.gz, .vcf.bgz, or .bcf"
        )),
        (OutputKind::Bcf, OutputIndexKind::Tbi) => {
            Err("--write-index=tbi is only valid for BGZF VCF; BCF output requires CSI".to_string())
        }
        (OutputKind::BgzfVcf | OutputKind::Bcf, OutputIndexKind::Csi)
        | (OutputKind::BgzfVcf, OutputIndexKind::Tbi) => Ok(()),
    }
}

/// Resolve the effective index kind for an output path and policy.
///
/// This is the safe resolver for callers: explicit invalid requests return a
/// friendly error before htslib indexing is attempted.
pub fn resolve_output_index(
    output_path: Option<&str>,
    policy: OutputIndexPolicy,
) -> Result<Option<OutputIndexKind>, String> {
    validate_index_request(output_path, policy)?;
    Ok(unchecked_resolved_output_index(output_path, policy))
}

fn index_error_message(code: i32) -> &'static str {
    match code {
        -1 => "indexing failed",
        -2 => "opening output failed",
        -3 => "format not indexable",
        -4 => "failed to create or save the index",
        _ => "unknown indexing error",
    }
}

/// Build a CSI/TBI index for a completed VCF.GZ/VCF.BGZ/BCF output file.
pub fn build_index(path: &str, index_kind: OutputIndexKind, threads: usize) -> Result<(), String> {
    let c_path = CString::new(path.as_bytes())
        .map_err(|_| format!("output path contains NUL byte: {path}"))?;
    let threads = threads.min(i32::MAX as usize) as i32;
    let ret = unsafe {
        htslib::bcf_index_build3(
            c_path.as_ptr(),
            std::ptr::null(),
            index_kind.min_shift(),
            threads,
        )
    };
    if ret == 0 {
        Ok(())
    } else {
        Err(format!(
            "failed to build {} index for '{}': {} (htslib code {})",
            index_kind.as_str(),
            path,
            index_error_message(ret),
            ret
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn infers_output_kind_from_suffix() {
        assert_eq!(infer_output_kind(None), OutputKind::PlainVcf);
        assert_eq!(infer_output_kind(Some("-")), OutputKind::PlainVcf);
        assert_eq!(infer_output_kind(Some("out.vcf")), OutputKind::PlainVcf);
        assert_eq!(infer_output_kind(Some("out.vcf.gz")), OutputKind::BgzfVcf);
        assert_eq!(infer_output_kind(Some("out.vcf.bgz")), OutputKind::BgzfVcf);
        assert_eq!(infer_output_kind(Some("out.bcf")), OutputKind::Bcf);
    }

    #[test]
    fn resolves_auto_index_policy_only_for_indexable_files() {
        assert_eq!(
            resolve_output_index(Some("out.vcf.gz"), OutputIndexPolicy::Auto).unwrap(),
            Some(OutputIndexKind::Csi)
        );
        assert_eq!(
            resolve_output_index(Some("out.bcf"), OutputIndexPolicy::Auto).unwrap(),
            Some(OutputIndexKind::Csi)
        );
        assert_eq!(
            resolve_output_index(Some("out.vcf"), OutputIndexPolicy::Auto).unwrap(),
            None
        );
        assert_eq!(
            resolve_output_index(None, OutputIndexPolicy::Auto).unwrap(),
            None
        );
    }

    #[test]
    fn explicit_policy_validates_bad_requests() {
        assert!(resolve_output_index(None, OutputIndexPolicy::On(OutputIndexKind::Csi)).is_err());
        assert!(
            resolve_output_index(Some("out.vcf"), OutputIndexPolicy::On(OutputIndexKind::Csi))
                .is_err()
        );
        assert!(
            resolve_output_index(Some("out.bcf"), OutputIndexPolicy::On(OutputIndexKind::Tbi))
                .is_err()
        );
        assert_eq!(
            resolve_output_index(
                Some("out.vcf.gz"),
                OutputIndexPolicy::On(OutputIndexKind::Tbi)
            )
            .unwrap(),
            Some(OutputIndexKind::Tbi)
        );
    }

    #[test]
    fn off_policy_disables_indexing() {
        assert_eq!(
            resolve_output_index(Some("out.vcf.gz"), OutputIndexPolicy::Off).unwrap(),
            None
        );
    }

    #[test]
    fn build_index_rejects_nul_path_before_htslib_call() {
        let err = build_index("bad\0path.vcf.gz", OutputIndexKind::Csi, 1).unwrap_err();
        assert!(err.contains("NUL byte"));
    }
}
