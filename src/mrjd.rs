//! Initial multi-region joint-detection kernels.
//!
//! This module is a small, auditable foundation for high-homology / multi-region
//! candidate evidence. It is not a DRAGEN-equivalent caller: it groups candidate
//! SNV evidence by region-group offset so downstream code can inspect homologous
//! support across multiple loci before a full posterior haplotype model exists.

use crate::io::fasta::Fai;
use rust_htslib::bam::{self, Read};
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub group: String,
    pub chrom: String,
    pub start1: i64,
    pub end1: i64,
    pub copy: String,
}

impl Region {
    pub fn len(&self) -> i64 {
        self.end1 - self.start1 + 1
    }

    fn contains_pos0(&self, pos0: i64) -> bool {
        (self.start1 - 1) <= pos0 && pos0 < self.end1
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CandidateConfig {
    pub min_mapq: u8,
    pub min_baseq: u8,
    pub min_alt_count: u32,
    pub min_alt_fraction: f64,
}

impl Default for CandidateConfig {
    fn default() -> Self {
        Self {
            min_mapq: 0,
            min_baseq: 13,
            min_alt_count: 2,
            min_alt_fraction: 0.20,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RegionObservation {
    pub copy: String,
    pub chrom: String,
    pub pos1: i64,
    pub ref_base: u8,
    pub depth: u32,
    pub alt_count: u32,
    pub alt_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct JointCandidate {
    pub group: String,
    pub offset1: i64,
    pub alt_base: u8,
    pub alt_positive_depth: u32,
    pub alt_positive_alt_count: u32,
    pub regions_with_alt: usize,
    pub region_count: usize,
    pub observations: Vec<RegionObservation>,
}

fn is_regions_header(fields: &[&str]) -> bool {
    fields.len() >= 4
        && fields[0].eq_ignore_ascii_case("group")
        && fields[1].eq_ignore_ascii_case("chrom")
        && fields[2].eq_ignore_ascii_case("start")
        && fields[3].eq_ignore_ascii_case("end")
        && fields
            .get(4)
            .map_or(true, |value| value.eq_ignore_ascii_case("copy"))
}

pub fn read_regions_tsv(path: &str) -> Result<Vec<Region>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open regions TSV '{path}': {e}"))?;
    let reader = BufReader::new(file);
    let mut regions = Vec::new();
    let mut first_record = true;
    for (line_no, line) in reader.lines().enumerate() {
        let line_no = line_no + 1;
        let line = line.map_err(|e| format!("failed reading regions TSV '{path}': {e}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let fields = trimmed.split('\t').collect::<Vec<_>>();
        if first_record && is_regions_header(&fields) {
            first_record = false;
            continue;
        }
        first_record = false;
        if fields.len() < 4 || fields.len() > 5 {
            return Err(format!(
                "regions TSV line {line_no} must have 4 or 5 tab-separated fields: group chrom start end [copy]"
            ));
        }
        let group = fields[0].to_string();
        let chrom = fields[1].to_string();
        let start1 = fields[2]
            .parse::<i64>()
            .map_err(|_| format!("regions TSV line {line_no} has invalid start"))?;
        let end1 = fields[3]
            .parse::<i64>()
            .map_err(|_| format!("regions TSV line {line_no} has invalid end"))?;
        if group.is_empty() || chrom.is_empty() {
            return Err(format!(
                "regions TSV line {line_no} has empty group or chrom"
            ));
        }
        if start1 < 1 || end1 < start1 {
            return Err(format!(
                "regions TSV line {line_no} has invalid interval {chrom}:{start1}-{end1}"
            ));
        }
        let copy = fields
            .get(4)
            .filter(|value| !value.is_empty())
            .map(|value| (*value).to_string())
            .unwrap_or_else(|| format!("{}:{}-{}", chrom, start1, end1));
        regions.push(Region {
            group,
            chrom,
            start1,
            end1,
            copy,
        });
    }
    if regions.is_empty() {
        return Err(format!("regions TSV '{path}' did not contain any regions"));
    }
    Ok(regions)
}

pub fn validate_candidate_config(cfg: CandidateConfig) -> Result<(), String> {
    if cfg.min_alt_count == 0 {
        return Err("min_alt_count must be >= 1".to_string());
    }
    if !(0.0..=1.0).contains(&cfg.min_alt_fraction) {
        return Err("min_alt_fraction must be between 0 and 1".to_string());
    }
    Ok(())
}

fn base_index(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn base_from_index(index: usize) -> u8 {
    [b'A', b'C', b'G', b'T'][index]
}

fn region_count_for_offset(regions: &[Region], group: &str, offset1: i64) -> usize {
    regions
        .iter()
        .filter(|region| region.group == group && offset1 >= 1 && offset1 <= region.len())
        .count()
}

fn usable_record(record: &bam::Record, min_mapq: u8) -> bool {
    !record.is_unmapped()
        && !record.is_quality_check_failed()
        && !record.is_duplicate()
        && !record.is_secondary()
        && !record.is_supplementary()
        && (min_mapq == 0 || (record.mapq() != 255 && record.mapq() >= min_mapq))
}

pub fn detect_snv_candidates(
    bam_path: &str,
    reference_path: &str,
    fai: &Fai,
    regions: &[Region],
    cfg: CandidateConfig,
    threads: usize,
) -> Result<Vec<JointCandidate>, String> {
    validate_candidate_config(cfg)?;
    if regions.is_empty() {
        return Ok(Vec::new());
    }
    let mut bam = bam::IndexedReader::from_path(bam_path)
        .map_err(|e| format!("cannot open indexed BAM/CRAM '{bam_path}': {e}"))?;
    bam.set_reference(reference_path)
        .map_err(|e| format!("failed to set CRAM reference '{reference_path}': {e}"))?;
    if threads > 1 {
        bam.set_threads(threads)
            .map_err(|e| format!("failed to enable BAM/CRAM threads for '{bam_path}': {e}"))?;
    }

    let mut grouped: BTreeMap<(String, i64, u8), Vec<RegionObservation>> = BTreeMap::new();
    for region in regions {
        let ref_seq = fai.fetch_seq(&region.chrom, region.start1, region.end1)?;
        if ref_seq.len() != region.len() as usize {
            return Err(format!(
                "reference interval {}:{}-{} returned length {}, expected {}",
                region.chrom,
                region.start1,
                region.end1,
                ref_seq.len(),
                region.len()
            ));
        }
        let start0 = region.start1 - 1;
        let end0 = region.end1;
        bam.fetch((region.chrom.as_bytes(), start0, end0))
            .map_err(|e| {
                format!(
                    "failed to fetch {}:{}-{} from '{bam_path}': {e}",
                    region.chrom, region.start1, region.end1
                )
            })?;
        for pileup in bam.pileup() {
            let pileup =
                pileup.map_err(|e| format!("failed to read pileup from '{bam_path}': {e}"))?;
            let pos0 = pileup.pos() as i64;
            if !region.contains_pos0(pos0) {
                continue;
            }
            let offset0 = pos0 - start0;
            let Some(&ref_base) = ref_seq.get(offset0 as usize) else {
                continue;
            };
            let Some(ref_index) = base_index(ref_base) else {
                continue;
            };
            let mut counts = [0u32; 4];
            for alignment in pileup.alignments() {
                let record = alignment.record();
                if !usable_record(&record, cfg.min_mapq) {
                    continue;
                }
                let Some(qpos) = alignment.qpos() else {
                    continue;
                };
                if record.qual().get(qpos).copied().unwrap_or(0) < cfg.min_baseq {
                    continue;
                }
                let base = record.seq()[qpos];
                if let Some(index) = base_index(base) {
                    counts[index] += 1;
                }
            }
            let depth = counts.iter().sum::<u32>();
            if depth == 0 {
                continue;
            }
            for (index, &alt_count) in counts.iter().enumerate() {
                if index == ref_index || alt_count < cfg.min_alt_count {
                    continue;
                }
                let alt_fraction = alt_count as f64 / depth as f64;
                if alt_fraction < cfg.min_alt_fraction {
                    continue;
                }
                let offset1 = offset0 + 1;
                grouped
                    .entry((region.group.clone(), offset1, base_from_index(index)))
                    .or_default()
                    .push(RegionObservation {
                        copy: region.copy.clone(),
                        chrom: region.chrom.clone(),
                        pos1: pos0 + 1,
                        ref_base: ref_base.to_ascii_uppercase(),
                        depth,
                        alt_count,
                        alt_fraction,
                    });
            }
        }
    }

    let mut out = Vec::new();
    for ((group, offset1, alt_base), observations) in grouped {
        let alt_positive_depth = observations.iter().map(|obs| obs.depth).sum();
        let alt_positive_alt_count = observations.iter().map(|obs| obs.alt_count).sum();
        let region_count = region_count_for_offset(regions, &group, offset1);
        out.push(JointCandidate {
            group,
            offset1,
            alt_base,
            alt_positive_depth,
            alt_positive_alt_count,
            regions_with_alt: observations.len(),
            region_count,
            observations,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_region_manifest() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!("phase_tools_regions_{}.tsv", std::process::id()));
        std::fs::write(
            &path,
            "group\tchrom\tstart\tend\tcopy\nG1\tchr1\t10\t20\tcopy1\n",
        )
        .unwrap();
        let regions = read_regions_tsv(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].group, "G1");
        assert_eq!(regions[0].len(), 11);
    }

    #[test]
    fn headerless_group_named_group_is_not_dropped() {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "phase_tools_regions_group_{}.tsv",
            std::process::id()
        ));
        std::fs::write(&path, "group\tchr1\t10\t20\tcopy1\n").unwrap();
        let regions = read_regions_tsv(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(regions.len(), 1);
        assert_eq!(regions[0].group, "group");
        assert_eq!(regions[0].chrom, "chr1");
    }
}
