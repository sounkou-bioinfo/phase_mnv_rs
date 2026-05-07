use libc::c_void;
use rust_htslib::htslib;
use std::ffi::CString;

/// RAII wrapper around an htslib FASTA index.
pub struct Fai(*mut htslib::faidx_t);

impl Fai {
    /// Load a FASTA index from `path`.
    pub fn from_path(path: &str) -> Result<Self, String> {
        let c_path = CString::new(path.as_bytes()).map_err(|_| "reference path contains NUL")?;
        let fai = unsafe { htslib::fai_load(c_path.as_ptr()) };
        if fai.is_null() {
            Err(format!("cannot load FASTA index for '{path}'"))
        } else {
            Ok(Self(fai))
        }
    }

    /// Fetch a 1-based inclusive interval as uppercase ASCII bases.
    pub fn fetch_seq(&self, chrom: &str, start1: i64, end1: i64) -> Result<Vec<u8>, String> {
        if start1 < 1 || end1 < start1 {
            return Err(format!("invalid FASTA interval {chrom}:{start1}-{end1}"));
        }
        let c_chrom = CString::new(chrom.as_bytes()).map_err(|_| "contig contains NUL")?;
        let mut len: htslib::hts_pos_t = 0;
        let ptr = unsafe {
            htslib::faidx_fetch_seq64(
                self.0,
                c_chrom.as_ptr(),
                (start1 - 1) as htslib::hts_pos_t,
                (end1 - 1) as htslib::hts_pos_t,
                &mut len,
            )
        };
        if ptr.is_null() || len <= 0 {
            unsafe {
                if !ptr.is_null() {
                    libc::free(ptr as *mut c_void);
                }
            }
            return Err(format!(
                "failed to fetch reference interval {chrom}:{start1}-{end1}"
            ));
        }
        let seq = unsafe { std::slice::from_raw_parts(ptr as *const u8, len as usize) }
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect::<Vec<_>>();
        unsafe { libc::free(ptr as *mut c_void) };
        Ok(seq)
    }

    /// Fetch a single 1-based reference base as uppercase ASCII.
    pub fn fetch_base(&self, chrom: &str, pos1: i64) -> Result<u8, String> {
        if pos1 < 1 {
            return Err(format!("invalid FASTA position {chrom}:{pos1}"));
        }
        let seq = self.fetch_seq(chrom, pos1, pos1)?;
        if seq.len() != 1 {
            return Err(format!("failed to fetch reference base {chrom}:{pos1}"));
        }
        Ok(seq[0])
    }
}

impl Drop for Fai {
    fn drop(&mut self) {
        unsafe { htslib::fai_destroy(self.0) };
    }
}
