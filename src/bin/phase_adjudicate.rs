use phase_tools::assembly::fermi_lite::{assemble_reads, AssembleOptions, AssemblyRead, Unitig};
use phase_tools::io::fasta::Fai;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::{self, Read as BcfRead};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

fn usage() -> &'static str {
    "usage: phase_adjudicate --reference ref.fa --bam reads.bam --variants variants.vcf|bcf --pair-tsv pairs.tsv [options]\n\n\
Experimental read-evidence adjudicator for phase_compare pair TSV rows. The\n\
first implementation is deliberately narrow: it adjudicates biallelic SNV pairs\n\
by counting reads spanning both sites and comparing observed read allele parity\n\
with the truth/query phased GT patterns from the pair TSV. No MAPQ or baseQ\n\
filter is applied by default; optional thresholds are explicit.\n\n\
options:\n\
  -r, --reference FILE       Required FASTA reference (used for CRAM decoding)\n\
      --bam FILE             Indexed BAM/CRAM read evidence\n\
      --variants FILE        VCF/BCF containing the pair sites and alleles\n\
      --pair-tsv FILE        phase_compare --pair-tsv output\n\
  -o, --output FILE          Output TSV (default: stdout)\n\
  -@, --threads N            htslib reader threads (default: 1)\n\
      --min-mapq N           Optional MAPQ cutoff (default: 0; no cutoff)\n\
      --min-baseq N          Optional per-site baseQ cutoff (default: 0; no cutoff)\n\
      --include-duplicates   Include duplicate reads\n\
      --include-secondary    Include secondary alignments\n\
      --include-supplementary Include supplementary alignments\n\
      --assembly-fasta FILE  Experimental fermi-lite local assembly sidecar FASTA\n\
      --assembly-tsv FILE    Experimental per-unitig assembly parity sidecar TSV\n\
      --use-assembly-decision Use assembly evidence to break otherwise ambiguous decisions\n\
      --assembly-window N    Bases of padding around each pair for assembly (default: 100)\n\
      --assembly-context N   Bases around pair used for unitig parity scoring (default: 10)\n\
      --assembly-max-reads N Maximum reads per pair assembly (default: 200)\n\
      --assembly-min-asm-ovlp N fermi-lite minimum assembly overlap (default: 21)\n\
  -h, --help                 Show this help\n"
}

#[derive(Debug)]
struct Config {
    reference: String,
    bam: String,
    variants: String,
    pair_tsv: String,
    output: Option<String>,
    threads: usize,
    min_mapq: u8,
    min_baseq: u8,
    include_duplicates: bool,
    include_secondary: bool,
    include_supplementary: bool,
    assembly_fasta: Option<String>,
    assembly_tsv: Option<String>,
    use_assembly_decision: bool,
    assembly_window: i64,
    assembly_context: i64,
    assembly_max_reads: usize,
    assembly_min_asm_overlap: i32,
}

#[derive(Debug, Clone)]
struct PairRecord {
    chrom: String,
    prev_pos: i64,
    pos: i64,
    prev_gt_truth: String,
    gt_truth: String,
    prev_gt_query: String,
    gt_query: String,
}

#[derive(Debug, Clone)]
struct Variant {
    ref_base: u8,
    alt_base: u8,
}

#[derive(Debug, Default)]
struct PairAssembly {
    input_reads: usize,
    unitigs: Vec<Unitig>,
}

#[derive(Debug, Default)]
struct AssemblyEvidence {
    informative_unitigs: u64,
    ambiguous_unitigs: u64,
    truth_support: u64,
    query_support: u64,
    other_support: u64,
}

#[derive(Debug)]
struct AssemblyCall {
    start1: i64,
    end1: i64,
    best_prev_allele: Option<u8>,
    best_allele: Option<u8>,
    best_parity: Option<u8>,
    best_distance: usize,
    second_distance: Option<usize>,
    status: &'static str,
}

#[derive(Debug, Default)]
struct PairEvidence {
    usable_reads: u64,
    spanning_reads: u64,
    informative_reads: u64,
    truth_support: u64,
    query_support: u64,
    other_support: u64,
    sum_min_baseq: u64,
    sum_mapq: u64,
    mapq_known_reads: u64,
    mapq_unknown_reads: u64,
    forward_reads: u64,
    reverse_reads: u64,
}

impl PairEvidence {
    fn mean_min_baseq(&self) -> String {
        if self.informative_reads == 0 {
            "NA".to_string()
        } else {
            format!(
                "{:.3}",
                self.sum_min_baseq as f64 / self.informative_reads as f64
            )
        }
    }

    fn mean_mapq(&self) -> String {
        if self.mapq_known_reads == 0 {
            "NA".to_string()
        } else {
            format!("{:.3}", self.sum_mapq as f64 / self.mapq_known_reads as f64)
        }
    }
}

fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    std::process::exit(1);
}

fn parse_i64(s: &str, name: &str) -> i64 {
    s.parse::<i64>()
        .unwrap_or_else(|_| die(&format!("{name} must be an integer")))
}

fn parse_u8(s: &str, name: &str) -> u8 {
    let value = parse_i64(s, name);
    if !(0..=255).contains(&value) {
        die(&format!("{name} must be between 0 and 255"));
    }
    value as u8
}

fn parse_args() -> Config {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{}", usage());
        std::process::exit(0);
    }

    let mut reference = None;
    let mut bam = None;
    let mut variants = None;
    let mut pair_tsv = None;
    let mut output = None;
    let mut threads = 1usize;
    let mut min_mapq = 0u8;
    let mut min_baseq = 0u8;
    let mut include_duplicates = false;
    let mut include_secondary = false;
    let mut include_supplementary = false;
    let mut assembly_fasta = None;
    let mut assembly_tsv = None;
    let mut use_assembly_decision = false;
    let mut assembly_window = 100i64;
    let mut assembly_context = 10i64;
    let mut assembly_max_reads = 200usize;
    let mut assembly_min_asm_overlap = 21i32;
    let mut positional = Vec::new();

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-r" | "--reference" => {
                i += 1;
                if i >= args.len() {
                    die("--reference requires an argument");
                }
                reference = Some(args[i].clone());
            }
            "--bam" => {
                i += 1;
                if i >= args.len() {
                    die("--bam requires an argument");
                }
                bam = Some(args[i].clone());
            }
            "--variants" | "--vcf" => {
                i += 1;
                if i >= args.len() {
                    die("--variants requires an argument");
                }
                variants = Some(args[i].clone());
            }
            "--pair-tsv" | "--pairs" => {
                i += 1;
                if i >= args.len() {
                    die("--pair-tsv requires an argument");
                }
                pair_tsv = Some(args[i].clone());
            }
            "-o" | "--output" => {
                i += 1;
                if i >= args.len() {
                    die("--output requires an argument");
                }
                output = Some(args[i].clone());
            }
            "-@" | "--threads" => {
                i += 1;
                if i >= args.len() {
                    die("--threads requires an argument");
                }
                threads = args[i]
                    .parse::<usize>()
                    .unwrap_or_else(|_| die("--threads must be an integer"));
                if threads == 0 {
                    die("--threads must be >= 1");
                }
            }
            "--min-mapq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-mapq requires an argument");
                }
                min_mapq = parse_u8(&args[i], "--min-mapq");
            }
            "--min-baseq" => {
                i += 1;
                if i >= args.len() {
                    die("--min-baseq requires an argument");
                }
                min_baseq = parse_u8(&args[i], "--min-baseq");
            }
            "--include-duplicates" => include_duplicates = true,
            "--include-secondary" => include_secondary = true,
            "--include-supplementary" => include_supplementary = true,
            "--assembly-fasta" => {
                i += 1;
                if i >= args.len() {
                    die("--assembly-fasta requires an argument");
                }
                assembly_fasta = Some(args[i].clone());
            }
            "--assembly-tsv" => {
                i += 1;
                if i >= args.len() {
                    die("--assembly-tsv requires an argument");
                }
                assembly_tsv = Some(args[i].clone());
            }
            "--use-assembly-decision" => use_assembly_decision = true,
            "--assembly-window" => {
                i += 1;
                if i >= args.len() {
                    die("--assembly-window requires an argument");
                }
                assembly_window = parse_i64(&args[i], "--assembly-window");
                if assembly_window < 0 {
                    die("--assembly-window must be >= 0");
                }
            }
            "--assembly-max-reads" => {
                i += 1;
                if i >= args.len() {
                    die("--assembly-max-reads requires an argument");
                }
                assembly_max_reads = args[i]
                    .parse::<usize>()
                    .unwrap_or_else(|_| die("--assembly-max-reads must be an integer"));
                if assembly_max_reads == 0 {
                    die("--assembly-max-reads must be >= 1");
                }
            }
            "--assembly-context" => {
                i += 1;
                if i >= args.len() {
                    die("--assembly-context requires an argument");
                }
                assembly_context = parse_i64(&args[i], "--assembly-context");
                if assembly_context < 0 {
                    die("--assembly-context must be >= 0");
                }
            }
            "--assembly-min-asm-ovlp" => {
                i += 1;
                if i >= args.len() {
                    die("--assembly-min-asm-ovlp requires an argument");
                }
                assembly_min_asm_overlap = args[i]
                    .parse::<i32>()
                    .unwrap_or_else(|_| die("--assembly-min-asm-ovlp must be an integer"));
                if assembly_min_asm_overlap < 1 {
                    die("--assembly-min-asm-ovlp must be >= 1");
                }
            }
            x if x.starts_with('-') => die(&format!("unknown option: {x}")),
            _ => positional.push(args[i].clone()),
        }
        i += 1;
    }

    if !positional.is_empty() {
        die("unexpected positional arguments; use --bam, --variants, and --pair-tsv");
    }
    if use_assembly_decision && assembly_tsv.is_none() {
        die("--use-assembly-decision requires --assembly-tsv so assembly-supported decisions are auditable");
    }

    Config {
        reference: reference.unwrap_or_else(|| die("--reference is required")),
        bam: bam.unwrap_or_else(|| die("--bam is required")),
        variants: variants.unwrap_or_else(|| die("--variants is required")),
        pair_tsv: pair_tsv.unwrap_or_else(|| die("--pair-tsv is required")),
        output,
        threads,
        min_mapq,
        min_baseq,
        include_duplicates,
        include_secondary,
        include_supplementary,
        assembly_fasta,
        assembly_tsv,
        use_assembly_decision,
        assembly_window,
        assembly_context,
        assembly_max_reads,
        assembly_min_asm_overlap,
    }
}

fn is_acgt(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn read_variants(path: &str, threads: usize) -> Result<HashMap<(String, i64), Variant>, String> {
    let mut reader =
        bcf::Reader::from_path(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    if threads > 1 {
        reader
            .set_threads(threads)
            .map_err(|e| format!("failed to set threads for {path}: {e}"))?;
    }
    let header = reader.header().clone();
    let mut variants = HashMap::new();
    for result in reader.records() {
        let record = result.map_err(|e| format!("failed to read {path}: {e}"))?;
        let Some(rid) = record.rid() else {
            continue;
        };
        let chrom = String::from_utf8_lossy(
            header
                .rid2name(rid)
                .map_err(|e| format!("failed to resolve RID {rid}: {e}"))?,
        )
        .into_owned();
        let alleles = record.alleles();
        if alleles.len() != 2 || alleles[0].len() != 1 || alleles[1].len() != 1 {
            continue;
        }
        let ref_base = alleles[0][0].to_ascii_uppercase();
        let alt_base = alleles[1][0].to_ascii_uppercase();
        if !is_acgt(ref_base) || !is_acgt(alt_base) || ref_base == alt_base {
            continue;
        }
        let pos = record.pos() + 1;
        if variants
            .insert((chrom.clone(), pos), Variant { ref_base, alt_base })
            .is_some()
        {
            return Err(format!(
                "duplicate biallelic SNV records at {chrom}:{pos}; phase_adjudicate pair TSV does not carry alleles for disambiguation"
            ));
        }
    }
    Ok(variants)
}

fn required_column(header: &[&str], name: &str) -> Result<usize, String> {
    header
        .iter()
        .position(|&col| col == name)
        .ok_or_else(|| format!("pair TSV is missing required column '{name}'"))
}

fn read_pairs(path: &str) -> Result<Vec<PairRecord>, String> {
    let file = File::open(path).map_err(|e| format!("cannot open pair TSV '{path}': {e}"))?;
    let mut lines = BufReader::new(file).lines();
    let header_line = lines
        .next()
        .ok_or_else(|| format!("pair TSV '{path}' is empty"))?
        .map_err(|e| format!("failed to read pair TSV header: {e}"))?;
    let header = header_line.split('\t').collect::<Vec<_>>();
    let chrom_col = required_column(&header, "chrom")?;
    let prev_pos_col = required_column(&header, "prev_pos")?;
    let pos_col = required_column(&header, "pos")?;
    let prev_gt_truth_col = required_column(&header, "prev_gt_truth")?;
    let gt_truth_col = required_column(&header, "gt_truth")?;
    let prev_gt_query_col = required_column(&header, "prev_gt_query")?;
    let gt_query_col = required_column(&header, "gt_query")?;

    let mut pairs = Vec::new();
    for (line_no, line) in lines.enumerate() {
        let line = line.map_err(|e| format!("failed to read pair TSV: {e}"))?;
        if line.trim().is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let needed = *[
            chrom_col,
            prev_pos_col,
            pos_col,
            prev_gt_truth_col,
            gt_truth_col,
            prev_gt_query_col,
            gt_query_col,
        ]
        .iter()
        .max()
        .expect("non-empty columns");
        if fields.len() <= needed {
            return Err(format!("pair TSV line {} has too few columns", line_no + 2));
        }
        pairs.push(PairRecord {
            chrom: fields[chrom_col].to_string(),
            prev_pos: fields[prev_pos_col]
                .parse::<i64>()
                .map_err(|_| format!("invalid prev_pos on pair TSV line {}", line_no + 2))?,
            pos: fields[pos_col]
                .parse::<i64>()
                .map_err(|_| format!("invalid pos on pair TSV line {}", line_no + 2))?,
            prev_gt_truth: fields[prev_gt_truth_col].to_string(),
            gt_truth: fields[gt_truth_col].to_string(),
            prev_gt_query: fields[prev_gt_query_col].to_string(),
            gt_query: fields[gt_query_col].to_string(),
        });
    }
    Ok(pairs)
}

fn phased_het_first_allele(gt: &str) -> Option<u8> {
    if gt.contains('/') {
        return None;
    }
    let parts = gt.split('|').collect::<Vec<_>>();
    if parts.len() != 2 {
        return None;
    }
    let a0 = match parts[0] {
        "0" => 0,
        "1" => 1,
        _ => return None,
    };
    let a1 = match parts[1] {
        "0" => 0,
        "1" => 1,
        _ => return None,
    };
    if a0 == a1 {
        None
    } else {
        Some(a0)
    }
}

fn pair_parity(prev_gt: &str, current_gt: &str) -> Option<u8> {
    Some(phased_het_first_allele(prev_gt)? ^ phased_het_first_allele(current_gt)?)
}

fn usable_record(record: &bam::Record, cfg: &Config) -> bool {
    !record.is_unmapped()
        && !record.is_quality_check_failed()
        && (cfg.min_mapq == 0 || (record.mapq() != 255 && record.mapq() >= cfg.min_mapq))
        && (cfg.include_duplicates || !record.is_duplicate())
        && (cfg.include_secondary || !record.is_secondary())
        && (cfg.include_supplementary || !record.is_supplementary())
}

fn base_at(record: &bam::Record, pos0: i64) -> Option<(u8, u8)> {
    let seq = record.seq();
    let qual = record.qual();
    let mut ref_pos = record.pos();
    let mut read_pos = 0usize;
    for cigar in record.cigar().iter() {
        match *cigar {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                for _ in 0..len {
                    if ref_pos == pos0 && read_pos < seq.len() {
                        return Some((
                            seq[read_pos].to_ascii_uppercase(),
                            qual.get(read_pos).copied().unwrap_or(255),
                        ));
                    }
                    ref_pos += 1;
                    read_pos += 1;
                }
            }
            Cigar::Ins(len) | Cigar::SoftClip(len) => read_pos += len as usize,
            Cigar::Del(len) | Cigar::RefSkip(len) => {
                if pos0 >= ref_pos && pos0 < ref_pos + len as i64 {
                    return None;
                }
                ref_pos += len as i64;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

fn record_sequence(record: &bam::Record) -> String {
    let seq = record.seq();
    (0..seq.len())
        .map(|i| seq[i].to_ascii_uppercase() as char)
        .collect()
}

fn assemble_pair(
    bam: &mut bam::IndexedReader,
    cfg: &Config,
    pair: &PairRecord,
) -> Result<PairAssembly, String> {
    let start0 = (pair.prev_pos - 1 - cfg.assembly_window).max(0);
    let end0 = pair.pos.max(pair.prev_pos) + cfg.assembly_window;
    bam.fetch((pair.chrom.as_bytes(), start0, end0))
        .map_err(|e| {
            format!(
                "failed to fetch assembly window {}:{}-{}: {e}",
                pair.chrom,
                start0 + 1,
                end0
            )
        })?;

    let mut reads = Vec::new();
    for result in bam.records() {
        let record = result.map_err(|e| format!("failed to read BAM/CRAM record: {e}"))?;
        if !usable_record(&record, cfg) {
            continue;
        }
        let seq = record_sequence(&record);
        if !seq.is_empty() {
            reads.push(AssemblyRead::sequence(seq));
            if reads.len() >= cfg.assembly_max_reads {
                break;
            }
        }
    }

    let mut options = AssembleOptions::default();
    options.threads = cfg.threads as i32;
    options.min_asm_overlap = cfg.assembly_min_asm_overlap;
    let unitigs = assemble_reads(&reads, &options)?;
    Ok(PairAssembly {
        input_reads: reads.len(),
        unitigs,
    })
}

fn write_assembly_fasta(
    pair: &PairRecord,
    assembly: &PairAssembly,
    writer: &mut Option<BufWriter<File>>,
) -> Result<(), String> {
    let Some(writer) = writer.as_mut() else {
        return Ok(());
    };
    for (idx, unitig) in assembly.unitigs.iter().enumerate() {
        writeln!(
            writer,
            ">{}:{}-{}|unitig={}|len={}|supporting_reads={}|input_reads={}\n{}",
            pair.chrom,
            pair.prev_pos,
            pair.pos,
            idx + 1,
            unitig.len,
            unitig.supporting_reads,
            assembly.input_reads,
            unitig.seq
        )
        .map_err(|e| format!("failed to write assembly FASTA: {e}"))?;
    }
    Ok(())
}

fn allele_base(variant: &Variant, allele: u8) -> u8 {
    if allele == 0 {
        variant.ref_base
    } else {
        variant.alt_base
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|base| match base.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut prev = (0..=b.len()).collect::<Vec<_>>();
    let mut curr = vec![0usize; b.len() + 1];
    for (i, ca) in a.iter().enumerate() {
        curr[0] = i + 1;
        for (j, cb) in b.iter().enumerate() {
            let sub_cost = usize::from(ca != cb);
            curr[j + 1] = (prev[j + 1] + 1).min(curr[j] + 1).min(prev[j] + sub_cost);
        }
        std::mem::swap(&mut prev, &mut curr);
    }
    prev[b.len()]
}

fn best_window_edit_distance(query: &[u8], haplotype: &[u8]) -> usize {
    if query.is_empty() || haplotype.is_empty() {
        return usize::MAX / 4;
    }
    if query.len() >= haplotype.len() {
        query
            .windows(haplotype.len())
            .map(|window| edit_distance(window, haplotype))
            .min()
            .unwrap_or(usize::MAX / 4)
    } else {
        edit_distance(query, haplotype)
    }
}

fn haplotype_for_pair(
    fai: &Fai,
    cfg: &Config,
    pair: &PairRecord,
    prev_variant: &Variant,
    variant: &Variant,
    prev_allele: u8,
    allele: u8,
) -> Result<(i64, i64, Vec<u8>), String> {
    let left = pair.prev_pos.min(pair.pos);
    let right = pair.prev_pos.max(pair.pos);
    let start1 = (left - cfg.assembly_context).max(1);
    let end1 = right + cfg.assembly_context;
    let mut seq = fai.fetch_seq(&pair.chrom, start1, end1)?;
    let prev_offset = (pair.prev_pos - start1) as usize;
    let offset = (pair.pos - start1) as usize;
    if prev_offset >= seq.len() || offset >= seq.len() {
        return Err(format!(
            "reference interval {}:{}-{} does not cover pair {}-{}",
            pair.chrom, start1, end1, pair.prev_pos, pair.pos
        ));
    }
    seq[prev_offset] = allele_base(prev_variant, prev_allele);
    seq[offset] = allele_base(variant, allele);
    Ok((start1, end1, seq))
}

fn assembly_call_for_unitig(
    fai: &Fai,
    cfg: &Config,
    pair: &PairRecord,
    prev_variant: &Variant,
    variant: &Variant,
    unitig: &Unitig,
) -> Result<AssemblyCall, String> {
    let unitig_seq = unitig.seq.as_bytes();
    let unitig_rc = reverse_complement(unitig_seq);
    let mut scored = Vec::<(usize, u8, u8)>::new();
    let mut interval = None;
    for prev_allele in 0..=1u8 {
        for allele in 0..=1u8 {
            let (start1, end1, haplotype) =
                haplotype_for_pair(fai, cfg, pair, prev_variant, variant, prev_allele, allele)?;
            interval = Some((start1, end1));
            let forward = best_window_edit_distance(unitig_seq, &haplotype);
            let reverse = best_window_edit_distance(&unitig_rc, &haplotype);
            scored.push((forward.min(reverse), prev_allele, allele));
        }
    }
    scored.sort_by_key(|(distance, prev_allele, allele)| (*distance, *prev_allele, *allele));
    let (start1, end1) = interval.expect("four haplotypes were scored");
    let best_distance = scored[0].0;
    let second_distance = scored.get(1).map(|entry| entry.0);
    let best = scored
        .iter()
        .filter(|(distance, _, _)| *distance == best_distance)
        .copied()
        .collect::<Vec<_>>();
    let mut best_parities = best
        .iter()
        .map(|(_, prev_allele, allele)| prev_allele ^ allele)
        .collect::<Vec<_>>();
    best_parities.sort_unstable();
    best_parities.dedup();

    let (best_prev_allele, best_allele) = if best.len() == 1 {
        (Some(best[0].1), Some(best[0].2))
    } else {
        (None, None)
    };
    let best_parity = if best_parities.len() == 1 {
        Some(best_parities[0])
    } else {
        None
    };
    let status = if best_parity.is_some() {
        "informative"
    } else {
        "ambiguous"
    };
    Ok(AssemblyCall {
        start1,
        end1,
        best_prev_allele,
        best_allele,
        best_parity,
        best_distance,
        second_distance,
        status,
    })
}

fn format_option_u8(value: Option<u8>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| "NA".to_string())
}

fn format_option_usize(value: Option<usize>) -> String {
    value
        .map(|v| v.to_string())
        .unwrap_or_else(|| "NA".to_string())
}

fn write_assembly_tsv(
    fai: &Fai,
    cfg: &Config,
    pair: &PairRecord,
    prev_variant: &Variant,
    variant: &Variant,
    truth_parity: u8,
    query_parity: u8,
    assembly: &PairAssembly,
    writer: &mut Option<BufWriter<File>>,
) -> Result<AssemblyEvidence, String> {
    let mut evidence = AssemblyEvidence::default();
    for (idx, unitig) in assembly.unitigs.iter().enumerate() {
        let call = assembly_call_for_unitig(fai, cfg, pair, prev_variant, variant, unitig)?;
        let supports_truth = call.best_parity == Some(truth_parity);
        let supports_query = call.best_parity == Some(query_parity);
        if call.best_parity.is_some() {
            evidence.informative_unitigs += 1;
            if supports_truth {
                evidence.truth_support += 1;
            }
            if supports_query {
                evidence.query_support += 1;
            }
            if !supports_truth && !supports_query {
                evidence.other_support += 1;
            }
        } else {
            evidence.ambiguous_unitigs += 1;
        }
        if let Some(writer) = writer.as_mut() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                pair.chrom,
                pair.prev_pos,
                pair.pos,
                idx + 1,
                assembly.input_reads,
                unitig.len,
                unitig.supporting_reads,
                call.start1,
                call.end1,
                format_option_u8(call.best_prev_allele),
                format_option_u8(call.best_allele),
                format_option_u8(call.best_parity),
                call.best_distance,
                format_option_usize(call.second_distance),
                call.status,
                supports_truth,
                supports_query
            )
            .map_err(|e| format!("failed to write assembly TSV: {e}"))?;
        }
    }
    Ok(evidence)
}

fn allele_for_base(base: u8, variant: &Variant) -> Option<u8> {
    let base = base.to_ascii_uppercase();
    if base == variant.ref_base {
        Some(0)
    } else if base == variant.alt_base {
        Some(1)
    } else {
        None
    }
}

fn adjudicate_pair(
    bam: &mut bam::IndexedReader,
    cfg: &Config,
    pair: &PairRecord,
    prev_variant: &Variant,
    variant: &Variant,
    truth_parity: u8,
    query_parity: u8,
) -> Result<PairEvidence, String> {
    let start0 = (pair.prev_pos - 1).max(0);
    let end0 = pair.pos.max(pair.prev_pos);
    bam.fetch((pair.chrom.as_bytes(), start0, end0))
        .map_err(|e| {
            format!(
                "failed to fetch {}:{}-{}: {e}",
                pair.chrom,
                start0 + 1,
                end0
            )
        })?;

    let mut evidence = PairEvidence::default();
    for result in bam.records() {
        let record = result.map_err(|e| format!("failed to read BAM/CRAM record: {e}"))?;
        if !usable_record(&record, cfg) {
            continue;
        }
        evidence.usable_reads += 1;
        let Some((prev_base, prev_q)) = base_at(&record, pair.prev_pos - 1) else {
            continue;
        };
        let Some((base, q)) = base_at(&record, pair.pos - 1) else {
            continue;
        };
        if prev_q < cfg.min_baseq || q < cfg.min_baseq {
            continue;
        }
        evidence.spanning_reads += 1;
        let Some(prev_allele) = allele_for_base(prev_base, prev_variant) else {
            evidence.other_support += 1;
            continue;
        };
        let Some(allele) = allele_for_base(base, variant) else {
            evidence.other_support += 1;
            continue;
        };
        evidence.informative_reads += 1;
        evidence.sum_min_baseq += u64::from(prev_q.min(q));
        if record.mapq() == 255 {
            evidence.mapq_unknown_reads += 1;
        } else {
            evidence.mapq_known_reads += 1;
            evidence.sum_mapq += u64::from(record.mapq());
        }
        if record.is_reverse() {
            evidence.reverse_reads += 1;
        } else {
            evidence.forward_reads += 1;
        }
        let observed_parity = prev_allele ^ allele;
        if observed_parity == truth_parity {
            evidence.truth_support += 1;
        }
        if observed_parity == query_parity {
            evidence.query_support += 1;
        }
        if observed_parity != truth_parity && observed_parity != query_parity {
            evidence.other_support += 1;
        }
    }
    Ok(evidence)
}

fn decision(
    evidence: &PairEvidence,
    truth_parity: u8,
    query_parity: u8,
) -> (&'static str, bool, &'static str) {
    if truth_parity == query_parity {
        return ("both", false, "same_phase");
    }
    if evidence.informative_reads == 0 {
        return ("none", true, "no_informative_reads");
    }
    if evidence.truth_support > evidence.query_support {
        ("truth", false, "evidence")
    } else if evidence.query_support > evidence.truth_support {
        ("query", false, "evidence")
    } else {
        ("tie", true, "equal_support")
    }
}

fn assembly_decision(
    evidence: &AssemblyEvidence,
    truth_parity: u8,
    query_parity: u8,
) -> Option<(&'static str, bool, &'static str)> {
    if truth_parity == query_parity || evidence.informative_unitigs == 0 {
        return None;
    }
    if evidence.truth_support > evidence.query_support {
        Some(("truth", false, "assembly_evidence"))
    } else if evidence.query_support > evidence.truth_support {
        Some(("query", false, "assembly_evidence"))
    } else if evidence.truth_support > 0 {
        Some(("tie", true, "assembly_equal_support"))
    } else {
        None
    }
}

fn open_output(path: &Option<String>) -> Result<Box<dyn Write>, String> {
    if let Some(path) = path {
        let file = File::create(path).map_err(|e| format!("cannot create output '{path}': {e}"))?;
        Ok(Box::new(BufWriter::new(file)))
    } else {
        Ok(Box::new(BufWriter::new(std::io::stdout())))
    }
}

fn run() -> Result<(), String> {
    let cfg = parse_args();
    let variants = read_variants(&cfg.variants, cfg.threads)?;
    let pairs = read_pairs(&cfg.pair_tsv)?;

    let mut bam = bam::IndexedReader::from_path(&cfg.bam)
        .map_err(|e| format!("cannot open indexed BAM/CRAM '{}': {e}", cfg.bam))?;
    bam.set_reference(&cfg.reference)
        .map_err(|e| format!("failed to set CRAM reference '{}': {e}", cfg.reference))?;
    if cfg.threads > 1 {
        bam.set_threads(cfg.threads)
            .map_err(|e| format!("failed to enable BAM/CRAM threads for '{}': {e}", cfg.bam))?;
    }

    let mut out = open_output(&cfg.output)?;
    let mut assembly_out = match cfg.assembly_fasta.as_deref() {
        Some(path) => {
            Some(BufWriter::new(File::create(path).map_err(|e| {
                format!("cannot create assembly FASTA '{path}': {e}")
            })?))
        }
        None => None,
    };
    let mut assembly_tsv = match cfg.assembly_tsv.as_deref() {
        Some(path) => {
            Some(BufWriter::new(File::create(path).map_err(|e| {
                format!("cannot create assembly TSV '{path}': {e}")
            })?))
        }
        None => None,
    };
    let score_assembly = cfg.assembly_tsv.is_some() || cfg.use_assembly_decision;
    let fai = if score_assembly {
        Some(Fai::from_path(&cfg.reference)?)
    } else {
        None
    };
    if let Some(writer) = assembly_tsv.as_mut() {
        writeln!(
            writer,
            "chrom\tprev_pos\tpos\tunitig\tinput_reads\tunitig_len\tsupporting_reads\tassembly_start\tassembly_end\tbest_prev_allele\tbest_allele\tbest_parity\tbest_distance\tsecond_distance\tstatus\tsupports_truth\tsupports_query"
        )
        .map_err(|e| format!("failed to write assembly TSV header: {e}"))?;
    }
    let write_assembly = assembly_out.is_some() || score_assembly;
    writeln!(
        out,
        "chrom\tprev_pos\tpos\ttruth_parity\tquery_parity\tusable_reads\tspanning_reads\tinformative_reads\ttruth_support\tquery_support\tother_support\tmean_min_baseq\tmean_mapq\tmapq_unknown_reads\tforward_reads\treverse_reads\twinner\tambiguous\treason"
    )
    .map_err(|e| format!("failed to write output header: {e}"))?;

    for pair in &pairs {
        let truth_parity = pair_parity(&pair.prev_gt_truth, &pair.gt_truth);
        let query_parity = pair_parity(&pair.prev_gt_query, &pair.gt_query);
        let Some(truth_parity) = truth_parity else {
            writeln!(
                out,
                "{}\t{}\t{}\tNA\tNA\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tunsupported_truth_gt",
                pair.chrom, pair.prev_pos, pair.pos
            )
            .map_err(|e| format!("failed to write output: {e}"))?;
            continue;
        };
        let Some(query_parity) = query_parity else {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\tNA\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tunsupported_query_gt",
                pair.chrom, pair.prev_pos, pair.pos, truth_parity
            )
            .map_err(|e| format!("failed to write output: {e}"))?;
            continue;
        };
        let prev_key = (pair.chrom.clone(), pair.prev_pos);
        let key = (pair.chrom.clone(), pair.pos);
        let Some(prev_variant) = variants.get(&prev_key) else {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tmissing_prev_variant",
                pair.chrom, pair.prev_pos, pair.pos, truth_parity, query_parity
            )
            .map_err(|e| format!("failed to write output: {e}"))?;
            continue;
        };
        let Some(variant) = variants.get(&key) else {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tmissing_variant",
                pair.chrom, pair.prev_pos, pair.pos, truth_parity, query_parity
            )
            .map_err(|e| format!("failed to write output: {e}"))?;
            continue;
        };
        let evidence = adjudicate_pair(
            &mut bam,
            &cfg,
            pair,
            prev_variant,
            variant,
            truth_parity,
            query_parity,
        )?;
        let mut assembly_evidence = AssemblyEvidence::default();
        if write_assembly {
            let assembly = assemble_pair(&mut bam, &cfg, pair)?;
            write_assembly_fasta(pair, &assembly, &mut assembly_out)?;
            if let Some(fai) = fai.as_ref() {
                assembly_evidence = write_assembly_tsv(
                    fai,
                    &cfg,
                    pair,
                    prev_variant,
                    variant,
                    truth_parity,
                    query_parity,
                    &assembly,
                    &mut assembly_tsv,
                )?;
            }
        }
        let (mut winner, mut ambiguous, mut reason) =
            decision(&evidence, truth_parity, query_parity);
        if cfg.use_assembly_decision && ambiguous {
            if let Some(assembly_call) =
                assembly_decision(&assembly_evidence, truth_parity, query_parity)
            {
                (winner, ambiguous, reason) = assembly_call;
            }
        }
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            pair.chrom,
            pair.prev_pos,
            pair.pos,
            truth_parity,
            query_parity,
            evidence.usable_reads,
            evidence.spanning_reads,
            evidence.informative_reads,
            evidence.truth_support,
            evidence.query_support,
            evidence.other_support,
            evidence.mean_min_baseq(),
            evidence.mean_mapq(),
            evidence.mapq_unknown_reads,
            evidence.forward_reads,
            evidence.reverse_reads,
            winner,
            ambiguous,
            reason
        )
        .map_err(|e| format!("failed to write output: {e}"))?;
    }
    out.flush()
        .map_err(|e| format!("failed to flush output: {e}"))?;
    if let Some(writer) = assembly_out.as_mut() {
        writer
            .flush()
            .map_err(|e| format!("failed to flush assembly FASTA: {e}"))?;
    }
    if let Some(writer) = assembly_tsv.as_mut() {
        writer
            .flush()
            .map_err(|e| format!("failed to flush assembly TSV: {e}"))?;
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        die(&e);
    }
}
