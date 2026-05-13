# Validation notes

This Rust-only repository validates `phase_tools-rs` library-backed binaries
with explicit tracked fixtures. No private/local paths are embedded in tests or
documentation.

## Positive behavior fixtures

Run:

```bash
make test
```

CI installs `samtools` and `bcftools` before fixture tests so generated-BAM,
compressed-output, and normalization regressions run in the hosted build. Local
runs without those tools skip only the tool-dependent checks and still exercise
tracked BAM/VCF fixtures.

Behavior fixtures cover:

- adjacent phased SNVs that become `TYPE=MNV`
- `--max-gap` behavior
- mixed SNV/indel blocks that become `TYPE=COMPLEX`
- selected ALT handling at multi-allelic sites
- symbolic/non-DNA ALT skipping semantics
- `--warn-on-n` warnings while preserving `N` as plain DNA
- native `unphase_vcf` conversion of phased GT separators to unphased GT
  while dropping phase-specific FORMAT tags by default on VCF, stdin VCF, BCF,
  and BGZF VCF output paths
- experimental Rust `--phase-from-bam` read-backed phasing on a tiny tracked
  BAM/BAI fixture before MNV construction, using the default exact
  single-sample MEC dynamic-programming algorithm; additional generated-BAM
  phasing regressions run when `samtools` is available
- Rust output format inference for plain VCF, BGZF-compressed VCF, and BCF,
  including `bcftools index` checks for constructed MNV VCF.GZ/VCF.BGZ/BCF and
  BAM-phased all-sites VCF.GZ output when `bcftools` is available
- `bcftools norm -f ... -c x` checks that emitted MNV/COMPLEX records are not
  further realigned or mismatch-removed on tracked fixtures
- Rust `--emit all-sites` header preservation: the original VCF header is kept
  and `phase_mnv` metadata is appended while BAM-backed `GT:PS` updates are
  applied to all input records in the tiny tracked fixture
- Rust `--mnv-algorithm nirvana-codon` same-codon SNV recomposition on a tracked
  BED-like codon-map fixture
- native `phase_compare` switch-error, phase-match, and blockwise-Hamming stats
  on a tiny tracked truth/query fixture
- bindgen-backed fermi-lite FFI assembly coverage through
  `fermi_lite_assemble`, including FASTQ/base-quality passthrough
- empirical BAM error-model summary, exact-Q event composition TSV,
  mapchk-like high-nonref site guard, and per-read-position TSV coverage through
  `bam_error_model` on tracked BAM/SAM fixtures, with no MAPQ filter by default
- experimental `phase_adjudicate` pair-level read-evidence coverage on tracked
  VCF/BAM fixtures, including same-phase rows, switched truth/query parity, and
  explicit baseQ filtering; generated-BAM fermi-lite assembly FASTA/TSV sidecar
  and guarded `--use-assembly-decision` fallback checks run when `samtools` is
  available
- experimental `bam_contamination` anchor-site contamination probe coverage on
  tracked BAM fixtures, including homozygous-alt reference infiltration,
  optional CHARR-like allele-frequency adjustment, explicit baseQ filtering, and
  unsupported genotype rejection
- experimental `bam_ancestry` Summix-style ancestry mixture probe coverage on
  tracked BAM fixtures, including population-column ordering, FASTA REF
  validation, explicit baseQ filtering, and no-observation failure behavior

## Negative/failure-mode fixtures

Run:

```bash
make negative-test
```

Negative fixtures and generated checks cover:

- missing required reference option
- missing input VCF/BCF
- missing FASTA reference
- unknown sample
- invalid negative `--max-gap`
- REF/FASTA mismatch
- `--unsupported-alleles fail` on selected unsupported ALT alleles
- truncated gzipped VCF input

The truncated-input fixture is tracked explicitly:

```text
tests/fixtures/truncated.vcf.gz
```

## WhatsHap + native `phase_compare`

The external conformance comparison is WhatsHap-based and no longer uses
hap.py. It uses the in-repository `phase_compare` binary, which is a narrow,
fast phase-concordance comparator.

Run:

```bash
make compare-whatshap-phase
```

The script compares two paths on the tracked tiny BAM/VCF fixture by default:

1. unphase the input VCF with the native `unphase_vcf` binary;
2. run external `whatshap phase` on the unphased VCF and BAM to create the truth
   all-sites phased VCF;
3. run Rust `phase_mnv_rs --emit all-sites --phase-from-bam` directly on the
   input VCF and BAM to create the query all-sites phased VCF;
4. run `phase_compare` on truth/query all-sites VCFs.

Default fixture inputs:

```text
tests/fixtures/read_phase.vcf
tests/fixtures/read_phase.bam
tests/fixtures/ref.fa
```

`phase_compare` reports exact shared variant records, common heterozygous sites,
phased sites with PS in both files, intersection PS blocks, assessed adjacent
pairs, switch errors, switch rate, blockwise Hamming distance, and blockwise
Hamming rate.

Important limitation: `phase_compare` is not a generic hap.py replacement. It is
for exact-site phasing/block concordance after both paths have been normalized to
the same input records. It does not perform variant representation matching,
ROC/stratification, decompose/atomize, or truth-query callset scoring.

The comparison script accepts thresholds by environment variable:

```bash
MAX_SWITCH_ERRORS=0 MAX_SWITCH_RATE=0 make compare-whatshap-phase
```

For exploratory local runs where non-perfect concordance is expected:

```bash
ALLOW_NONPERFECT=1 KEEP_TMP=1 make compare-whatshap-phase
```

The script sanitizes generated VCF headers before comparison to remove command
lines and local path-bearing records.

## Current phasing-quality benchmark frame

Tracked CI quality gates are deliberately small and deterministic. They require
zero switch errors on the tiny WhatsHap comparison fixture and exercise both the
MEC and greedy Rust phasing paths on synthetic BAM examples.

Larger private replicate checks currently use three complementary metrics:

1. within-run concordance against WhatsHap-phased output;
2. pairwise replicate reproducibility across independent runs;
3. switched-pair read/assembly adjudication for rows emitted by
   `phase_compare --pair-tsv`.

Recent private 13-run WES replicate summaries on two chromosomes were:

| chromosome | comparison | method | assessed pairs | switch errors | switch rate |
| --- | --- | --- | ---: | ---: | ---: |
| 22 | within-run vs WhatsHap | `rust_greedy` | 2,639 | 21 | 0.007958 |
| 22 | within-run vs WhatsHap | `rust_mec` | 2,639 | 25 | 0.009473 |
| 1 | within-run vs WhatsHap | `rust_greedy` | 8,364 | 138 | 0.016499 |
| 1 | within-run vs WhatsHap | `rust_mec` | 8,367 | 158 | 0.018884 |
| 22 | pairwise replicate | WhatsHap | 13,225 | 77 | 0.005822 |
| 22 | pairwise replicate | `rust_greedy` | 13,301 | 60 | 0.004511 |
| 22 | pairwise replicate | `rust_mec` | 13,301 | 73 | 0.005488 |
| 1 | pairwise replicate | WhatsHap | 42,395 | 899 | 0.021205 |
| 1 | pairwise replicate | `rust_greedy` | 42,472 | 680 | 0.016011 |
| 1 | pairwise replicate | `rust_mec` | 42,493 | 950 | 0.022357 |

Interpretation: `rust_greedy` is currently the strongest reproducibility
baseline on these private chromosomes, while `rust_mec` is closer to a
WhatsHap-style exact objective and remains the target for read-selection and
weighting improvements. These private numbers guide tuning; they are not a
public benchmark claim.

## Reproducibility matrix harness

`scripts/phase_reproducibility_matrix.py` runs `phase_compare` over all
within-method pairs in a manifest and writes both per-pair summaries and
aggregate method-level metrics. The input manifest is a TSV with at least
`method`, `run`, and `vcf` columns plus an optional `sample` column:

```text
method	run	vcf	sample
whatshap	runA	runA.whatshap.vcf.gz	S1
rust_greedy	runA	runA.rust_greedy.vcf.gz	S1
rust_greedy	runB	runB.rust_greedy.vcf.gz	S1
```

Example invocation:

```bash
python3 scripts/phase_reproducibility_matrix.py \
  --phase-compare target/release/phase_compare \
  --manifest local_runs/repro/manifest.tsv \
  --out-dir local_runs/repro/matrix \
  --sample S1 \
  --only-snvs \
  --write-pairs
```

Relative `vcf` paths are resolved against the manifest directory. The harness
writes `pairwise_long.tsv`, `summary_by_method.tsv`, individual
`pairwise/*.summary.tsv` files, and optional `pairs/*.pairs.tsv` files for
adjudication follow-up. It is intentionally data-agnostic and should be run with
privacy-safe paths under ignored local directories.

## Local private replicate checks

Larger non-committed checks are kept under ignored local output directories. One
recent local check used a private 13-run WES replicate panel on chromosomes 1 and
22. For each run, WhatsHap-phased VCFs were compared against Rust BAM-phased
outputs with `phase_compare --only-snvs --pair-tsv`; switch rows were then sent
to `phase_adjudicate` with `--assembly-fasta`, `--assembly-tsv`, and guarded
`--use-assembly-decision`.

Aggregate switched-pair adjudication results from that local run were:

| chromosome | method | switch pairs | read-evidence decisions | assembly-rescued decisions | no-informative-read rows |
| --- | --- | ---: | ---: | ---: | ---: |
| 1 | `rust_mec` | 118 | 101 | 7 | 9 |
| 1 | `rust_greedy` | 102 | 88 | 6 | 7 |
| 22 | `rust_mec` | 9 | 5 | 0 | 4 |
| 22 | `rust_greedy` | 8 | 4 | 0 | 4 |

Assembly fallback triggered only on chromosome 1 in this run: 13 total rescued
ambiguous decisions, with 9 supporting the WhatsHap-oriented parity and 4
supporting the Rust-oriented parity. These private validation artifacts are not
tracked and are intended as tuning evidence, not as a public benchmark claim.

## Local/private data policy

Default tests must use tracked fixtures only. Larger validation runs should be
launched with environment overrides, for example:

```bash
WHATSHAP_BIN=whatshap make compare-whatshap-phase
WHATSHAP_ENV=my-whatshap-env make compare-whatshap-phase
REF=ref.fa VCF=input.vcf.gz BAM=reads.bam SAMPLE=S1 ALLOW_NONPERFECT=1 make compare-whatshap-phase
```

Do not commit private paths, sample names, references, BAMs, or generated local
outputs. Use ignored directories such as `local_runs/` or `resources/` for local
validation runs.
