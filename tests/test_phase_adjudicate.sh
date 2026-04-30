#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_adjudicate_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/pairs.tsv" <<'EOF'
chrom	prev_pos	pos	truth_ps	query_ps	prev_orientation	orientation	status	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chr1	1	2	1	1	0	1	switch	0|1	0|1	0|1	1|0
chr1	2	4	1	1	0	1	switch	0|1	1|0	0|1	0|1
chr1	4	5	1	1	0	0	match	1|0	1|0	0|1	0|1
EOF

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/pairs.tsv" > "$tmp/out.tsv"

grep -qx $'chrom\tprev_pos\tpos\ttruth_parity\tquery_parity\tusable_reads\tspanning_reads\tinformative_reads\ttruth_support\tquery_support\tother_support\tmean_min_baseq\tmean_mapq\tmapq_unknown_reads\tforward_reads\treverse_reads\twinner\tambiguous\treason' "$tmp/out.tsv"
grep -qx $'chr1\t1\t2\t0\t1\t4\t4\t4\t4\t0\t0\t40.000\t60.000\t0\t4\t0\ttruth\tfalse\tevidence' "$tmp/out.tsv"
grep -qx $'chr1\t2\t4\t1\t0\t4\t4\t4\t4\t0\t0\t40.000\t60.000\t0\t4\t0\ttruth\tfalse\tevidence' "$tmp/out.tsv"
grep -qx $'chr1\t4\t5\t0\t0\t4\t4\t4\t4\t4\t0\t40.000\t60.000\t0\t4\t0\tboth\tfalse\tsame_phase' "$tmp/out.tsv"

"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/pairs.tsv" \
  --min-baseq 41 > "$tmp/baseq.tsv"
grep -qx $'chr1\t1\t2\t0\t1\t4\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tno_informative_reads' "$tmp/baseq.tsv"
grep -qx $'chr1\t4\t5\t0\t0\t4\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tboth\tfalse\tsame_phase' "$tmp/baseq.tsv"

cat > "$tmp/unsupported_pairs.tsv" <<'EOF'
chrom	prev_pos	pos	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chr1	1	2	0/1	0|1	0|1	0|1
chr1	1	2	0|1	0|1	0|0	0|1
EOF
"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/unsupported_pairs.tsv" > "$tmp/unsupported.tsv"
grep -qx $'chr1\t1\t2\tNA\tNA\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tunsupported_truth_gt' "$tmp/unsupported.tsv"
grep -qx $'chr1\t1\t2\t0\tNA\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tunsupported_query_gt' "$tmp/unsupported.tsv"

cat > "$tmp/bad_pairs.tsv" <<'EOF'
chrom	prev_pos	pos	prev_gt_truth	gt_truth	prev_gt_query	gt_query
chr1	1	3	0|1	0|1	0|1	0|1
EOF
"$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$fixtures/read_phase.vcf" \
  --pair-tsv "$tmp/bad_pairs.tsv" > "$tmp/missing.tsv"
grep -qx $'chr1\t1\t3\t0\t0\t0\t0\t0\t0\t0\t0\tNA\tNA\t0\t0\t0\tnone\ttrue\tmissing_variant' "$tmp/missing.tsv"

cat > "$tmp/duplicate_pos.vcf" <<'EOF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	A	G	.	PASS	.
chr1	1	.	A	C	.	PASS	.
chr1	2	.	C	T	.	PASS	.
EOF
if "$bin" \
  --reference "$fixtures/ref.fa" \
  --bam "$fixtures/read_phase.bam" \
  --variants "$tmp/duplicate_pos.vcf" \
  --pair-tsv "$tmp/pairs.tsv" > "$tmp/duplicate.out" 2> "$tmp/duplicate.err"; then
  echo "phase_adjudicate unexpectedly accepted duplicate-position SNVs" >&2
  exit 1
fi
grep -q 'duplicate biallelic SNV records at chr1:1' "$tmp/duplicate.err"

echo "phase_adjudicate tests passed"
