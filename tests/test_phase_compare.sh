#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <phase_compare_binary>}
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$repo_root/tests/fixtures"
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

"$bin" \
  --sample S1 \
  -r "$fixtures/ref.fa" \
  --engine xcmp \
  --no-roc \
  --no-decompose \
  --switch-error-bed "$tmp/switches.bed" \
  --tsv-pairwise "$tmp/pairs.tsv" \
  --summary-tsv "$tmp/summary.copy.tsv" \
  "$fixtures/phase_compare_truth.vcf" \
  "$fixtures/phase_compare_query.vcf" > "$tmp/summary.tsv"

cmp "$tmp/summary.tsv" "$tmp/summary.copy.tsv"

"$bin" \
  --sample S1 \
  -o "$tmp/prefix" \
  "$fixtures/phase_compare_truth.vcf" \
  "$fixtures/phase_compare_query.vcf" > "$tmp/prefix.stdout.tsv"
test -s "$tmp/prefix.summary.tsv"

header=$(head -1 "$tmp/summary.tsv")
[[ "$header" == chrom$'\t'truth_records* ]]

total=$(awk -F'\t' '$1=="TOTAL" {print}' "$tmp/summary.tsv")
[[ -n "$total" ]]

# Columns:
# 1 chrom, 2 truth_records, 3 query_records, 4 common_records,
# 7 common_het, 11 both_phased_het_with_ps, 12 intersection_blocks,
# 13 intersection_variants, 14 assessed_pairs, 15 phase_match_pairs,
# 16 switch_errors, 17 switch_rate, 18 blockwise_hamming, 19 hamming_rate.
# Singleton PS blocks are counted as phased sites but not as intersection blocks.
awk -F'\t' '$1=="TOTAL" {
  if ($2 != 5 || $3 != 5 || $4 != 5 || $7 != 5 || $11 != 4 ||
      $12 != 1 || $13 != 3 || $14 != 2 || $15 != 1 || $16 != 1 ||
      $17 != "0.500000" || $18 != 1 || $19 != "0.333333") exit 1
}' "$tmp/summary.tsv"

[[ $(awk 'END {print NR + 0}' "$tmp/switches.bed") == 1 ]]
grep -qx $'chr1\t1\t3' "$tmp/switches.bed"

[[ $(awk 'END {print NR + 0}' "$tmp/pairs.tsv") == 3 ]]
grep -q $'chr1\t2\t3\t1\t1\t0\t1\tswitch' "$tmp/pairs.tsv"

cat > "$tmp/hp_truth.vcf" <<'VCF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	G	.	PASS	.	GT:PS	0|1:1
chr1	2	.	C	T	.	PASS	.	GT:PS	0|1:1
VCF
cat > "$tmp/hp_query.vcf" <<'VCF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=HP,Number=.,Type=String,Description="Phasing haplotype identifier">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1	.	A	G	.	PASS	.	GT:HP	0/1:1-2,1-1
chr1	2	.	C	T	.	PASS	.	GT:HP	0/1:1-2,1-1
VCF
"$bin" --sample S1 "$tmp/hp_truth.vcf" "$tmp/hp_query.vcf" > "$tmp/hp.summary.tsv"
awk -F'\t' '$1=="TOTAL" {
  if ($11 != 2 || $12 != 1 || $13 != 2 || $14 != 1 || $15 != 1 || $16 != 0 || $17 != "0.000000") exit 1
}' "$tmp/hp.summary.tsv"

cat > "$tmp/indel_complex_truth.vcf" <<'VCF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	10	.	A	G	.	PASS	.	GT:PS	0|1:10
chr1	20	.	AT	A	.	PASS	.	GT:PS	0|1:10
chr1	30	.	C	CA	.	PASS	.	GT:PS	0|1:10
chr1	40	.	AG	TC	.	PASS	.	GT:PS	0|1:10
VCF
cat > "$tmp/indel_complex_query.vcf" <<'VCF'
##fileformat=VCFv4.3
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	10	.	A	G	.	PASS	.	GT:PS	0|1:10
chr1	20	.	AT	A	.	PASS	.	GT:PS	0|1:10
chr1	30	.	C	CA	.	PASS	.	GT:PS	1|0:10
chr1	40	.	AG	TC	.	PASS	.	GT:PS	1|0:10
VCF
"$bin" --sample S1 --pair-tsv "$tmp/indel_complex.pairs.tsv" \
  "$tmp/indel_complex_truth.vcf" "$tmp/indel_complex_query.vcf" > "$tmp/indel_complex.summary.tsv"
awk -F'\t' '$1=="TOTAL" {
  if ($2 != 4 || $3 != 4 || $4 != 4 || $7 != 4 || $11 != 4 ||
      $12 != 1 || $13 != 4 || $14 != 3 || $15 != 2 || $16 != 1 ||
      $17 != "0.333333" || $18 != 2 || $19 != "0.500000") exit 1
}' "$tmp/indel_complex.summary.tsv"
grep -q $'chr1\t20\t30\t10\t10\t0\t1\tswitch' "$tmp/indel_complex.pairs.tsv"

"$bin" --sample S1 --only-snvs "$tmp/indel_complex_truth.vcf" \
  "$tmp/indel_complex_query.vcf" > "$tmp/indel_complex.only_snvs.summary.tsv"
awk -F'\t' '$1=="TOTAL" {
  if ($2 != 1 || $3 != 1 || $4 != 1 || $7 != 1 || $11 != 1 ||
      $12 != 0 || $13 != 0 || $14 != 0 || $16 != 0 || $18 != 0) exit 1
}' "$tmp/indel_complex.only_snvs.summary.tsv"

mkdir -p "$tmp/repro_inputs"
cp "$fixtures/phase_compare_truth.vcf" "$tmp/repro_inputs/truth.vcf"
cp "$fixtures/phase_compare_query.vcf" "$tmp/repro_inputs/query.vcf"
cat > "$tmp/repro_inputs/manifest.tsv" <<'EOF'
method	run	vcf	sample
synthetic	run/a	truth.vcf	S1
synthetic	run a	query.vcf	S1
synthetic	run_a	truth.vcf	S1
EOF
python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro" \
  --sample S1 \
  --only-snvs \
  --write-pairs > "$tmp/repro.stdout"
grep -q 'wrote 3 pairwise comparison' "$tmp/repro.stdout"
[[ $(awk 'END {print NR + 0}' "$tmp/repro/pairwise_long.tsv") == 4 ]]
[[ $(find "$tmp/repro/pairwise" -name '*.summary.tsv' | wc -l | tr -d ' ') == 3 ]]
grep -q '^synthetic' "$tmp/repro/summary_by_method.tsv"
grep -q $'synthetic\t3\tblockwise_hamming_rate' "$tmp/repro/summary_by_method.tsv"
ls "$tmp/repro/pairs"/*.pairs.tsv >/dev/null

python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro" \
  --sample S1 \
  --limit-pairs 1 \
  --write-pairs > "$tmp/repro_rerun.stdout"
[[ $(find "$tmp/repro/pairwise" -name '*.summary.tsv' | wc -l | tr -d ' ') == 1 ]]
[[ $(find "$tmp/repro/pairs" -name '*.pairs.tsv' | wc -l | tr -d ' ') == 1 ]]
[[ $(awk 'END {print NR + 0}' "$tmp/repro/pairwise_long.tsv") == 2 ]]

cat > "$tmp/repro_mixed_samples.tsv" <<EOF
method	run	vcf	sample
synthetic	run1	$fixtures/phase_compare_truth.vcf	S1
synthetic	run2	$fixtures/phase_compare_query.vcf	S2
EOF
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/repro_mixed_samples.tsv" \
  --out-dir "$tmp/repro_mixed" > "$tmp/repro_mixed.stdout" 2> "$tmp/repro_mixed.stderr"; then
  echo "mixed-sample reproducibility manifest unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'multiple samples for method' "$tmp/repro_mixed.stderr"

if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$tmp/no_such_phase_compare" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro_missing_bin" > "$tmp/repro_missing_bin.stdout" 2> "$tmp/repro_missing_bin.stderr"; then
  echo "missing phase_compare unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'phase_compare not found' "$tmp/repro_missing_bin.stderr"
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$tmp/no_such_phase_compare" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro" > "$tmp/repro_missing_bin_existing.stdout" 2> "$tmp/repro_missing_bin_existing.stderr"; then
  echo "missing phase_compare into existing out-dir unexpectedly succeeded" >&2
  exit 1
fi
[[ $(find "$tmp/repro/pairwise" -name '*.summary.tsv' | wc -l | tr -d ' ') == 1 ]]
[[ $(awk 'END {print NR + 0}' "$tmp/repro/pairwise_long.tsv") == 2 ]]

if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/no_manifest.tsv" \
  --out-dir "$tmp/repro_missing_manifest" > "$tmp/repro_missing_manifest.stdout" 2> "$tmp/repro_missing_manifest.stderr"; then
  echo "missing manifest unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'manifest not found' "$tmp/repro_missing_manifest.stderr"

cat > "$tmp/repro_bad_fields.tsv" <<EOF
method	run	vcf	sample
synthetic	run1	$fixtures/phase_compare_truth.vcf	S1	extra
synthetic	run2	$fixtures/phase_compare_query.vcf	S1
EOF
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/repro_bad_fields.tsv" \
  --out-dir "$tmp/repro_bad_fields" > "$tmp/repro_bad_fields.stdout" 2> "$tmp/repro_bad_fields.stderr"; then
  echo "malformed manifest unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'more fields than the header' "$tmp/repro_bad_fields.stderr"
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/repro_bad_fields.tsv" \
  --out-dir "$tmp/repro" > "$tmp/repro_bad_fields_existing.stdout" 2> "$tmp/repro_bad_fields_existing.stderr"; then
  echo "malformed manifest into existing out-dir unexpectedly succeeded" >&2
  exit 1
fi
[[ $(find "$tmp/repro/pairwise" -name '*.summary.tsv' | wc -l | tr -d ' ') == 1 ]]
[[ $(awk 'END {print NR + 0}' "$tmp/repro/pairwise_long.tsv") == 2 ]]

if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$bin" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro_negative_limit" \
  --limit-pairs -1 > "$tmp/repro_negative_limit.stdout" 2> "$tmp/repro_negative_limit.stderr"; then
  echo "negative --limit-pairs unexpectedly succeeded" >&2
  exit 1
fi
grep -q -- '--limit-pairs must be >= 0' "$tmp/repro_negative_limit.stderr"

cat > "$tmp/nonexec_phase_compare" <<'SH'
#!/usr/bin/env bash
exit 0
SH
chmod -x "$tmp/nonexec_phase_compare"
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$tmp/nonexec_phase_compare" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro_nonexec" > "$tmp/repro_nonexec.stdout" 2> "$tmp/repro_nonexec.stderr"; then
  echo "non-executable phase_compare unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'cannot execute phase_compare' "$tmp/repro_nonexec.stderr"

cat > "$tmp/failing_phase_compare" <<'SH'
#!/usr/bin/env bash
echo 'intentional failure' >&2
exit 42
SH
chmod +x "$tmp/failing_phase_compare"
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$tmp/failing_phase_compare" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro_failing" > "$tmp/repro_failing.stdout" 2> "$tmp/repro_failing.stderr"; then
  echo "failing phase_compare unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'phase_compare failed' "$tmp/repro_failing.stderr"

cat > "$tmp/fake_phase_compare" <<'SH'
#!/usr/bin/env bash
printf 'not_a_phase_compare_summary\n'
SH
chmod +x "$tmp/fake_phase_compare"
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$tmp/fake_phase_compare" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro_malformed" > "$tmp/repro_malformed.stdout" 2> "$tmp/repro_malformed.stderr"; then
  echo "malformed phase_compare output unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'invalid phase_compare summary' "$tmp/repro_malformed.stderr"

cat > "$tmp/truncated_phase_compare" <<'SH'
#!/usr/bin/env bash
printf 'chrom\ttruth_records\tquery_records\nTOTAL\t1\n'
SH
chmod +x "$tmp/truncated_phase_compare"
if python3 "$repo_root/scripts/phase_reproducibility_matrix.py" \
  --phase-compare "$tmp/truncated_phase_compare" \
  --manifest "$tmp/repro_inputs/manifest.tsv" \
  --out-dir "$tmp/repro_truncated" > "$tmp/repro_truncated.stdout" 2> "$tmp/repro_truncated.stderr"; then
  echo "truncated phase_compare output unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'TOTAL row has' "$tmp/repro_truncated.stderr"

echo "phase_compare tests passed"
