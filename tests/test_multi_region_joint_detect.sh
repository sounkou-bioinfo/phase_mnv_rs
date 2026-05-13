#!/usr/bin/env bash
set -euo pipefail

bin=${1:?usage: $0 <multi_region_joint_detect_binary>}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

if ! command -v samtools >/dev/null 2>&1; then
  echo "samtools not found; skipping multi_region_joint_detect regression test"
  exit 0
fi

python3 - <<'PY' > "$tmp/ref.fa"
print('>chr1')
print('A' * 80)
print('>chr2')
print('A' * 80)
PY
samtools faidx "$tmp/ref.fa"

cat > "$tmp/regions.tsv" <<'TSV'
group	chrom	start	end	copy
G1	chr1	10	20	copy1
G1	chr2	30	40	copy2
TSV

cat > "$tmp/reads.sam" <<'SAM'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:80
@SQ	SN:chr2	LN:80
c1_alt_1	0	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_alt_2	0	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_alt_3	0	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_ref_1	0	chr1	10	60	11M	*	0	0	AAAAAAAAAAA	IIIIIIIIIII
c1_ref_2	0	chr1	10	60	11M	*	0	0	AAAAAAAAAAA	IIIIIIIIIII
c1_dup_alt	1024	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_secondary_alt	256	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_supp_alt	2048	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_qcfail_alt	512	chr1	10	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c1_mapq255_alt	0	chr1	10	255	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c2_alt_1	0	chr2	30	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c2_alt_2	0	chr2	30	60	11M	*	0	0	AAAAACAAAAA	IIIIIIIIIII
c2_ref_1	0	chr2	30	60	11M	*	0	0	AAAAAAAAAAA	IIIIIIIIIII
c2_ref_2	0	chr2	30	60	11M	*	0	0	AAAAAAAAAAA	IIIIIIIIIII
SAM
samtools view -bS "$tmp/reads.sam" | samtools sort -o "$tmp/reads.bam"
samtools index "$tmp/reads.bam"

"$bin" \
  --reference "$tmp/ref.fa" \
  --regions "$tmp/regions.tsv" \
  --min-mapq 20 \
  --min-alt-count 2 \
  --min-alt-fraction 0.25 \
  "$tmp/reads.bam" > "$tmp/out.tsv"

cat > "$tmp/expected.tsv" <<'TSV'
group	offset1	alt	alt_positive_depth	alt_positive_alt_count	regions_with_alt	region_count	per_region
G1	6	C	9	5	2	2	copy1|chr1:15|A|5|3|0.600000;copy2|chr2:35|A|4|2|0.500000
TSV

diff -u "$tmp/expected.tsv" "$tmp/out.tsv"

"$bin" \
  --reference "$tmp/ref.fa" \
  --regions "$tmp/regions.tsv" \
  --min-mapq 20 \
  --min-alt-count 4 \
  "$tmp/reads.bam" > "$tmp/empty.tsv"
[[ $(awk 'END {print NR + 0}' "$tmp/empty.tsv") == 1 ]]

"$bin" \
  --reference "$tmp/ref.fa" \
  --regions "$tmp/regions.tsv" \
  --min-alt-count 4 \
  "$tmp/reads.bam" > "$tmp/mapq255.tsv"
cat > "$tmp/mapq255.expected.tsv" <<'TSV'
group	offset1	alt	alt_positive_depth	alt_positive_alt_count	regions_with_alt	region_count	per_region
G1	6	C	6	4	1	2	copy1|chr1:15|A|6|4|0.666667
TSV
diff -u "$tmp/mapq255.expected.tsv" "$tmp/mapq255.tsv"

cat > "$tmp/regions.headerless.tsv" <<'TSV'
group	chr1	10	20	copy1
group	chr2	30	40	copy2
TSV
"$bin" \
  --reference "$tmp/ref.fa" \
  --regions "$tmp/regions.headerless.tsv" \
  --min-mapq 20 \
  --min-alt-count 2 \
  "$tmp/reads.bam" > "$tmp/headerless.tsv"
grep -q $'^group\t6\tC\t9\t5\t2\t2\tcopy1' "$tmp/headerless.tsv"

cat > "$tmp/regions.out_of_bounds.tsv" <<'TSV'
group	chrom	start	end	copy
Gbad	chr1	75	90	bad
TSV
if "$bin" \
  --reference "$tmp/ref.fa" \
  --regions "$tmp/regions.out_of_bounds.tsv" \
  "$tmp/reads.bam" > "$tmp/out_of_bounds.tsv" 2> "$tmp/out_of_bounds.err"; then
  echo "out-of-bounds region unexpectedly succeeded" >&2
  exit 1
fi
grep -q 'returned length' "$tmp/out_of_bounds.err"

echo "multi_region_joint_detect tests passed"
