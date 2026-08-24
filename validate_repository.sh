#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "$0")" && pwd)
cd "$repo_dir"

fail() {
  printf 'ERROR: %s\n' "$1" >&2
  exit 1
}

test -f data/stacks_exports/populations.snps.vcf.gz || fail "VCF is missing"
gzip -t data/stacks_exports/populations.snps.vcf.gz || fail "VCF gzip test failed"

vcf_samples=$(gzip -dc data/stacks_exports/populations.snps.vcf.gz | awk -F '\t' '/^#CHROM/{n=NF-9} END{print n+0}')
vcf_variants=$(gzip -dc data/stacks_exports/populations.snps.vcf.gz | awk '!/^#/{n++} END{print n+0}')
popmap_samples=$(awk 'NF>=2{n++} END{print n+0}' data/metadata/popmap_excluding_MVZ216210.tsv)
coord_samples=$(awk 'NR>1 && NF>=3{n++} END{print n+0}' data/metadata/sample_geographic_coordinates.tsv)
unlinked_samples=$(wc -l < data/unlinked_snps/unlinked.fam | tr -d ' ')
unlinked_snps=$(wc -l < data/unlinked_snps/unlinked.bim | tr -d ' ')

test "$vcf_samples" -eq 99 || fail "expected 99 VCF samples; found $vcf_samples"
test "$vcf_variants" -eq 23482 || fail "expected 23,482 VCF variants; found $vcf_variants"
test "$popmap_samples" -eq 99 || fail "expected 99 popmap samples; found $popmap_samples"
test "$coord_samples" -eq 98 || fail "expected 98 coordinate samples; found $coord_samples"
test "$unlinked_samples" -eq 99 || fail "expected 99 unlinked samples; found $unlinked_samples"
test "$unlinked_snps" -eq 2557 || fail "expected 2,557 unlinked SNPs; found $unlinked_snps"

if gzip -dc data/stacks_exports/populations.snps.vcf.gz | awk '/MVZ216210/{found=1} END{exit !found}'; then
  fail "MVZ216210 remains in the final VCF"
fi
if grep -E -q '(^|[[:space:]])MVZ216210([[:space:];,()]|$)' \
  data/metadata/popmap_excluding_MVZ216210.tsv \
  data/unlinked_snps/unlinked.fam \
  data/phylogeny/trow_individuals_biallelic_iqtree.varsites.phy \
  data/phylogeny/trow_individuals_biallelic_iqtree.treefile \
  data/phylogeny/trow_individuals_biallelic_iqtree.contree; then
  fail "MVZ216210 remains in a final analysis input"
fi

largest=$(find . -path './.git' -prune -o -type f -exec stat -f '%z %N' {} \; | sort -nr | sed -n '1p')
largest_bytes=${largest%% *}
test "$largest_bytes" -lt 100000000 || fail "a file exceeds GitHub's 100 MB limit: $largest"

if test -f SHA256SUMS.txt; then
  shasum -a 256 -c SHA256SUMS.txt >/dev/null || fail "checksum validation failed"
fi

printf 'Repository validation passed.\n'
printf 'VCF: %s samples, %s variants\n' "$vcf_samples" "$vcf_variants"
printf 'Unlinked dataset: %s samples, %s SNPs\n' "$unlinked_samples" "$unlinked_snps"
printf 'Georeferenced samples: %s\n' "$coord_samples"
printf 'Largest file: %s\n' "$largest"
