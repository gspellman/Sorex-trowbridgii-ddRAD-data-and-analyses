#!/usr/bin/env bash
set -euo pipefail
ROOT='/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean'
S='MVZ216210'
if grep -q "^${S}[[:space:]]" "$ROOT/popmap.tsv"; then
  echo "ERROR: $S present in popmap.tsv"; exit 2
fi
python3 - <<'PY'
vcf='/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/populations.snps.vcf'
s='MVZ216210'
with open(vcf) as f:
    for ln in f:
        if ln.startswith('#CHROM'):
            h=ln.rstrip('\n').split('\t')
            if s in h: raise SystemExit('ERROR: sample present in VCF header')
            print('Preflight OK: sample excluded from VCF and popmap.')
            break
PY
