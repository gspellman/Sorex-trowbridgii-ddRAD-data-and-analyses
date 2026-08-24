#!/usr/bin/env bash
set -euo pipefail

ROOT="${DOWNSTREAM_ROOT:-/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23}"
WORK="$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/work"
ALIGN="$WORK/biallelic_individuals_relaxed.phy"
PREFIX="$WORK/ml_biallelic_individuals"

if [[ ! -f "$ALIGN" ]]; then
  echo "ERROR: missing alignment: $ALIGN" >&2
  exit 2
fi

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE_BIN="$(command -v iqtree2)"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE_BIN="$(command -v iqtree)"
else
  echo "ERROR: iqtree2/iqtree not found" >&2
  exit 2
fi

"$IQTREE_BIN" \
  -s "$ALIGN" \
  -m GTR+F+R4 \
  -B 1000 \
  --bnni \
  -T AUTO \
  --prefix "$PREFIX"

for ext in treefile contree iqtree log mldist; do
  [[ -f "$PREFIX.$ext" ]] || { echo "ERROR: missing IQ-TREE output $PREFIX.$ext" >&2; exit 2; }
done

echo "IQ-TREE run complete: $PREFIX"
