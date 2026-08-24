#!/usr/bin/env bash
set -euo pipefail

ROOT="${DOWNSTREAM_ROOT:-/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23}"
BASE="${STACKS_BASE_DIR:-$(cd "$ROOT/.." && pwd)}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
TREEMIX_PLOTTING_FUNCS="${TREEMIX_PLOTTING_FUNCS:-$(cd "$BASE/.." && pwd)/.envs/treemix_x86/bin/plotting_funcs.R}"

fail() { echo "ERROR: $*" >&2; exit 2; }

[[ -f "$BASE/popmap.tsv" ]] || fail "Missing popmap.tsv at $BASE/popmap.tsv"
[[ -f "$BASE/populations.snps.vcf" ]] || fail "Missing populations.snps.vcf at $BASE/populations.snps.vcf"
[[ -f "$TREEMIX_PLOTTING_FUNCS" ]] || fail "Missing TreeMix plotting helper at $TREEMIX_PLOTTING_FUNCS"

command -v "$PYTHON_BIN" >/dev/null 2>&1 || fail "Python executable not found: $PYTHON_BIN"

"$PYTHON_BIN" - <<'PY'
mods = ["numpy", "pandas", "scipy", "matplotlib", "cartopy", "dadi"]
missing = []
for m in mods:
    try:
        __import__(m)
    except Exception:
        missing.append(m)
if missing:
    raise SystemExit("Missing Python modules: " + ", ".join(missing))
print("Python dependency check: OK")
PY

Rscript - <<'RS'
req <- c("ape")
miss <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Missing R packages: ", paste(miss, collapse = ", "))
cat("R dependency check: OK\n")
RS

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE_BIN="$(command -v iqtree2)"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE_BIN="$(command -v iqtree)"
else
  fail "Missing IQ-TREE executable (iqtree2 or iqtree)"
fi

echo "IQ-TREE executable: $IQTREE_BIN"
echo "Preflight runtime dependency check: OK"
