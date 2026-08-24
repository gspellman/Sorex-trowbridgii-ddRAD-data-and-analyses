#!/usr/bin/env bash
set -u

ROOT='/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean/downstream_analysis_2026-02-23'
BASE='/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean'
LOG_DIR="$ROOT/logs/rerun_excluding_MVZ216210_clean"
STATUS="$LOG_DIR/rerun_status.tsv"
RUNNER_LOG="$LOG_DIR/runner.log"
mkdir -p "$LOG_DIR"

export DOWNSTREAM_ROOT="$ROOT"
export STACKS_BASE_DIR="$BASE"
export TREEMIX_PLOTTING_FUNCS='/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/.envs/treemix_x86/bin/plotting_funcs.R'
if [[ -x /tmp/demomodel_env/bin/python ]]; then
  export PYTHON_BIN="${PYTHON_BIN:-/tmp/demomodel_env/bin/python}"
else
  export PYTHON_BIN="${PYTHON_BIN:-python3}"
fi
export MPLCONFIGDIR="$ROOT/.mplconfig"
export XDG_CACHE_HOME="$ROOT/.cache"
mkdir -p "$MPLCONFIGDIR"
mkdir -p "$XDG_CACHE_HOME"

log() {
  local ts
  ts="$(date '+%Y-%m-%d %H:%M:%S')"
  echo "[$ts] $*" | tee -a "$RUNNER_LOG"
}

check_reqs() {
  local req
  for req in "$@"; do
    if [[ ! -f "$req" ]]; then
      log "MISSING required file: $req"
      return 1
    fi
  done
  return 0
}

run_step() {
  local label="$1"; shift
  local cmd="$1"; shift
  local reqs=("$@")
  local ts_start ts_end rc logf

  ts_start="$(date '+%Y-%m-%d %H:%M:%S')"
  logf="$LOG_DIR/${label}.log"
  log "START $label"

  if [[ ${#reqs[@]} -gt 0 ]]; then
    check_reqs "${reqs[@]}" >"$logf" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then
      ts_end="$(date '+%Y-%m-%d %H:%M:%S')"
      echo -e "$label\t$ts_start\t$ts_end\t$rc\t$logf" >> "$STATUS"
      log "END $label (rc=$rc; missing prerequisite)"
      return $rc
    fi
  fi

  bash -lc "$cmd" >"$logf" 2>&1
  rc=$?
  ts_end="$(date '+%Y-%m-%d %H:%M:%S')"
  echo -e "$label\t$ts_start\t$ts_end\t$rc\t$logf" >> "$STATUS"
  log "END $label (rc=$rc)"
  return $rc
}

echo -e "step\tstart\tend\texit_code\tlog" > "$STATUS"
: > "$RUNNER_LOG"

run_or_stop() {
  local label="$1"; shift
  local cmd="$1"; shift
  run_step "$label" "$cmd" "$@" || {
    log "Pipeline halted at $label"
    exit 1
  }
}

# Preflight checks (sample exclusion + runtime dependencies)
run_or_stop 00_preflight_exclusion_check "'$ROOT/scripts/preflight_exclusion_check_2026-02-26.sh'"
run_or_stop 00b_preflight_runtime "'$ROOT/scripts/preflight_runtime_dependencies_2026-02-27.sh'"

# Core downstream tables/figures used by many later steps
run_or_stop 01_run_downstream "Rscript '$ROOT/scripts/run_downstream.R'" \
  "$BASE/populations.snps.vcf" \
  "$BASE/popmap.tsv"
run_or_stop 01b_run_treemix_models "Rscript '$ROOT/scripts/run_treemix_models_2026-02-27.R'" \
  "$BASE/populations.snps.vcf" \
  "$BASE/popmap.tsv"

# Core visual updates and summary stats
run_or_stop 02_redraw_pca "Rscript '$ROOT/scripts/redraw_pca_matrix.R'" \
  "$ROOT/tables/pca_variance_explained.tsv"
run_or_stop 03_run_fst_fis "Rscript '$ROOT/scripts/run_fst_fis_permutations.R'" \
  "$BASE/populations.snps.vcf" \
  "$BASE/popmap.tsv"
run_or_stop 04_redraw_fst "Rscript '$ROOT/scripts/redraw_fst_matrix.R'" \
  "$ROOT/tables/pairwise_fst_weir_cockerham.tsv"
run_or_stop 05_treemix_figs "Rscript '$ROOT/scripts/run_treemix_figs.R'" \
  "$ROOT/work/treemix/tm_m4.treeout.gz"
run_or_stop 06_splitstree_proxy "Rscript '$ROOT/scripts/run_splitstree_proxy.R'" \
  "$BASE/populations.snps.vcf"

# ML phylogeny prerequisites and plotting
run_or_stop 07_build_ml_alignment "${PYTHON_BIN} '$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/scripts/build_biallelic_alignment.py'" \
  "$BASE/populations.snps.vcf" \
  "$BASE/popmap.tsv"
run_or_stop 08_run_iqtree_ml "'$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/scripts/run_iqtree_ml_bootstrap.sh'" \
  "$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/work/biallelic_individuals_relaxed.phy"
run_or_stop 09_plot_ml_unrooted "Rscript '$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/scripts/plot_ml_tree_unrooted.R'" \
  "$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/work/ml_biallelic_individuals.contree"
run_or_stop 10_phylogeny_compoplot "Rscript '$ROOT/scripts/run_individual_phylogeny_compoplot.R'" \
  "$ROOT/tables/admixture_nmf_Qmatrix_bestK.tsv" \
  "$ROOT/additional_analyses_2026-02-24/phylogeny_ml_biallelic_2026-02-24/work/ml_biallelic_individuals.contree"

# Remaining phylogenetic/network analyses
run_or_stop 11_astral "Rscript '$ROOT/astral_analysis_2026-02-24/scripts/run_astral_from_vcf.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 12_unrooted_network "Rscript '$ROOT/all_samples_unrooted_network_2026-02-24/scripts/run_unrooted_sample_network.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 13_neighbornet "Rscript '$ROOT/all_samples_unrooted_network_2026-02-24/scripts/run_neighbornet_all_samples.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 14_prepare_ml_network_inputs "'$ROOT/ml_network_2026-02-24/scripts/prepare_ml_network_inputs.sh'" \
  "$ROOT/work/treemix/tm_m4.treeout.gz"
run_or_stop 15_ml_network "Rscript '$ROOT/ml_network_2026-02-24/scripts/plot_treemix_ml_network.R'" \
  "$ROOT/ml_network_2026-02-24/work/treemix_m4.treeout.gz"

# Demography
run_or_stop 16_demographic_summary "Rscript '$ROOT/scripts/run_demographic_history_summary.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 17_effective_size "Rscript '$ROOT/scripts/run_effective_size_skyline_proxy.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 18_beast_skyline "Rscript '$ROOT/scripts/run_beast_bayesian_skyline.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 19_psmc_proxy "Rscript '$ROOT/psmc_analysis_2026-02-24/scripts/run_psmc_rad_proxy.R'" \
  "$BASE/populations.snps.vcf"
run_or_stop 20_beast_species_tree "Rscript '$ROOT/beast_species_tree_2026-02-24/scripts/run_beast_species_tree.R'" \
  "$BASE/populations.snps.vcf"

# Additional analyses
run_or_stop 21_additional_1236 "${PYTHON_BIN} '$ROOT/additional_analyses_2026-02-24/scripts/run_additional_analyses_1_2_3_6.py'" \
  "$BASE/populations.snps.vcf"
run_or_stop 22_fst_manhattan_overlap "${PYTHON_BIN} '$ROOT/additional_analyses_2026-02-24/scripts/run_pairwise_fst_manhattan_and_overlap.py'" \
  "$BASE/populations.snps.vcf"
run_or_stop 23_locus_discordance "Rscript '$ROOT/additional_analyses_2026-02-24/scripts/run_locus_tree_discordance_7.R'" \
  "$ROOT/tables/individual_phylogeny_iqtree_midpoint.newick"
run_or_stop 24_maps_ibd "${PYTHON_BIN} '$ROOT/additional_analyses_2026-02-24/scripts/run_maps_and_ibd_from_coordinates.py'" \
  "$ROOT/tables/admixture_nmf_Qmatrix_bestK.tsv"

# Tables and missingness
run_or_stop 25_stacks_tables "Rscript '$ROOT/scripts/make_stacks_publication_tables.R'"
run_or_stop 26_pop_snp_table "Rscript '$ROOT/scripts/make_population_snp_diversity_table.R'"
run_or_stop 27_pop_snp_table_render "Rscript '$ROOT/scripts/render_population_snp_diversity_table_figure.R'" \
  "$ROOT/tables/Table3_population_SNP_diversity_publication.tsv"
run_or_stop 28_missingness_outliers "Rscript '$ROOT/scripts/run_sample_missingness_outliers_2026-02-26.R'" \
  "$BASE/populations.snps.vcf"

# Species delimitation figure chain
run_or_stop 29_speede_rf_delimitr "${PYTHON_BIN} '$ROOT/scripts/run_speede_rf_delimitr_unlinked_2026-02-26.py'" \
  "$BASE/populations.snps.vcf"
run_or_stop 30_species_structure_figure "${PYTHON_BIN} '$ROOT/scripts/make_species_structure_summary_figure_2026-02-26.py'" \
  "$ROOT/tables/admixture_nmf_Qmatrix_bestK.tsv"
run_or_stop 31_bpp_reconfigured "Rscript '$ROOT/scripts/run_bpp_reconfigured_2026-02-26.R'"
run_or_stop 32_bpp_fixedtree "Rscript '$ROOT/scripts/run_bpp_fixedtree_ess_2026-02-26.R'"
run_or_stop 33_bfdstar_bpp_unlinked "Rscript '$ROOT/scripts/run_bfdstar_bpp_unlinked_ess_2026-02-26.R'" \
  "$ROOT/species_delimitation_speede_rf_delimitr_2026-02-26/tables/unlinked_one_snp_per_locus_metadata.tsv"

log "Rerun complete: $(date '+%Y-%m-%d %H:%M:%S')"
