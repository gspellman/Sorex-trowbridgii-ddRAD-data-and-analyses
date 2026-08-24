# Trowbridgii ddRAD analysis code

Publication-release code archive for the reference-mapped ddRAD analysis of
the dataset excluding the missingness outlier `MVZ216210`.

## Scope

This archive contains the custom analysis and figure-production code used for
the final reported analyses. Raw reads, reference-genome files, intermediate
data, final result tables, and third-party software source trees are not
duplicated here. The scripts expect those files to be supplied separately.

The scripts preserve the paths and parameters used for the original analysis.
Before running them on another computer, replace the project-root paths near
the start of each script or expose the same directory structure.

## Folder organization

| Folder | Contents |
|---|---|
| `01_stacks_qc_snp_calling` | STACKS reporting, dependency checks, missingness screening, and verification of the `MVZ216210` exclusion |
| `02_population_structure` | PCA, DAPC, fastStructure, PopCluster, and core SNP clustering analyses |
| `03_phylogenetics` | Biallelic alignment preparation, IQ-TREE maximum-likelihood inference, bootstrap analysis, and tree plotting |
| `04_species_delimitation` | SPEEDEMON, delimitR, and BPP workflows using unlinked SNPs and the IQ-TREE guide topology |
| `05_population_statistics_sfs_fst_fis` | Population diversity, site-frequency spectra, FST/FIS permutations, genomic scans, and divergence-overlap tests |
| `06_geography_ibd_maps` | Geographic ancestry maps and corrected isolation-by-distance analyses |
| `07_eems` | Complete primary-to-Extension-5 EEMS restart lineage, diagnostics, and final publication plotting |
| `08_tables_figures` | Final publication table and composite-figure rendering |
| `09_configuration_files` | BPP control files, BEAST XML files, and final EEMS parameter files |
| `99_workflow` | Top-level downstream and clean-exclusion workflow entry points |

Each analysis folder has an `INDEX.txt`. `MANIFEST.tsv` records every released
file, and `SHA256SUMS.txt` permits integrity checks after transfer.

## Final dataset

- Reference genome: *Sorex ornatus*, NCBI assembly GCA_041430635.1.
- A priori populations: North, North_Coast, Sierra_1, Sierra_2, Sierra_3,
  and South_Coast.
- `MVZ216210` was removed as the extreme missingness outlier before final
  downstream analyses.
- The core final downstream dataset contains 99 individuals; analyses needing
  complete geographic coordinates contain 98 individuals.

## Important interpretation notes

- BFD* is not included because its runs did not meet the required ESS target
  and it was removed from the final reported results.
- NeighborNet is not included because it was removed from the final Methods.
- Standalone Random Forest species-delimitation results are not included in
  final summaries or figures. `run_speedemon_delimitr_unlinked.py` retains the
  original random-forest helper because the delimitR-style model-choice step
  uses random-forest classification internally.
- EEMS primary and Extensions 1-4 scripts are retained because they document
  the restart lineage leading to the accepted Extension 5 estimate. Superseded
  MCMC output directories are not included.
- Third-party programs must be obtained from their original distributions and
  cited according to the manuscript Methods.

## Suggested execution order

1. Run the scripts in `01_stacks_qc_snp_calling` and verify that `MVZ216210` is
   absent from the final popmap and SNP exports.
2. Run `99_workflow/run_downstream.R` for core population-genomic outputs.
3. Run analysis-specific scripts in folders `02` through `06`.
4. For EEMS, follow `07_eems/EEMS_extension5_workflow_README.txt` in numerical
   restart order, ending with the Extension 5 diagnostic and plotting scripts.
5. Generate final tables and composite figures with folder `08`.

## Reproducibility

Random seeds and analysis settings are retained in the scripts and control
files. Software versions should be recorded from the publication environment
when the archive is deposited. Codex was used to help correct code and ensure
that analyses ran smoothly; all analytical decisions, input validation, and
interpretation remain the responsibility of the authors.

