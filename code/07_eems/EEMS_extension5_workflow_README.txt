EEMS extension-5 workflow (dataset excluding MVZ216210)

Final execution order:
1. run_eems_ess200_primary_2026-07-29.R
2. run_eems_ess200_extension1_2026-07-29.R
3. run_eems_ess200_extension2_2026-07-30.R
4. run_eems_ess200_extension3_2026-07-30.R
5. run_eems_ess200_extension4_2026-08-04.R
6. run_eems_ess200_extension5_2026-08-05.R
7. diagnose_eems_ess200_extension5_2026-08-05.R
8. plot_eems_extension5_publication_2026-08-10.R

The extension-5 runner performs a 2-million-iteration tuning phase followed by
a 20-million-iteration production phase for each of four chains. Production
uses qVoronoiPr=0.30, mEffctProposalS2=0.15, 2 million burn-in iterations,
thin=999, and approximately 18,000 retained states per chain.

Final diagnostics: minimum ESS=318.1076; maximum rank-normalized split
R-hat=1.047568. 
