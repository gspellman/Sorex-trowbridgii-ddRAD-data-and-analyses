# EEMS Methods and Results Summary

## Methods

Spatial variation in effective migration and effective diversity was estimated with the SNP implementation of estimated effective migration surfaces (`runeems_snps`; Petkova et al., 2016). The analysis used 98 georeferenced individuals and 14,578 filtered SNPs from the final reference-mapped ddRAD dataset. Sample MVZ216210 was excluded before constructing the genetic dissimilarity matrix because it was identified as an extreme missing-data outlier. Input validation confirmed that the 98 x 98 dissimilarity matrix was symmetric, had a zero diagonal, and matched the sample and coordinate order.

The prediction domain was restricted to terrestrial habitat. Sample coordinates were transformed to NAD83/Conus Albers (EPSG:5070), and a convex hull around sampled localities was buffered by 100 km. This polygon was intersected with 1:10-million Natural Earth land polygons for the United States and Canada, buffered by 1 km to retain coastal samples, and simplified at 1-km tolerance while preserving topology. The resulting polygon included every sample but excluded marine and major freshwater areas from prediction.

Three independent diploid EEMS chains were run with random seeds 20260729, 20260730, and 20260731. To assess sensitivity to the population grid, chains used 200, 250, and 300 demes, respectively. Each chain comprised 1,000,000 MCMC iterations, of which the first 500,000 were discarded as burn-in; every 499th post-burn-in state was retained, yielding 1,000 posterior states per chain and 3,000 states overall. Posterior migration and diversity surfaces were averaged across chains and deme grids with `rEEMSplots`. Trace effective sample sizes (ESS) for the log prior and log likelihood were estimated from paired autocorrelations using an initial-positive-sequence truncation rule.

Publication maps were generated from the posterior-averaged raster grids. Surfaces were explicitly masked to the terrestrial prediction polygon, plotted at the correct geographic aspect, and overlaid with coastlines, state or provincial boundaries, and population-coloured sample locations. Migration and diversity are shown as centred log10 relative rates, such that blue indicates values above the spatial mean, orange indicates values below the spatial mean, and white indicates values near the mean. Vector PDF and 600-dpi PNG versions were produced.

## Results

All three EEMS chains completed successfully and each produced 1,000 retained posterior states. Runtime was 26.5 min for the 200-deme chain, 35.3 min for the 250-deme chain, and 39.4 min for the 300-deme chain. Proposal acceptance rates were generally moderate for migration-rate and tile-movement updates, whereas migration birth-death proposals had lower acceptance (8.8-9.5%).

The posterior mean migration surface indicated pronounced spatial heterogeneity. Higher-than-average effective migration was inferred across the northern portion of the sampled range and in several localized regions centred on sampled population clusters. Broad regions of lower-than-average migration separated portions of the northern, central, and southern sampling distribution, consistent with spatially variable connectivity rather than a single homogeneous isolation-by-distance surface. The effective-diversity surface was comparatively smooth, with higher relative values across much of central California and lower relative values toward parts of the northern and southern range margins.

Trace autocorrelation remained appreciable. Per-chain ESS values were 41.4-47.3 for the log prior and 19.7-56.3 for the log likelihood. Consequently, broad spatial patterns were more repeatable and interpretable than fine-scale local features, and individual small migration peaks or troughs should be treated as provisional unless confirmed by longer chains or independent data.

## Reference

Petkova, D., Novembre, J. & Stephens, M. 2016. Visualizing spatial population structure with estimated effective migration surfaces. *Nature Genetics*, 48, 94-100. doi:10.1038/ng.3464
