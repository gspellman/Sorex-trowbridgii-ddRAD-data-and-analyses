#!/usr/bin/env python3
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path('/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean')
OUT = BASE / 'downstream_analysis_2026-02-23'
TAB = OUT / 'tables'

VCF = BASE / 'populations.snps.vcf'
POPMAP = BASE / 'popmap.tsv'
TABLE2 = TAB / 'Table2_population_SNP_variation_summary.tsv'
DEMOG = TAB / 'demography_population_summary.tsv'
FSTMAT = TAB / 'pairwise_fst_weir_cockerham.tsv'
FIS = TAB / 'population_fis_permutation.tsv'

out_tsv = TAB / 'Table3_population_SNP_diversity_publication.tsv'
out_md = TAB / 'Table3_population_SNP_diversity_publication.md'
out_dxy = TAB / 'pairwise_dxy_matrix.tsv'

pop_order = ['North', 'North_Coast', 'Sierra_1', 'Sierra_2', 'Sierra_3', 'South_Coast']


def parse_gt(gt):
    g = gt.split(':', 1)[0].replace('|', '/')
    if g in ('./.', '.', '.|.'):
        return np.nan
    if g == '0/0':
        return 0.0
    if g in ('0/1', '1/0'):
        return 1.0
    if g == '1/1':
        return 2.0
    return np.nan


# metadata
popmap = pd.read_csv(POPMAP, sep='\t', header=None, names=['Sample', 'Population'])
popmap = popmap[popmap['Population'].isin(pop_order)].copy()
pop_lookup = dict(zip(popmap['Sample'], popmap['Population']))

# read VCF as dosage matrix
samples = None
rows = []
with open(VCF, 'r') as f:
    for ln in f:
        if ln.startswith('##'):
            continue
        if ln.startswith('#CHROM'):
            parts = ln.rstrip('\n').split('\t')
            all_s = parts[9:]
            samples = [s for s in all_s if s in pop_lookup]
            keep_idx = [all_s.index(s) for s in samples]
            continue
        p = ln.rstrip('\n').split('\t')
        ref, alt = p[3], p[4]
        if len(ref) != 1 or len(alt) != 1 or ',' in alt:
            continue
        g = p[9:]
        rows.append([parse_gt(g[j]) for j in keep_idx])

G = np.array(rows, dtype=float)

# global analysis filters used across downstream analyses
site_missing = np.mean(np.isnan(G), axis=1)
maf = np.nanmean(G, axis=1) / 2.0
maf = np.minimum(maf, 1 - maf)
mask = (site_missing <= 0.20) & np.isfinite(maf) & (maf >= 0.05)
G = G[mask]

# per-population sample indices
pidx = {p: np.array([i for i, s in enumerate(samples) if pop_lookup[s] == p], dtype=int) for p in pop_order}

# pairwise Dxy matrix
n_pop = len(pop_order)
dxy_mat = pd.DataFrame(np.nan, index=pop_order, columns=pop_order)

for i, p1 in enumerate(pop_order):
    for j, p2 in enumerate(pop_order):
        if i == j:
            dxy_mat.loc[p1, p2] = np.nan
            continue
        if j < i:
            dxy_mat.loc[p1, p2] = dxy_mat.loc[p2, p1]
            continue

        g1 = G[:, pidx[p1]]
        g2 = G[:, pidx[p2]]

        p1_alt = np.nanmean(g1, axis=1) / 2.0
        p2_alt = np.nanmean(g2, axis=1) / 2.0

        valid = np.isfinite(p1_alt) & np.isfinite(p2_alt)
        if np.sum(valid) == 0:
            dxy = np.nan
        else:
            dxy_site = p1_alt[valid] * (1.0 - p2_alt[valid]) + p2_alt[valid] * (1.0 - p1_alt[valid])
            dxy = float(np.nanmean(dxy_site))

        dxy_mat.loc[p1, p2] = dxy
        dxy_mat.loc[p2, p1] = dxy

# mean pairwise Dxy per population
mean_dxy = dxy_mat.mean(axis=1, skipna=True).rename('Dxy_mean_to_other_populations')

# mean pairwise Fst per population from existing Weir-Cockerham matrix
fst = pd.read_csv(FSTMAT, sep='\t', index_col=0)
fst = fst.reindex(index=pop_order, columns=pop_order)
for p in pop_order:
    fst.loc[p, p] = np.nan
mean_fst = fst.mean(axis=1, skipna=True).rename('Fst_mean_to_other_populations')

# base per-pop summary data
t2 = pd.read_csv(TABLE2, sep='\t')
dem = pd.read_csv(DEMOG, sep='\t')
fis = pd.read_csv(FIS, sep='\t')

m = (
    pd.DataFrame({'Population': pop_order})
    .merge(t2[['Population', 'N_samples', 'Total_biallelic_SNP_sites', 'Polymorphic_SNP_sites', 'Private', 'Mean_MAF_polymorphic_sites']], on='Population', how='left')
    .merge(dem[['Population', 'Segregating_sites', 'Mean_MAF', 'TajimasD_ascertained', 'Pi_per_site']], on='Population', how='left')
    .merge(fis[['Population', 'Fis']], on='Population', how='left')
)

# derive requested fields
m['SNPs_total'] = m['Segregating_sites']
m['SNPs_biallelic'] = m['Polymorphic_SNP_sites']

# allelic diversity represented as expected heterozygosity (gene diversity): 2p(1-p)
m['Allelic_diversity_total'] = 2.0 * m['Mean_MAF'] * (1.0 - m['Mean_MAF'])
m['Allelic_diversity_biallelic'] = 2.0 * m['Mean_MAF_polymorphic_sites'] * (1.0 - m['Mean_MAF_polymorphic_sites'])

m = m.set_index('Population')
m['Dxy_mean_to_other_populations'] = mean_dxy
m['Fst_mean_to_other_populations'] = mean_fst
m = m.reset_index()

final = m[[
    'Population',
    'N_samples',
    'SNPs_total',
    'SNPs_biallelic',
    'Allelic_diversity_total',
    'Allelic_diversity_biallelic',
    'Private',
    'TajimasD_ascertained',
    'Dxy_mean_to_other_populations',
    'Fst_mean_to_other_populations',
    'Fis',
    'Pi_per_site'
]].copy()

final = final.rename(columns={
    'N_samples': 'Sample_size',
    'Private': 'Private_alleles',
    'TajimasD_ascertained': "Tajimas_D",
    'Fis': 'Fis',
    'Pi_per_site': 'Pi_nucleotide_diversity'
})

# rounding for publication readability
for c in ['Allelic_diversity_total', 'Allelic_diversity_biallelic', 'Tajimas_D', 'Dxy_mean_to_other_populations', 'Fst_mean_to_other_populations', 'Fis', 'Pi_nucleotide_diversity']:
    final[c] = final[c].astype(float).round(6)

final.to_csv(out_tsv, sep='\t', index=False)
dxy_mat.to_csv(out_dxy, sep='\t', index=True)

# markdown table
md = final.copy()
md.columns = [
    'Population', 'Sample size', 'SNPs (total)', 'SNPs (biallelic)',
    'Allelic diversity (total)', 'Allelic diversity (biallelic)',
    'Private alleles', "Tajima's D", 'Dxy (mean to other pops)',
    'Fst (mean to other pops)', 'Fis', 'Pi (nucleotide diversity)'
]

with open(out_md, 'w') as f:
    f.write('# Table 3. SNP diversity summary by population\n\n')
    f.write(md.to_markdown(index=False))
    f.write('\n\n')
    f.write('Notes:\n')
    f.write('- `SNPs (total)` are segregating sites from the demographic summary table.\n')
    f.write('- `SNPs (biallelic)` are polymorphic biallelic SNP sites in the filtered panel.\n')
    f.write('- `Allelic diversity` is expected heterozygosity (2p(1-p)) based on the corresponding mean MAF.\n')
    f.write('- `Dxy` and `Fst` values are means of pairwise values between each focal population and the other five populations.\n')

print(f'Wrote: {out_tsv}')
print(f'Wrote: {out_md}')
print(f'Wrote: {out_dxy}')
