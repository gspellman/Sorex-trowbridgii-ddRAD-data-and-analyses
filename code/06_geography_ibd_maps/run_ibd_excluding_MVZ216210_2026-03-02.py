#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

BASE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
DOWN = os.path.join(BASE, "downstream_analysis_2026-02-23")
FIG = os.path.join(DOWN, "figures")
TAB = os.path.join(DOWN, "tables")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

COORD = "/Users/gspellman/Trowbridgii_analyses/Sample_geographic_coordinants.txt"
POPMAP = os.path.join(BASE, "popmap.tsv")
VCF = os.path.join(BASE, "populations.snps.vcf")

OUT_STATS = os.path.join(TAB, "IBD_statistics_summary_excl_MVZ216210.tsv")
OUT_PAIRS = os.path.join(TAB, "IBD_pairwise_distances_excl_MVZ216210.tsv")
OUT_BINS = os.path.join(TAB, "IBD_distance_bin_summary_excl_MVZ216210.tsv")
OUT_SUMMARY = os.path.join(TAB, "IBD_analysis_summary_excl_MVZ216210.txt")

OUT_SCATTER_PNG = os.path.join(FIG, "IBD_scatter_regression_excl_MVZ216210_publication.png")
OUT_SCATTER_PDF = os.path.join(FIG, "IBD_scatter_regression_excl_MVZ216210_publication.pdf")
OUT_BIN_PNG = os.path.join(FIG, "IBD_distance_bin_trend_excl_MVZ216210_publication.png")
OUT_BIN_PDF = os.path.join(FIG, "IBD_distance_bin_trend_excl_MVZ216210_publication.pdf")


def parse_gt(gt_field):
    gt = gt_field.split(':', 1)[0].replace('|', '/')
    if gt in ("./.", "."):
        return np.nan
    if gt == "0/0":
        return 0.0
    if gt in ("0/1", "1/0"):
        return 1.0
    if gt == "1/1":
        return 2.0
    return np.nan


def read_vcf_matrix(vcf_path, sample_set):
    chrom, pos, vid, rows = [], [], [], []
    samples = None
    with open(vcf_path, 'r') as fh:
        sample_idx = None
        for ln in fh:
            if ln.startswith('##'):
                continue
            if ln.startswith('#CHROM'):
                h = ln.rstrip('\n').split('\t')
                all_samples = h[9:]
                samples = [s for s in all_samples if s in sample_set]
                sample_idx = [all_samples.index(s) for s in samples]
                continue
            p = ln.rstrip('\n').split('\t')
            ref, alt = p[3], p[4]
            if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                continue
            g = p[9:]
            row = [parse_gt(g[j]) for j in sample_idx]
            chrom.append(p[0])
            pos.append(int(p[1]))
            vid.append(p[2] if p[2] != '.' else f"{p[0]}:{p[1]}")
            rows.append(row)
    return samples, np.array(chrom), np.array(pos), np.array(vid), np.array(rows, dtype=float)


def haversine_km(lat1, lon1, lat2, lon2):
    r = 6371.0088
    p1 = np.radians(lat1)
    p2 = np.radians(lat2)
    dphi = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.0) ** 2
    return 2.0 * r * np.arcsin(np.sqrt(a))


def mantel_test(d1, d2, perms=4999, seed=20260302):
    rng = np.random.default_rng(seed)
    n = d1.shape[0]
    triu = np.triu_indices(n, k=1)
    x = d1[triu]
    y = d2[triu]
    r_obs = np.corrcoef(x, y)[0, 1]
    c = 0
    for _ in range(perms):
        p = rng.permutation(n)
        yp = d2[p][:, p][triu]
        r = np.corrcoef(x, yp)[0, 1]
        if abs(r) >= abs(r_obs):
            c += 1
    pval = (c + 1) / (perms + 1)
    return r_obs, pval


def main():
    for req in (COORD, POPMAP, VCF):
        if not os.path.exists(req):
            raise FileNotFoundError(req)

    coords = pd.read_csv(COORD, sep='\t')
    coords = coords.dropna(subset=['Sample', 'Latitude', 'Longitude']).copy()
    coords['Sample'] = coords['Sample'].astype(str)

    pop = pd.read_csv(POPMAP, sep='\t', header=None, names=['Sample', 'Population'])
    pop['Sample'] = pop['Sample'].astype(str)

    if 'MVZ216210' in set(pop['Sample']):
        raise RuntimeError('MVZ216210 still present in popmap for excluded dataset.')

    meta = pop.merge(coords, on='Sample', how='inner').copy()
    meta = meta.sort_values(['Population', 'Sample']).reset_index(drop=True)

    sample_set = set(meta['Sample'])
    samples, _chrom, _pos, _vid, G = read_vcf_matrix(VCF, sample_set)
    if samples is None or len(samples) < 3:
        raise RuntimeError('Insufficient samples parsed from VCF for IBD analysis.')

    # align metadata to genotype sample order
    meta2 = meta.set_index('Sample').loc[samples].reset_index()

    # standard SNP filters for distance analysis
    site_missing = np.mean(np.isnan(G), axis=1)
    p = np.nanmean(G, axis=1) / 2.0
    maf = np.minimum(p, 1.0 - p)
    keep = (site_missing <= 0.20) & np.isfinite(maf) & (maf >= 0.05)
    G = G[keep]
    if G.shape[0] < 10:
        raise RuntimeError('Too few SNPs remain after filtering for IBD analysis.')

    n = len(samples)
    geo = np.zeros((n, n), dtype=float)
    gen = np.zeros((n, n), dtype=float)
    pair_rows = []

    for i in range(n - 1):
        gi = G[:, i]
        for j in range(i + 1, n):
            gj = G[:, j]
            ok = np.isfinite(gi) & np.isfinite(gj)
            if np.any(ok):
                dgen = float(np.mean(np.abs(gi[ok] - gj[ok]) / 2.0))
            else:
                dgen = np.nan

            dgeo = float(haversine_km(
                meta2.loc[i, 'Latitude'], meta2.loc[i, 'Longitude'],
                meta2.loc[j, 'Latitude'], meta2.loc[j, 'Longitude']
            ))

            gen[i, j] = gen[j, i] = dgen
            geo[i, j] = geo[j, i] = dgeo

            pair_rows.append({
                'Sample1': samples[i],
                'Sample2': samples[j],
                'Population1': meta2.loc[i, 'Population'],
                'Population2': meta2.loc[j, 'Population'],
                'Geographic_distance_km': dgeo,
                'Genetic_distance': dgen,
            })

    mgen = np.nanmedian(gen[np.isfinite(gen)])
    gen[~np.isfinite(gen)] = mgen

    triu = np.triu_indices(n, k=1)
    x = geo[triu]
    y = gen[triu]

    pear_r, pear_p = stats.pearsonr(x, y)
    spear_r, spear_p = stats.spearmanr(x, y)
    mantel_r, mantel_p = mantel_test(gen, geo, perms=4999)

    lm = stats.linregress(x, y)

    stats_df = pd.DataFrame([{
        'Dataset': 'stacks_refmap_sorex_excl_MVZ216210',
        'Excluded_sample': 'MVZ216210',
        'N_samples': n,
        'N_pairwise_comparisons': len(x),
        'N_SNPs_after_filtering': int(G.shape[0]),
        'Pearson_r': float(pear_r),
        'Pearson_p': float(pear_p),
        'Spearman_rho': float(spear_r),
        'Spearman_p': float(spear_p),
        'Mantel_r': float(mantel_r),
        'Mantel_p_permutation': float(mantel_p),
        'Mantel_permutations': 4999,
        'Linear_slope': float(lm.slope),
        'Linear_intercept': float(lm.intercept),
        'Linear_r2': float(lm.rvalue ** 2),
        'Linear_p': float(lm.pvalue),
        'Linear_slope_stderr': float(lm.stderr),
    }])
    stats_df.to_csv(OUT_STATS, sep='\t', index=False)

    pairs_df = pd.DataFrame(pair_rows)
    pairs_df.to_csv(OUT_PAIRS, sep='\t', index=False)

    bins = np.quantile(x, np.linspace(0, 1, 11))
    bins = np.unique(bins)
    idx = np.digitize(x, bins[1:-1], right=True)
    b_rows = []
    for b in range(len(bins) - 1):
        m = idx == b
        if np.sum(m) == 0:
            continue
        b_rows.append({
            'Bin': b + 1,
            'Geo_km_min': float(bins[b]),
            'Geo_km_max': float(bins[b + 1]),
            'N_pairs': int(np.sum(m)),
            'Mean_genetic_distance': float(np.mean(y[m])),
            'SE_genetic_distance': float(np.std(y[m], ddof=1) / np.sqrt(np.sum(m))) if np.sum(m) > 1 else np.nan,
        })
    bin_df = pd.DataFrame(b_rows)
    bin_df.to_csv(OUT_BINS, sep='\t', index=False)

    for out in (OUT_SCATTER_PNG, OUT_SCATTER_PDF):
        plt.figure(figsize=(8.2, 6.4), dpi=400 if out.endswith('.png') else None)
        ax = plt.gca()
        ax.scatter(x, y, s=16, c='#4C78A8', alpha=0.40, edgecolors='none')
        xs = np.linspace(np.nanmin(x), np.nanmax(x), 300)
        ax.plot(xs, lm.intercept + lm.slope * xs, color='#D62728', linewidth=2.1)
        txt = (
            f"Pearson r = {pear_r:.3f} (p = {pear_p:.2e})\\n"
            f"Spearman rho = {spear_r:.3f} (p = {spear_p:.2e})\\n"
            f"Mantel r = {mantel_r:.3f} (p = {mantel_p:.4f}, 4,999 perms)\\n"
            f"Linear slope = {lm.slope:.3e}"
        )
        ax.text(
            0.02, 0.98, txt, transform=ax.transAxes, va='top', ha='left', fontsize=10,
            bbox=dict(facecolor='white', alpha=0.88, edgecolor='#BFBFBF')
        )
        ax.set_xlabel('Geographic distance (km)', fontsize=11)
        ax.set_ylabel('Genetic distance (mean allele-dosage difference / 2)', fontsize=11)
        ax.set_title('Isolation by distance (MVZ216210 excluded)', fontsize=12)
        ax.grid(color='#E0E0E0', linestyle=':', linewidth=0.9)
        plt.tight_layout()
        plt.savefig(out)
        plt.close()

    for out in (OUT_BIN_PNG, OUT_BIN_PDF):
        plt.figure(figsize=(7.8, 5.6), dpi=400 if out.endswith('.png') else None)
        xc = (bin_df['Geo_km_min'] + bin_df['Geo_km_max']) / 2.0
        plt.errorbar(
            xc, bin_df['Mean_genetic_distance'], yerr=bin_df['SE_genetic_distance'],
            fmt='o-', color='#1F77B4', ecolor='#1F77B4', elinewidth=1.4, capsize=3.0, linewidth=1.6
        )
        plt.xlabel('Geographic distance bin midpoint (km)', fontsize=11)
        plt.ylabel('Mean genetic distance', fontsize=11)
        plt.title('IBD distance-class trend (MVZ216210 excluded)', fontsize=12)
        plt.grid(color='#E0E0E0', linestyle=':', linewidth=0.9)
        plt.tight_layout()
        plt.savefig(out)
        plt.close()

    with open(OUT_SUMMARY, 'w') as fh:
        fh.write('Isolation-by-distance analysis summary (MVZ216210 excluded)\n')
        fh.write(f'Samples analyzed: {n}\n')
        fh.write(f'SNPs after filters: {G.shape[0]}\n')
        fh.write(f'Pairwise comparisons: {len(x)}\n')
        fh.write(f'Pearson r={pear_r:.6f}, p={pear_p:.6e}\n')
        fh.write(f'Spearman rho={spear_r:.6f}, p={spear_p:.6e}\n')
        fh.write(f'Mantel r={mantel_r:.6f}, p={mantel_p:.6f}, permutations=4999\n')
        fh.write(f'Scatter figure: {OUT_SCATTER_PNG}\n')
        fh.write(f'Binned figure: {OUT_BIN_PNG}\n')
        fh.write(f'Stats table: {OUT_STATS}\n')
        fh.write(f'Pairwise table: {OUT_PAIRS}\n')
        fh.write(f'Binned table: {OUT_BINS}\n')

    print('IBD analysis complete', flush=True)


if __name__ == '__main__':
    main()
