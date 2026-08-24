#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.path import Path
from scipy import stats
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

BASE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
OUT = os.path.join(BASE, "downstream_analysis_2026-02-23", "additional_analyses_2026-02-24")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
for d in (FIG, TAB):
    os.makedirs(d, exist_ok=True)

COORD = "/Users/gspellman/Trowbridgii_analyses/Sample_geographic_coordinants.txt"
POPMAP = os.path.join(BASE, "popmap.tsv")
QMAT = os.path.join(BASE, "downstream_analysis_2026-02-23", "tables", "admixture_nmf_Qmatrix_bestK.tsv")
CLMAP = os.path.join(BASE, "downstream_analysis_2026-02-23", "tables", "admixture_cluster_to_population_color_map.tsv")
VCF = os.path.join(BASE, "populations.snps.vcf")
TOPO_TIF = os.path.join(OUT, "work", "topography", "HYP_50M_SR_W.tif")

pop_cols = {
    "North": "#4169E1",
    "North_Coast": "#4CBB17",
    "Sierra_1": "#800080",
    "Sierra_2": "#00FFFF",
    "Sierra_3": "#FF0000",
    "South_Coast": "#000000",
}
pop_order = ["North", "North_Coast", "Sierra_1", "Sierra_2", "Sierra_3", "South_Coast"]


def parse_gt(gt):
    g = gt.split(':', 1)[0].replace('|', '/')
    if g in ("./.", "."):
        return np.nan
    if g == "0/0":
        return 0.0
    if g in ("0/1", "1/0"):
        return 1.0
    if g == "1/1":
        return 2.0
    return np.nan


def read_vcf_matrix(vcf_path, sample_set):
    chrom, pos, vid, rows = [], [], [], []
    with open(vcf_path, 'r') as f:
        sample_idx = None
        samples = None
        for ln in f:
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
    a = np.sin(dphi / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(a))


def _wedge_marker(theta1, theta2, n=36):
    th = np.linspace(np.deg2rad(theta1), np.deg2rad(theta2), n)
    verts = np.vstack((
        [0.0, 0.0],
        np.column_stack((np.cos(th), np.sin(th))),
        [0.0, 0.0],
    ))
    codes = [Path.MOVETO] + [Path.LINETO] * len(th) + [Path.CLOSEPOLY]
    return Path(verts, codes)


def draw_pie(ax, x, y, fracs, colors, size_pts2):
    start = 0.0
    for frac, col in zip(fracs, colors):
        if frac <= 0:
            continue
        theta1 = 360.0 * start
        theta2 = 360.0 * (start + frac)
        ax.scatter(
            [x], [y], marker=_wedge_marker(theta1, theta2), s=size_pts2, c=[col],
            edgecolors='white', linewidths=0.25, transform=ccrs.PlateCarree(), zorder=5
        )
        start += frac
    ax.scatter(
        [x], [y], marker='o', s=size_pts2, facecolors='none', edgecolors='white',
        linewidths=0.25, transform=ccrs.PlateCarree(), zorder=6
    )


def map_extent(df, pad_x=0.8, pad_y=0.6):
    xmin = float(df['Longitude'].min() - pad_x)
    xmax = float(df['Longitude'].max() + pad_x)
    ymin = float(df['Latitude'].min() - pad_y)
    ymax = float(df['Latitude'].max() + pad_y)
    return [xmin, xmax, ymin, ymax]


def add_geo_background(ax):
    if os.path.exists(TOPO_TIF):
        topo = mpimg.imread(TOPO_TIF)
        ax.imshow(
            topo, extent=[-180, 180, -90, 90], origin='upper',
            transform=ccrs.PlateCarree(), interpolation='bilinear', zorder=0, alpha=0.95
        )
    else:
        ax.stock_img()
    ax.add_feature(cfeature.LAND, facecolor=(0.95, 0.95, 0.92, 0.18), zorder=1)
    ax.add_feature(cfeature.OCEAN, facecolor=(0.80, 0.88, 0.95, 0.20), zorder=1)
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=0.8, edgecolor='#303030', zorder=2)
    ax.add_feature(cfeature.STATES.with_scale('10m'), linewidth=0.55, edgecolor='#4f4f4f', zorder=2)


def mantel_test(d1, d2, perms=1999, seed=20260224):
    rng = np.random.default_rng(seed)
    n = d1.shape[0]
    triu = np.triu_indices(n, k=1)
    x = d1[triu]
    y = d2[triu]
    r_obs = np.corrcoef(x, y)[0, 1]
    more = 0
    for _ in range(perms):
        p = rng.permutation(n)
        yp = d2[p][:, p][triu]
        r = np.corrcoef(x, yp)[0, 1]
        if abs(r) >= abs(r_obs):
            more += 1
    pval = (more + 1) / (perms + 1)
    return r_obs, pval


# ---------- Load and merge metadata ----------
coords = pd.read_csv(COORD, sep='\t')
coords = coords.dropna(subset=['Sample', 'Latitude', 'Longitude']).copy()
coords['Sample'] = coords['Sample'].astype(str)

pop = pd.read_csv(POPMAP, sep='\t', header=None, names=['Sample', 'Population'])
meta = pop.merge(coords, on='Sample', how='inner')
meta = meta[meta['Population'].isin(pop_order)].copy()
meta['Population'] = pd.Categorical(meta['Population'], categories=pop_order, ordered=True)
meta = meta.sort_values(['Population', 'Sample'])

meta.to_csv(os.path.join(TAB, 'sample_coordinates_with_population.tsv'), sep='\t', index=False)

# ---------- Map 1: sample locations ----------
extent = map_extent(meta)
proj = ccrs.LambertConformal(central_longitude=-121.0, central_latitude=41.5, standard_parallels=(33, 45))
for ext in ('png', 'pdf'):
    fig = plt.figure(figsize=(8.6, 7.2), dpi=320 if ext == 'png' else None)
    ax = plt.axes(projection=proj)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_aspect('equal', adjustable='box')
    add_geo_background(ax)
    for p in pop_order:
        d = meta[meta['Population'] == p]
        if len(d) == 0:
            continue
        ax.scatter(
            d['Longitude'], d['Latitude'], s=54, c=pop_cols[p], edgecolors='white', linewidths=0.45,
            label=p, alpha=0.96, transform=ccrs.PlateCarree(), zorder=4
        )

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Sample locations by population')
    gl = ax.gridlines(draw_labels=True, linewidth=0.45, color='#7a7a7a', alpha=0.25, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.x_inline = False
    gl.y_inline = False
    gl.rotate_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    gl.xformatter = LongitudeFormatter(number_format='.1f')
    gl.yformatter = LatitudeFormatter(number_format='.1f')
    ax.legend(loc='upper left', frameon=True, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f'Map_sample_locations_by_population_publication.{ext}'))
    plt.close()

# ---------- Map 2: admixture pies ----------
q = pd.read_csv(QMAT, sep='\t')
clmap = pd.read_csv(CLMAP, sep='\t')
cluster_cols = [c for c in q.columns if c.startswith('Cluster')]
cluster_color = dict(zip(clmap['Cluster'], clmap['Color']))
cluster_colors = [cluster_color.get(c, '#888888') for c in cluster_cols]

mq = meta.merge(q[['Sample'] + cluster_cols], on='Sample', how='left')
for c in cluster_cols:
    mq[c] = mq[c].fillna(0.0)

# radius based on map span
# adaptive pie size from nearest-neighbor spacing to reduce overlap
xy = mq[['Longitude', 'Latitude']].to_numpy(dtype=float)
if len(xy) > 1:
    lat0 = np.deg2rad(np.nanmean(xy[:, 1]))
    dlon = (xy[:, 0][:, None] - xy[:, 0][None, :]) * np.cos(lat0)
    dlat = (xy[:, 1][:, None] - xy[:, 1][None, :])
    dd = np.sqrt(dlon ** 2 + dlat ** 2)
    dd[dd == 0] = np.nan
    nn = np.nanmin(dd, axis=1)
    nn_med = float(np.nanmedian(nn))
else:
    nn_med = 0.15
pie_size = max(46.0, min(120.0, 1500.0 * nn_med))

for ext in ('png', 'pdf'):
    fig = plt.figure(figsize=(9.0, 7.4), dpi=320 if ext == 'png' else None)
    ax = plt.axes(projection=proj)
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.set_aspect('equal', adjustable='box')
    add_geo_background(ax)

    # light context points by population color
    for p in pop_order:
        d = mq[mq['Population'] == p]
        if len(d):
            ax.scatter(d['Longitude'], d['Latitude'], s=9, c=pop_cols[p], alpha=0.18, transform=ccrs.PlateCarree(), zorder=3)

    for _, r in mq.iterrows():
        fr = np.array([r[c] for c in cluster_cols], dtype=float)
        s = fr.sum()
        if s <= 0:
            continue
        fr = fr / s
        draw_pie(ax, r['Longitude'], r['Latitude'], fr, cluster_colors, pie_size)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title('Geographic distribution of admixture proportions (pie charts)')
    gl = ax.gridlines(draw_labels=True, linewidth=0.45, color='#7a7a7a', alpha=0.25, linestyle=':')
    gl.top_labels = False
    gl.right_labels = False
    gl.bottom_labels = True
    gl.left_labels = True
    gl.x_inline = False
    gl.y_inline = False
    gl.rotate_labels = False
    gl.xlabel_style = {'size': 9}
    gl.ylabel_style = {'size': 9}
    gl.xformatter = LongitudeFormatter(number_format='.1f')
    gl.yformatter = LatitudeFormatter(number_format='.1f')

    # cluster legend
    cl_handles = [plt.Line2D([0], [0], marker='s', linestyle='none', markerfacecolor=cluster_colors[i], markeredgecolor='white', markersize=7, label=cluster_cols[i])
                  for i in range(len(cluster_cols))]
    ax.legend(handles=cl_handles, title='Ancestry clusters', loc='lower left', fontsize=7, title_fontsize=8, frameon=True)

    fig.tight_layout()
    fig.savefig(os.path.join(FIG, f'Map_admixture_pies_by_sample_publication.{ext}'))
    plt.close()

# ---------- IBD analysis ----------
sample_set = set(meta['Sample'])
samples, chrom, pos, vid, G = read_vcf_matrix(VCF, sample_set)

# align meta to genotype order
meta2 = meta.set_index('Sample').loc[samples].reset_index()

# SNP filters for distance matrix
site_missing = np.mean(np.isnan(G), axis=1)
maf = np.minimum(np.nanmean(G, axis=1) / 2.0, 1 - np.nanmean(G, axis=1) / 2.0)
keep = (site_missing <= 0.20) & np.isfinite(maf) & (maf >= 0.05)
G = G[keep]

n = len(samples)
geo = np.zeros((n, n), dtype=float)
gen = np.zeros((n, n), dtype=float)
for i in range(n - 1):
    gi = G[:, i]
    for j in range(i + 1, n):
        gj = G[:, j]
        ok = np.isfinite(gi) & np.isfinite(gj)
        if np.any(ok):
            dij = np.mean(np.abs(gi[ok] - gj[ok]) / 2.0)
        else:
            dij = np.nan
        gen[i, j] = gen[j, i] = dij

        dgeo = haversine_km(meta2.loc[i, 'Latitude'], meta2.loc[i, 'Longitude'],
                            meta2.loc[j, 'Latitude'], meta2.loc[j, 'Longitude'])
        geo[i, j] = geo[j, i] = dgeo

# impute any nan genetic distances to median
mgen = np.nanmedian(gen[np.isfinite(gen)])
gen[~np.isfinite(gen)] = mgen

triu = np.triu_indices(n, k=1)
x = geo[triu]
y = gen[triu]
pear_r, pear_p = stats.pearsonr(x, y)
spear_r, spear_p = stats.spearmanr(x, y)
mantel_r, mantel_p = mantel_test(gen, geo, perms=1999)

ibd_stats = pd.DataFrame([{
    'N_samples': n,
    'N_pairwise_comparisons': len(x),
    'Pearson_r': pear_r,
    'Pearson_p': pear_p,
    'Spearman_rho': spear_r,
    'Spearman_p': spear_p,
    'Mantel_r': mantel_r,
    'Mantel_p_permutation': mantel_p,
    'Permutations': 1999,
}])
ibd_stats.to_csv(os.path.join(TAB, 'IBD_statistics_summary.tsv'), sep='\t', index=False)

pairs = []
for i in range(n - 1):
    for j in range(i + 1, n):
        pairs.append({
            'Sample1': samples[i],
            'Sample2': samples[j],
            'Population1': meta2.loc[i, 'Population'],
            'Population2': meta2.loc[j, 'Population'],
            'Geographic_distance_km': geo[i, j],
            'Genetic_distance': gen[i, j],
        })
pairs_df = pd.DataFrame(pairs)
pairs_df.to_csv(os.path.join(TAB, 'IBD_pairwise_distances.tsv'), sep='\t', index=False)

# IBD figure: scatter + trend
for ext in ('png', 'pdf'):
    plt.figure(figsize=(8.4, 6.8), dpi=320 if ext == 'png' else None)
    ax = plt.gca()
    ax.scatter(x, y, s=14, c='#4c78a8', alpha=0.45, edgecolors='none')

    # robust simple linear fit for visual trend
    slope, intercept, r, p, se = stats.linregress(x, y)
    xs = np.linspace(np.nanmin(x), np.nanmax(x), 300)
    ax.plot(xs, intercept + slope * xs, color='#d62728', linewidth=2.0)

    txt = (
        f"Pearson r = {pear_r:.3f} (p={pear_p:.2e})\\n"
        f"Spearman rho = {spear_r:.3f} (p={spear_p:.2e})\\n"
        f"Mantel r = {mantel_r:.3f} (p={mantel_p:.4f}, 1999 perms)"
    )
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va='top', ha='left', fontsize=9,
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='#cccccc'))

    ax.set_xlabel('Geographic distance (km)')
    ax.set_ylabel('Genetic distance (mean allele-dosage difference/2)')
    ax.set_title('Isolation by distance (individual-level)')
    ax.grid(color='#E3E3E3', linestyle=':', linewidth=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f'IBD_scatter_regression_publication.{ext}'))
    plt.close()

# IBD binned summary figure/table
bins = np.quantile(x, np.linspace(0, 1, 11))
# ensure strict monotonic bins
bins = np.unique(bins)
idx = np.digitize(x, bins[1:-1], right=True)
rows = []
for b in range(len(bins) - 1):
    m = idx == b
    if np.sum(m) == 0:
        continue
    rows.append({
        'Bin': b + 1,
        'Geo_km_min': bins[b],
        'Geo_km_max': bins[b + 1],
        'N_pairs': int(np.sum(m)),
        'Mean_genetic_distance': float(np.mean(y[m])),
        'SE_genetic_distance': float(np.std(y[m], ddof=1) / np.sqrt(np.sum(m))) if np.sum(m) > 1 else np.nan,
    })
ibd_bin = pd.DataFrame(rows)
ibd_bin.to_csv(os.path.join(TAB, 'IBD_distance_bin_summary.tsv'), sep='\t', index=False)

for ext in ('png', 'pdf'):
    plt.figure(figsize=(8.0, 5.6), dpi=320 if ext == 'png' else None)
    xc = (ibd_bin['Geo_km_min'] + ibd_bin['Geo_km_max']) / 2.0
    plt.errorbar(xc, ibd_bin['Mean_genetic_distance'], yerr=ibd_bin['SE_genetic_distance'],
                 fmt='o-', color='#1f77b4', ecolor='#1f77b4', elinewidth=1.2, capsize=2.5)
    plt.xlabel('Geographic distance bin midpoint (km)')
    plt.ylabel('Mean genetic distance')
    plt.title('IBD distance-class trend')
    plt.grid(color='#E3E3E3', linestyle=':', linewidth=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f'IBD_distance_bin_trend_publication.{ext}'))
    plt.close()

with open(os.path.join(TAB, 'maps_and_ibd_analysis_summary.txt'), 'w') as out:
    out.write('Maps and IBD analysis summary\n')
    out.write(f'Samples with coordinates and popmap match: {len(meta)}\n')
    out.write('Figures:\n')
    out.write('- Map_sample_locations_by_population_publication.pdf/png\n')
    out.write('- Map_admixture_pies_by_sample_publication.pdf/png\n')
    out.write('- IBD_scatter_regression_publication.pdf/png\n')
    out.write('- IBD_distance_bin_trend_publication.pdf/png\n')
    out.write('Tables:\n')
    out.write('- sample_coordinates_with_population.tsv\n')
    out.write('- IBD_statistics_summary.tsv\n')
    out.write('- IBD_pairwise_distances.tsv\n')
    out.write('- IBD_distance_bin_summary.tsv\n')

print('Done', flush=True)
