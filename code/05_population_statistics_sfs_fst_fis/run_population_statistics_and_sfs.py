#!/usr/bin/env python3
import os
import math
import itertools
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import dadi
import dadi.Demographics2D as D2

BASE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
OUT = os.path.join(BASE, "downstream_analysis_2026-02-23", "additional_analyses_2026-02-24")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
LOG = os.path.join(OUT, "logs")
WORK = os.path.join(OUT, "work")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)
os.makedirs(LOG, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

VCF = os.path.join(BASE, "populations.snps.vcf")
POPMAP = os.path.join(BASE, "popmap.tsv")

pop_cols = {
    "North": "#4169E1",
    "North_Coast": "#4CBB17",
    "Sierra_1": "#800080",
    "Sierra_2": "#00FFFF",
    "Sierra_3": "#FF0000",
    "South_Coast": "#000000",
}
pop_order = ["North", "North_Coast", "Sierra_1", "Sierra_2", "Sierra_3", "South_Coast"]


def log(msg):
    print(msg, flush=True)


def parse_gt_to_dosage(gt):
    if gt is None:
        return np.nan
    g = gt.split(':', 1)[0].replace('|', '/')
    if g in ("./.", ".", ".//."):
        return np.nan
    if g == "0/0":
        return 0.0
    if g in ("0/1", "1/0"):
        return 1.0
    if g == "1/1":
        return 2.0
    return np.nan


def read_vcf_matrix(vcf_path, sample_set):
    chrom = []
    pos = []
    vid = []
    samples = None
    geno_rows = []
    with open(vcf_path, 'r') as f:
        for ln in f:
            if ln.startswith('##'):
                continue
            if ln.startswith('#CHROM'):
                parts = ln.rstrip('\n').split('\t')
                all_samples = parts[9:]
                samples = [s for s in all_samples if s in sample_set]
                sample_idx = [all_samples.index(s) for s in samples]
                continue
            parts = ln.rstrip('\n').split('\t')
            ref, alt = parts[3], parts[4]
            if len(ref) != 1 or len(alt) != 1 or ',' in alt:
                continue
            row = []
            gts = parts[9:]
            for j in sample_idx:
                row.append(parse_gt_to_dosage(gts[j]))
            chrom.append(parts[0])
            pos.append(int(parts[1]))
            vid.append(parts[2] if parts[2] != '.' else f"{parts[0]}:{parts[1]}")
            geno_rows.append(row)
    G = np.array(geno_rows, dtype=np.float64)
    return samples, np.array(chrom), np.array(pos), np.array(vid), G


# ---------- Load data ----------
log("Loading popmap")
popmap = pd.read_csv(POPMAP, sep='\t', header=None, names=['Sample', 'Population'])
popmap = popmap[popmap['Population'].isin(pop_order)].copy()
sample_set = set(popmap['Sample'])

log("Parsing VCF into genotype matrix")
samples, chrom, pos, vid, G = read_vcf_matrix(VCF, sample_set)

popmap = popmap.set_index('Sample').loc[samples].reset_index()
sample_to_pop = dict(zip(popmap['Sample'], popmap['Population']))

pop_idx = {p: np.array([i for i, s in enumerate(samples) if sample_to_pop[s] == p], dtype=int)
           for p in pop_order}
pop_n = {p: len(pop_idx[p]) for p in pop_order}

# base biallelic filters
site_missing = np.mean(np.isnan(G), axis=1)
alt_mean = np.nanmean(G, axis=1) / 2.0
maf = np.minimum(alt_mean, 1 - alt_mean)

mask_20_005 = (site_missing <= 0.20) & np.isfinite(maf) & (maf >= 0.05)
Gf = G[mask_20_005]
chrom_f = chrom[mask_20_005]
pos_f = pos[mask_20_005]
vid_f = vid[mask_20_005]

log(f"Retained SNPs for main analyses (missing<=0.2, MAF>=0.05): {Gf.shape[0]}")

# per-pop counts on filtered set
alt_counts = {}
call_counts = {}
for p in pop_order:
    idx = pop_idx[p]
    sub = Gf[:, idx]
    alt_counts[p] = np.nansum(sub, axis=1)
    call_counts[p] = np.sum(~np.isnan(sub), axis=1)

# ---------- (1) Formal demographic model fitting with dadi ----------
log("Running demographic model fitting (dadi) on pairwise folded 2D SFS")
models = {
    "SI_no_mig": {
        "func": D2.no_mig,
        "p0": np.array([1.0, 1.0, 0.5]),
        "lb": [1e-3, 1e-3, 1e-4],
        "ub": [50, 50, 10],
    },
    "IM_sym_mig": {
        "func": D2.sym_mig,
        "p0": np.array([1.0, 1.0, 0.5, 0.5]),
        "lb": [1e-3, 1e-3, 1e-4, 1e-5],
        "ub": [50, 50, 10, 30],
    },
    "SC_sym_mig": {
        "func": D2.sec_contact_sym_mig,
        "p0": np.array([1.0, 1.0, 0.2, 0.2, 0.5]),
        "lb": [1e-3, 1e-3, 1e-5, 1e-4, 1e-4],
        "ub": [50, 50, 30, 10, 10],
    },
}

pair_rows = []
demog_comp = []
if os.environ.get("DADI_PAIR_MODE", "all").lower() == "core":
    pair_iter = [
        ("North", "North_Coast"),
        ("North_Coast", "Sierra_1"),
        ("Sierra_1", "Sierra_2"),
        ("Sierra_2", "Sierra_3"),
        ("Sierra_3", "South_Coast"),
        ("North", "Sierra_3"),
        ("North", "South_Coast"),
        ("North_Coast", "South_Coast"),
    ]
else:
    pair_iter = list(itertools.combinations(pop_order, 2))

for p1, p2 in pair_iter:
    log(f"  dadi pair: {p1} vs {p2}")
    i1 = pop_idx[p1]
    i2 = pop_idx[p2]
    n1 = len(i1)
    n2 = len(i2)
    # require complete calls within both populations for direct fixed-size SFS
    m = (call_counts[p1] == n1) & (call_counts[p2] == n2)
    if np.sum(m) < 500:
        continue

    a1 = alt_counts[p1][m].astype(int)
    a2 = alt_counts[p2][m].astype(int)
    nchr1 = 2 * n1
    nchr2 = 2 * n2

    # fold by minor allele in the pair
    tot_alt = a1 + a2
    flip = tot_alt > (nchr1 + nchr2) / 2
    a1 = np.where(flip, nchr1 - a1, a1)
    a2 = np.where(flip, nchr2 - a2, a2)

    fs_arr = np.zeros((nchr1 + 1, nchr2 + 1), dtype=float)
    for x, y in zip(a1, a2):
        fs_arr[x, y] += 1.0

    fs = dadi.Spectrum(fs_arr).fold()
    ns = fs.sample_sizes
    pts = [max(ns) + 10, max(ns) + 20, max(ns) + 30]

    local = []
    for mname, spec in models.items():
        func_ex = dadi.Numerics.make_extrap_log_func(spec["func"])
        try:
            # Fast randomized composite-likelihood search over parameter space.
            best_ll = -np.inf
            best_p = None
            for _ in range(8):
                draw = []
                for lo, hi in zip(spec["lb"], spec["ub"]):
                    lo2 = max(lo, 1e-8)
                    draw.append(np.exp(np.random.uniform(np.log(lo2), np.log(hi))))
                ptry = np.array(draw, dtype=float)
                mfs = func_ex(ptry, ns, pts)
                ll_try = dadi.Inference.ll_multinom(mfs, fs)
                if np.isfinite(ll_try) and ll_try > best_ll:
                    best_ll = ll_try
                    best_p = ptry
            if best_p is None:
                continue
            model_fs = func_ex(best_p, ns, pts)
            ll = dadi.Inference.ll_multinom(model_fs, fs)
            theta = dadi.Inference.optimal_sfs_scaling(model_fs, fs)
            k = len(best_p)
            aic = 2 * k - 2 * ll
            local.append((mname, ll, aic, theta, best_p))
        except Exception:
            continue

    if len(local) == 0:
        continue

    local = sorted(local, key=lambda z: z[2])
    best_aic = local[0][2]
    for mname, ll, aic, theta, popt in local:
        demog_comp.append({
            "Pop1": p1,
            "Pop2": p2,
            "Model": mname,
            "Sites_used": int(fs.S()),
            "LogLik": float(ll),
            "AIC": float(aic),
            "DeltaAIC": float(aic - best_aic),
            "Theta": float(theta),
            "Params": ",".join([f"{x:.4g}" for x in np.ravel(popt)]),
        })

    pair_rows.append({
        "Pop1": p1,
        "Pop2": p2,
        "Best_model": local[0][0],
        "Best_AIC": local[0][2],
        "Second_AIC": local[1][2] if len(local) > 1 else np.nan,
        "DeltaAIC_second": (local[1][2] - local[0][2]) if len(local) > 1 else np.nan,
        "Sites_used": int(fs.S()),
    })

pd.DataFrame(demog_comp).to_csv(os.path.join(TAB, "demographic_model_comparison_dadi.tsv"), sep='\t', index=False)
demog_best = pd.DataFrame(pair_rows)
demog_best.to_csv(os.path.join(TAB, "demographic_model_best_by_pair.tsv"), sep='\t', index=False)

# figure for analysis 1
if len(demog_best) > 0:
    mat = pd.DataFrame(index=pop_order, columns=pop_order, data="")
    for _, r in demog_best.iterrows():
        mat.loc[r['Pop1'], r['Pop2']] = r['Best_model']
        mat.loc[r['Pop2'], r['Pop1']] = r['Best_model']
    for p in pop_order:
        mat.loc[p, p] = "-"

    model_color = {
        "SI_no_mig": "#4daf4a",
        "IM_sym_mig": "#377eb8",
        "SC_sym_mig": "#e41a1c",
        "-": "#f0f0f0",
        "": "#ffffff",
    }
    cmat = np.zeros((len(pop_order), len(pop_order), 3))
    import matplotlib.colors as mcolors
    for i, a in enumerate(pop_order):
        for j, b in enumerate(pop_order):
            c = model_color.get(mat.loc[a, b], "#ffffff")
            cmat[i, j, :] = mcolors.to_rgb(c)

    for ext in ["pdf", "png"]:
        plt.figure(figsize=(8, 7), dpi=320 if ext == "png" else None)
        plt.imshow(cmat, aspect='equal')
        plt.xticks(range(len(pop_order)), pop_order, rotation=45, ha='right')
        plt.yticks(range(len(pop_order)), pop_order)
        for i, a in enumerate(pop_order):
            for j, b in enumerate(pop_order):
                txt = mat.loc[a, b]
                if txt == "SI_no_mig":
                    s = "SI"
                elif txt == "IM_sym_mig":
                    s = "IM"
                elif txt == "SC_sym_mig":
                    s = "SC"
                else:
                    s = txt
                plt.text(j, i, s, ha='center', va='center', fontsize=8)
        plt.title("Best-supported pairwise demographic model (dadi 2D-SFS)")
        handles = [plt.Line2D([0],[0], marker='s', markersize=10, linestyle='none', color=c, label=l)
                   for l, c in [("SI", "#4daf4a"), ("IM", "#377eb8"), ("SC", "#e41a1c")]]
        plt.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1.0), frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, f"Demography_dadi_best_model_matrix_publication.{ext}"))
        plt.close()

# ---------- (2) ABBA-BABA / f4 introgression tests ----------
log("Running ABBA-BABA/f4 introgression tests")
outgroup = "South_Coast"
in_pops = [p for p in pop_order if p != outgroup]

rows = []
for p1, p2, p3 in itertools.permutations(in_pops, 3):
    if len({p1, p2, p3}) < 3:
        continue

    # unique orientation only for p1<p2 lexicographically to avoid duplicates
    if p1 > p2:
        continue

    c1 = alt_counts[p1]
    c2 = alt_counts[p2]
    c3 = alt_counts[p3]
    c4 = alt_counts[outgroup]
    n1 = 2 * call_counts[p1]
    n2 = 2 * call_counts[p2]
    n3 = 2 * call_counts[p3]
    n4 = 2 * call_counts[outgroup]

    ok = (n1 > 0) & (n2 > 0) & (n3 > 0) & (n4 > 0)
    if np.sum(ok) < 500:
        continue

    p1f = np.clip(c1[ok] / n1[ok], 0, 1)
    p2f = np.clip(c2[ok] / n2[ok], 0, 1)
    p3f = np.clip(c3[ok] / n3[ok], 0, 1)
    p4f = np.clip(c4[ok] / n4[ok], 0, 1)
    ch = chrom_f[ok]

    abba = (1 - p1f) * p2f * p3f * (1 - p4f)
    baba = p1f * (1 - p2f) * p3f * (1 - p4f)

    denom = np.sum(abba + baba)
    if denom <= 0:
        continue
    D = np.sum(abba - baba) / denom

    # block-jackknife by chromosome/scaffold
    uniq = np.unique(ch)
    if len(uniq) < 5:
        continue
    d_leave = []
    for u in uniq:
        keep = ch != u
        dnm = np.sum((abba[keep] + baba[keep]))
        if dnm <= 0:
            continue
        d_leave.append(np.sum(abba[keep] - baba[keep]) / dnm)
    d_leave = np.array(d_leave)
    if len(d_leave) < 5:
        continue
    nblk = len(d_leave)
    dbar = np.mean(d_leave)
    se = np.sqrt((nblk - 1) / nblk * np.sum((d_leave - dbar) ** 2))
    z = D / se if se > 0 else np.nan
    pval = 2 * stats.norm.sf(abs(z)) if np.isfinite(z) else np.nan

    f4 = np.mean((p1f - p2f) * (p3f - p4f))

    rows.append({
        "P1": p1,
        "P2": p2,
        "P3": p3,
        "Outgroup": outgroup,
        "Sites_used": int(np.sum(ok)),
        "D": float(D),
        "SE": float(se),
        "Z": float(z),
        "P": float(pval),
        "f4": float(f4),
    })

abba = pd.DataFrame(rows).drop_duplicates(subset=["P1", "P2", "P3", "Outgroup"])
if len(abba):
    abba['FDR_BH'] = np.minimum(1.0, abba['P'] * len(abba) / np.maximum(1, abba['P'].rank(method='first')))
    abba = abba.sort_values('P')
abba.to_csv(os.path.join(TAB, "introgression_ABBA_BABA_f4_tests.tsv"), sep='\t', index=False)

if len(abba) > 0:
    for ext in ["pdf", "png"]:
        plt.figure(figsize=(9, 6), dpi=320 if ext == "png" else None)
        x = abba['D'].values
        y = -np.log10(np.clip(abba['P'].values, 1e-300, 1))
        sig = (abba['P'] < 0.05).values
        plt.scatter(x[~sig], y[~sig], c='#999999', s=35, alpha=0.8)
        plt.scatter(x[sig], y[sig], c='#d62728', s=42, alpha=0.9)
        for _, r in abba.head(10).iterrows():
            plt.text(r['D'], -np.log10(max(r['P'], 1e-300)), f"{r['P1']},{r['P2']}|{r['P3']}", fontsize=7)
        plt.axvline(0, color='black', lw=1, ls='--')
        plt.axhline(-np.log10(0.05), color='black', lw=1, ls=':')
        plt.xlabel("D statistic")
        plt.ylabel("-log10(p)")
        plt.title(f"ABBA-BABA introgression tests (outgroup={outgroup})")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, f"ABBA_BABA_introgression_volcano_publication.{ext}"))
        plt.close()

# ---------- (3) Selection scan ----------
log("Running genome scan for putative selected loci")
scan_rows = []
for i in range(Gf.shape[0]):
    ac = []
    nc = []
    pvec = []
    for p in pop_order:
        n = 2 * call_counts[p][i]
        a = alt_counts[p][i]
        ac.append(a)
        nc.append(n)
        pvec.append(a / n if n > 0 else np.nan)

    ac = np.array(ac, dtype=float)
    nc = np.array(nc, dtype=float)
    pvec = np.array(pvec, dtype=float)
    valid = nc > 0
    if np.sum(valid) < 2:
        continue

    # global Fst (Ht-Hs)/Ht
    w = nc[valid]
    pbar = np.sum(ac[valid]) / np.sum(nc[valid])
    hs = np.average(2 * pvec[valid] * (1 - pvec[valid]), weights=w)
    ht = 2 * pbar * (1 - pbar)
    fst = (ht - hs) / ht if ht > 0 else np.nan

    # heterogeneity chi-square on allele counts across populations
    ref = nc[valid] - ac[valid]
    table = np.vstack([ref, ac[valid]])
    try:
        chi2, pval, _, _ = stats.chi2_contingency(table)
    except Exception:
        pval = np.nan

    he = hs
    scan_rows.append({
        "CHROM": chrom_f[i],
        "POS": int(pos_f[i]),
        "ID": vid_f[i],
        "Fst_global": float(fst),
        "He_mean": float(he),
        "P_freq_diff": float(pval),
    })

scan = pd.DataFrame(scan_rows)
scan = scan[np.isfinite(scan['Fst_global']) & np.isfinite(scan['P_freq_diff'])].copy()
scan = scan.sort_values('P_freq_diff')
rank = np.arange(1, len(scan) + 1)
scan['FDR_BH'] = np.minimum(1.0, scan['P_freq_diff'] * len(scan) / rank)
fst_thr = np.quantile(scan['Fst_global'], 0.99)
scan['Outlier'] = (scan['Fst_global'] >= fst_thr) & (scan['FDR_BH'] <= 0.05)
scan.to_csv(os.path.join(TAB, "selection_scan_snp_stats.tsv"), sep='\t', index=False)
scan.sort_values(['Outlier', 'Fst_global', 'P_freq_diff'], ascending=[False, False, True]).head(200).to_csv(
    os.path.join(TAB, "selection_scan_top_candidates.tsv"), sep='\t', index=False
)

# Manhattan-like position index
chrom_order = pd.unique(scan['CHROM'])
offset = {}
cum = 0
for c in chrom_order:
    m = scan.loc[scan['CHROM'] == c, 'POS'].max()
    offset[c] = cum
    cum += m + 1
scan['POS_CUM'] = scan['POS'] + scan['CHROM'].map(offset)

for ext in ["pdf", "png"]:
    # Manhattan Fst
    plt.figure(figsize=(12, 5), dpi=320 if ext == "png" else None)
    colors = ['#4c78a8', '#9ecae9']
    for ii, c in enumerate(chrom_order):
        d = scan[scan['CHROM'] == c]
        plt.scatter(d['POS_CUM'], d['Fst_global'], s=8, c=colors[ii % 2], alpha=0.7, linewidths=0)
    out = scan[scan['Outlier']]
    plt.scatter(out['POS_CUM'], out['Fst_global'], s=16, c='#d62728', alpha=0.95, linewidths=0)
    plt.axhline(fst_thr, ls='--', lw=1, color='black')
    plt.ylabel('Global Fst')
    plt.xlabel('Genomic position (concatenated)')
    plt.title('Selection scan: global Fst outliers')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f"Selection_scan_Fst_manhattan_publication.{ext}"))
    plt.close()

    # Fst vs heterozygosity
    plt.figure(figsize=(7, 6), dpi=320 if ext == "png" else None)
    plt.scatter(scan['He_mean'], scan['Fst_global'], s=10, c='#999999', alpha=0.5)
    plt.scatter(out['He_mean'], out['Fst_global'], s=18, c='#d62728', alpha=0.9)
    plt.xlabel('Mean expected heterozygosity (He)')
    plt.ylabel('Global Fst')
    plt.title('Selection scan: Fst vs heterozygosity')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f"Selection_scan_Fst_vs_He_publication.{ext}"))
    plt.close()

# ---------- (6) Kinship/inbreeding ----------
log("Computing kinship and individual inbreeding summaries")
# matrix for relatedness
p = np.nanmean(Gf, axis=1) / 2.0
ok = np.isfinite(p) & (p > 0) & (p < 1)
X = Gf[ok].copy()
p = p[ok]

# inbreeding (per individual)
ind_rows = []
for j, s in enumerate(samples):
    gj = Gf[:, j]
    called = np.isfinite(gj)
    if np.sum(called) < 100:
        continue
    ho = np.mean(gj[called] == 1)
    pj = np.nanmean(Gf[called], axis=1) / 2.0
    he = np.mean(2 * pj * (1 - pj))
    fhat = 1 - ho / he if he > 0 else np.nan
    ind_rows.append({
        "Sample": s,
        "Population": sample_to_pop[s],
        "N_called": int(np.sum(called)),
        "Ho": float(ho),
        "He": float(he),
        "Fhat_inbreeding": float(fhat),
    })
ind_df = pd.DataFrame(ind_rows)
ind_df.to_csv(os.path.join(TAB, "individual_inbreeding_summary.tsv"), sep='\t', index=False)

# GRM / kinship
X_imp = X.copy()
for i in range(X_imp.shape[0]):
    miss = ~np.isfinite(X_imp[i])
    if np.any(miss):
        X_imp[i, miss] = 2 * p[i]

den = np.sqrt(2 * p * (1 - p))
Z = (X_imp - (2 * p[:, None])) / den[:, None]
K = (Z.T @ Z) / Z.shape[0] / 2.0  # kinship-like

kin_df = pd.DataFrame(K, index=samples, columns=samples)
kin_df.to_csv(os.path.join(TAB, "pairwise_kinship_matrix.tsv"), sep='\t')

pairs = []
for i in range(len(samples) - 1):
    for j in range(i + 1, len(samples)):
        pairs.append({
            "Sample1": samples[i],
            "Sample2": samples[j],
            "Population1": sample_to_pop[samples[i]],
            "Population2": sample_to_pop[samples[j]],
            "Kinship": float(K[i, j]),
        })
pairs_df = pd.DataFrame(pairs).sort_values('Kinship', ascending=False)
pairs_df.head(500).to_csv(os.path.join(TAB, "pairwise_kinship_top500.tsv"), sep='\t', index=False)

# figures for kinship/inbreeding
ord_samples = sorted(samples, key=lambda s: (pop_order.index(sample_to_pop[s]), s))
ord_idx = [samples.index(s) for s in ord_samples]
Kord = K[np.ix_(ord_idx, ord_idx)]

for ext in ["pdf", "png"]:
    plt.figure(figsize=(9, 8), dpi=320 if ext == "png" else None)
    im = plt.imshow(Kord, cmap='viridis', aspect='auto')
    plt.colorbar(im, fraction=0.046, pad=0.04, label='Kinship')
    plt.title('Pairwise kinship matrix (samples ordered by population)')
    plt.xticks([])
    plt.yticks([])
    # population separators
    breaks = []
    cur = 0
    for pnm in pop_order:
        n = sum(sample_to_pop[s] == pnm for s in ord_samples)
        cur += n
        breaks.append(cur)
    for b in breaks[:-1]:
        plt.axhline(b - 0.5, color='white', lw=0.6)
        plt.axvline(b - 0.5, color='white', lw=0.6)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f"Kinship_matrix_publication.{ext}"))
    plt.close()

    plt.figure(figsize=(8, 5), dpi=320 if ext == "png" else None)
    dfp = ind_df[ind_df['Population'].isin(pop_order)].copy()
    data = [dfp.loc[dfp['Population'] == p, 'Fhat_inbreeding'].dropna().values for p in pop_order]
    bp = plt.boxplot(data, patch_artist=True, labels=pop_order, showfliers=False)
    for patch, pnm in zip(bp['boxes'], pop_order):
        patch.set_facecolor(pop_cols[pnm])
        patch.set_alpha(0.7)
    plt.axhline(0, color='black', ls='--', lw=1)
    plt.ylabel('Individual inbreeding coefficient (Fhat)')
    plt.xticks(rotation=35, ha='right')
    plt.title('Inbreeding by population')
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, f"Inbreeding_by_population_publication.{ext}"))
    plt.close()

# summary text
with open(os.path.join(TAB, "additional_analyses_1_2_3_6_summary.txt"), 'w') as out:
    out.write("Additional analyses completed (1,2,3,6)\n")
    out.write(f"Input SNPs after filter (missing<=0.20, MAF>=0.05): {Gf.shape[0]}\n")
    out.write("1) Demographic model fitting: dadi pairwise 2D folded SFS; models SI, IM, SC\n")
    out.write("2) Introgression tests: ABBA-BABA D and f4; block-jackknife by scaffold\n")
    out.write("3) Selection scan: global Fst + allele-frequency heterogeneity test (BH FDR)\n")
    out.write("6) Kinship/inbreeding: GRM-derived kinship + individual Fhat\n")

log("Done")
