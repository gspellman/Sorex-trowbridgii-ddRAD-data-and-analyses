#!/usr/bin/env python3
import os
import sys
import math
import time
import random
from collections import defaultdict

import numpy as np
import pandas as pd
try:
    import matplotlib.pyplot as plt
except Exception:
    plt = None


BASE_DIR = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
DOWN_DIR = os.path.join(BASE_DIR, "downstream_analysis_2026-02-23")
OUT_DIR = os.path.join(DOWN_DIR, "species_delimitation_speede_rf_delimitr_2026-02-26")
TAB_DIR = os.path.join(OUT_DIR, "tables")
FIG_DIR = os.path.join(OUT_DIR, "figures")
LOG_DIR = os.path.join(OUT_DIR, "logs")

VCF_FILE = os.path.join(BASE_DIR, "populations.snps.vcf")
POPMAP_FILE = os.path.join(BASE_DIR, "popmap.tsv")
IQTREE_GUIDE_FILE = os.path.join(
    DOWN_DIR, "iqtree_analysis_2026-03-01", "trow_individuals_biallelic_iqtree.contree"
)

os.makedirs(TAB_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

RNG = np.random.default_rng(20260226)
random.seed(20260226)


class NewickNode:
    __slots__ = ("name", "length", "support", "children", "parent", "leafset")

    def __init__(self):
        self.name = None
        self.length = None
        self.support = None
        self.children = []
        self.parent = None
        self.leafset = None


def _parse_label_and_length(s, i):
    name = None
    length = None
    label = []
    n = len(s)
    while i < n and s[i] not in [":", ",", ")", ";"]:
        label.append(s[i])
        i += 1
    label_str = "".join(label).strip()
    if label_str:
        name = label_str
    if i < n and s[i] == ":":
        i += 1
        ln = []
        while i < n and s[i] not in [",", ")", ";"]:
            ln.append(s[i])
            i += 1
        try:
            length = float("".join(ln))
        except Exception:
            length = None
    return name, length, i


def parse_newick_tree(newick_text):
    s = "".join([c for c in newick_text.strip() if c not in "\n\r\t "])
    if not s.endswith(";"):
        s = s + ";"
    i = 0
    n = len(s)

    def parse_subtree(idx):
        node = NewickNode()
        if s[idx] == "(":
            idx += 1
            while True:
                child, idx = parse_subtree(idx)
                child.parent = node
                node.children.append(child)
                if idx < n and s[idx] == ",":
                    idx += 1
                    continue
                if idx < n and s[idx] == ")":
                    idx += 1
                    break
            label, length, idx = _parse_label_and_length(s, idx)
            if label is not None:
                try:
                    node.support = float(label)
                except Exception:
                    node.name = label
            node.length = length
        else:
            name, length, idx = _parse_label_and_length(s, idx)
            node.name = name
            node.length = length
        return node, idx

    root, i = parse_subtree(i)
    return root


def index_tree(root):
    leaves = {}
    nodes = []

    def postorder(node):
        nodes.append(node)
        if not node.children:
            node.leafset = set([node.name]) if node.name else set()
            if node.name:
                leaves[node.name] = node
            return node.leafset
        ls = set()
        for ch in node.children:
            ls |= postorder(ch)
        node.leafset = ls
        return ls

    postorder(root)
    return leaves, nodes


def mrca_of_taxa(leaves_index, taxa):
    taxa = [t for t in taxa if t in leaves_index]
    if not taxa:
        return None
    anc = []
    x = leaves_index[taxa[0]]
    while x is not None:
        anc.append(x)
        x = x.parent
    anc_set = set(anc)
    for t in taxa[1:]:
        y = leaves_index[t]
        path = []
        while y is not None:
            path.append(y)
            y = y.parent
        common = None
        for z in path:
            if z in anc_set:
                common = z
                break
        if common is None:
            return None
        anc = anc[: anc.index(common) + 1]
        anc_set = set(anc)
    return anc[-1] if anc else None


def monophyly_score_for_labels(samples, labels, leaves_index):
    df = pd.DataFrame({"Sample": samples, "Group": labels})
    by_group = df.groupby("Group")["Sample"].apply(list).to_dict()
    weights = []
    penalties = []
    missing = 0
    for grp, taxa in by_group.items():
        taxa_present = [t for t in taxa if t in leaves_index]
        missing += (len(taxa) - len(taxa_present))
        if len(taxa_present) <= 1:
            purity = 1.0
        else:
            mrca = mrca_of_taxa(leaves_index, taxa_present)
            if mrca is None or mrca.leafset is None or len(mrca.leafset) == 0:
                purity = 0.0
            else:
                purity = len(set(taxa_present)) / float(len(mrca.leafset))
        weights.append(len(taxa))
        penalties.append(1.0 - purity)
    if not weights:
        return 0.0, 1.0
    nonmono = float(np.average(np.array(penalties), weights=np.array(weights)))
    monophy = max(0.0, 1.0 - nonmono)
    if missing > 0:
        frac_missing = missing / float(len(samples))
        monophy = max(0.0, monophy - frac_missing)
    return nonmono, monophy


def log(msg):
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(os.path.join(LOG_DIR, "run.log"), "a") as fh:
        fh.write(line + "\n")


def parse_gt(gt):
    gt = gt.split(":", 1)[0].replace("|", "/")
    if gt in ("./.", ".", ".|."):
        return -1
    if gt == "0/0":
        return 0
    if gt in ("0/1", "1/0"):
        return 1
    if gt == "1/1":
        return 2
    return -1


def load_popmap(popmap_file):
    pop = pd.read_csv(popmap_file, sep="\t", header=0)
    pop.columns = ["Sample", "Population"]
    return pop


def read_unlinked_snps(vcf_file, sample_order, callrate_min=0.80, maf_min=0.01):
    log("Reading VCF and selecting one SNP per locus (unlinked dataset)")
    locus_best = {}
    kept = 0
    with open(vcf_file, "r") as fh:
        samples = None
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                continue
            if line.startswith("#"):
                continue
            p = line.rstrip().split("\t")
            if len(p) < 10:
                continue
            chrom, pos, vid, ref, alt = p[0], int(p[1]), p[2], p[3], p[4]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            locus_id = vid.split(":")[0]
            locus_key = f"{chrom}:{locus_id}"
            g = np.array([parse_gt(x) for x in p[9:]], dtype=np.int16)
            called = np.sum(g >= 0)
            if called == 0:
                continue
            alt_count = np.sum(g[g >= 0])
            af = alt_count / (2.0 * called)
            maf = min(af, 1.0 - af)
            callrate = called / g.size
            if callrate < callrate_min or maf < maf_min:
                continue
            rank = (callrate, maf)
            if locus_key not in locus_best or rank > locus_best[locus_key]["rank"]:
                locus_best[locus_key] = {
                    "chrom": chrom,
                    "pos": pos,
                    "id": vid,
                    "rank": rank,
                    "g": g,
                }
                kept += 1

    if samples is None:
        raise RuntimeError("VCF header with samples not found")

    # Reorder to popmap sample order
    idx = [samples.index(s) for s in sample_order]
    rows = []
    meta = []
    for lk, rec in locus_best.items():
        rows.append(rec["g"][idx])
        meta.append((lk, rec["chrom"], rec["pos"], rec["id"], rec["rank"][0], rec["rank"][1]))
    X = np.array(rows, dtype=np.int16)
    mdf = pd.DataFrame(meta, columns=["LocusKey", "CHROM", "POS", "ID", "CallRate", "MAF"])

    log(f"Selected {X.shape[0]} unlinked SNPs after filtering")
    return X, mdf


def impute_by_population_mode(X, pops):
    X2 = X.copy().astype(np.int16)
    pop_levels = sorted(set(pops))
    for j in range(X2.shape[0]):
        col = X2[j, :]
        miss = np.where(col < 0)[0]
        if miss.size == 0:
            continue
        global_nonmiss = col[col >= 0]
        global_mode = int(np.bincount(global_nonmiss).argmax()) if global_nonmiss.size else 0
        for i in miss:
            pp = pops[i]
            idx = np.where(np.array(pops) == pp)[0]
            vals = col[idx]
            vals = vals[vals >= 0]
            if vals.size:
                fill = int(np.bincount(vals).argmax())
            else:
                fill = global_mode
            col[i] = fill
        X2[j, :] = col
    return X2


def pca_scores(X, n_pc=10):
    # X is SNP x sample; convert to sample x SNP
    M = X.T.astype(float)
    M -= M.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    pcs = U[:, :n_pc] * S[:n_pc]
    return pcs


def model_labels(pop):
    p = np.array(pop)
    models = {}
    models["M1_6pop"] = p.copy()
    models["M2_NorthMerge"] = np.where(np.isin(p, ["North", "North_Coast"]), "North_Block", p)
    models["M3_SierraMerge"] = np.where(np.isin(p, ["Sierra_1", "Sierra_2", "Sierra_3"]), "Sierra_Block", p)
    out = []
    for x in p:
        if x in ("North", "North_Coast"):
            out.append("North_Block")
        elif x in ("Sierra_1", "Sierra_2", "Sierra_3"):
            out.append("Sierra_Block")
        else:
            out.append("South_Coast")
    models["M4_ThreeSpecies"] = np.array(out)
    out = []
    for x in p:
        if x in ("North", "North_Coast", "South_Coast"):
            out.append("Coastal")
        else:
            out.append("Inland")
    models["M5_CoastVsInland"] = np.array(out)
    return models


def speede_score(pc, labels):
    n = pc.shape[0]
    labs = np.array(labels)
    groups = sorted(set(labs))
    wss = 0.0
    bss = 0.0
    grand = pc.mean(axis=0)
    for g in groups:
        idx = np.where(labs == g)[0]
        Xg = pc[idx, :]
        mu = Xg.mean(axis=0)
        wss += np.sum((Xg - mu) ** 2)
        bss += len(idx) * np.sum((mu - grand) ** 2)
    k = len(groups)
    bic = n * np.log((wss / max(n, 1)) + 1e-12) + (k - 1) * np.log(max(n, 2))
    sep = bss / (wss + 1e-12)
    return bic, sep, k


class TreeNode:
    __slots__ = ("is_leaf", "pred", "proba", "feat", "thr", "left", "right")
    def __init__(self):
        self.is_leaf = True
        self.pred = None
        self.proba = None
        self.feat = None
        self.thr = None
        self.left = None
        self.right = None


def gini(y, n_class):
    if y.size == 0:
        return 0.0
    cnt = np.bincount(y, minlength=n_class).astype(float)
    p = cnt / cnt.sum()
    return 1.0 - np.sum(p * p)


class SimpleRF:
    def __init__(self, n_trees=80, max_depth=8, min_leaf=3, mtry=None, seed=20260226):
        self.n_trees = n_trees
        self.max_depth = max_depth
        self.min_leaf = min_leaf
        self.mtry = mtry
        self.rng = np.random.default_rng(seed)
        self.trees = []
        self.n_class = None
        self.n_feat = None

    def _best_split(self, X, y):
        n, p = X.shape
        mtry = self.mtry if self.mtry is not None else max(1, int(math.sqrt(p)))
        feats = self.rng.choice(p, size=min(mtry, p), replace=False)
        best = None
        base = gini(y, self.n_class)
        for f in feats:
            vals = X[:, f]
            uq = np.unique(vals)
            if uq.size <= 1:
                continue
            if uq.size > 10:
                qs = np.quantile(uq, [0.1, 0.3, 0.5, 0.7, 0.9])
                thrs = np.unique(qs)
            else:
                thrs = (uq[:-1] + uq[1:]) / 2.0
            for thr in thrs:
                left = vals <= thr
                nL = np.sum(left)
                nR = n - nL
                if nL < self.min_leaf or nR < self.min_leaf:
                    continue
                g = (nL / n) * gini(y[left], self.n_class) + (nR / n) * gini(y[~left], self.n_class)
                gain = base - g
                if (best is None) or (gain > best[0]):
                    best = (gain, f, thr, left)
        return best

    def _build(self, X, y, depth):
        node = TreeNode()
        cnt = np.bincount(y, minlength=self.n_class).astype(float)
        node.pred = int(np.argmax(cnt))
        node.proba = cnt / cnt.sum()
        if depth >= self.max_depth or X.shape[0] < (2 * self.min_leaf) or np.unique(y).size == 1:
            return node
        best = self._best_split(X, y)
        if best is None or best[0] <= 1e-10:
            return node
        _, f, thr, left = best
        node.is_leaf = False
        node.feat = int(f)
        node.thr = float(thr)
        node.left = self._build(X[left], y[left], depth + 1)
        node.right = self._build(X[~left], y[~left], depth + 1)
        return node

    def fit(self, X, y):
        n, p = X.shape
        self.n_feat = p
        self.n_class = int(np.max(y)) + 1
        self.trees = []
        self.boot_idx = []
        for _ in range(self.n_trees):
            idx = self.rng.integers(0, n, size=n)
            self.boot_idx.append(idx)
            tree = self._build(X[idx], y[idx], 0)
            self.trees.append(tree)
        return self

    def _predict_tree_proba(self, tree, X):
        out = np.zeros((X.shape[0], self.n_class), dtype=float)
        for i in range(X.shape[0]):
            node = tree
            while not node.is_leaf:
                if X[i, node.feat] <= node.thr:
                    node = node.left
                else:
                    node = node.right
            out[i, :] = node.proba
        return out

    def predict_proba(self, X):
        acc = np.zeros((X.shape[0], self.n_class), dtype=float)
        for tr in self.trees:
            acc += self._predict_tree_proba(tr, X)
        acc /= len(self.trees)
        return acc

    def oob_proba(self, X):
        n = X.shape[0]
        votes = np.zeros((n, self.n_class), dtype=float)
        counts = np.zeros(n, dtype=int)
        for tr, idx in zip(self.trees, self.boot_idx):
            inbag = np.zeros(n, dtype=bool)
            inbag[idx] = True
            oob = np.where(~inbag)[0]
            if oob.size == 0:
                continue
            pr = self._predict_tree_proba(tr, X[oob])
            votes[oob, :] += pr
            counts[oob] += 1
        counts2 = np.maximum(counts, 1)[:, None]
        out = votes / counts2
        return out, counts


def balanced_accuracy(y_true, y_pred):
    labs = sorted(set(y_true.tolist()))
    rec = []
    for k in labs:
        idx = (y_true == k)
        if np.sum(idx) == 0:
            continue
        rec.append(np.mean(y_pred[idx] == k))
    return float(np.mean(rec)) if rec else np.nan


def grouped_freqs(X, labels):
    labs = np.array(labels)
    groups = sorted(set(labs))
    freq = []
    nind = []
    for g in groups:
        idx = np.where(labs == g)[0]
        G = X[:, idx].astype(float) / 2.0
        p = np.nanmean(G, axis=1)
        freq.append(p)
        nind.append(len(idx))
    return np.array(freq), np.array(nind), groups


def summarize_from_freq(freq, nind):
    # freq: groups x snps
    g, m = freq.shape
    within = np.mean(2 * freq * (1 - freq), axis=1)
    pair_diffs = []
    pair_fst = []
    for i in range(g):
        for j in range(i + 1, g):
            pi = freq[i]
            pj = freq[j]
            d = np.mean(np.abs(pi - pj))
            pair_diffs.append(d)
            hs = (2 * pi * (1 - pi) + 2 * pj * (1 - pj)) / 2.0
            pbar = (pi + pj) / 2.0
            ht = 2 * pbar * (1 - pbar)
            fst = np.mean((ht - hs) / (ht + 1e-12))
            pair_fst.append(fst)
    pair_diffs = np.array(pair_diffs) if pair_diffs else np.array([0.0])
    pair_fst = np.array(pair_fst) if pair_fst else np.array([0.0])
    private = 0.0
    if g > 1:
        for i in range(g):
            pi = freq[i]
            others = np.delete(freq, i, axis=0)
            cond = (pi > 0.95) & (np.max(others, axis=0) < 0.05)
            private += np.mean(cond)
        private /= g
    maf = np.min(np.vstack([freq.mean(axis=0), 1.0 - freq.mean(axis=0)]), axis=0)
    vec = [
        float(g),
        float(np.mean(within)),
        float(np.std(within)),
        float(np.mean(pair_diffs)),
        float(np.std(pair_diffs)),
        float(np.mean(pair_fst)),
        float(np.std(pair_fst)),
        float(private),
        float(np.mean(maf)),
        float(np.std(maf)),
    ]
    # pad groupwise within values to max 6 groups
    within_pad = list(within[:6]) + [0.0] * max(0, 6 - len(within))
    vec.extend(within_pad[:6])
    return np.array(vec, dtype=float)


def simulate_freq(freq_obs, nind, conc=30.0):
    g, m = freq_obs.shape
    out = np.zeros_like(freq_obs)
    for i in range(g):
        p = freq_obs[i]
        a = p * conc + 1.0
        b = (1.0 - p) * conc + 1.0
        p2 = RNG.beta(a, b)
        alt = RNG.binomial(2 * nind[i], p2)
        out[i] = alt / (2.0 * nind[i])
    return out


def main():
    log("Starting SPEEDEMON + RandomForest + DelimitR workflows (offline implementation)")
    pop = load_popmap(POPMAP_FILE)
    samples = pop["Sample"].tolist()
    populations = pop["Population"].tolist()

    X_raw, meta = read_unlinked_snps(VCF_FILE, samples, callrate_min=0.8, maf_min=0.01)
    X_imp = impute_by_population_mode(X_raw, populations)

    # Save unlinked SNP metadata
    meta.to_csv(os.path.join(TAB_DIR, "unlinked_one_snp_per_locus_metadata.tsv"), sep="\t", index=False)
    pd.DataFrame({
        "Metric": ["Samples", "SNPs_unlinked", "Missing_rate_before_impute"],
        "Value": [X_imp.shape[1], X_imp.shape[0], float(np.mean(X_raw < 0))]
    }).to_csv(os.path.join(TAB_DIR, "unlinked_dataset_summary.tsv"), sep="\t", index=False)

    models = model_labels(populations)
    model_names = list(models.keys())
    Xs = X_imp.T.astype(float)  # sample x snp
    guide_nonmono = {}
    guide_mono = {}
    guide_penalty_weight = 0.25

    if os.path.exists(IQTREE_GUIDE_FILE):
        log(f"Loading IQ-TREE guide topology: {IQTREE_GUIDE_FILE}")
        with open(IQTREE_GUIDE_FILE, "r") as fh:
            nwk = fh.read().strip()
        guide_root = parse_newick_tree(nwk)
        guide_leaves, _ = index_tree(guide_root)
        log(f"Guide tree loaded with {len(guide_leaves)} tips")
        for mn in model_names:
            nn, mm = monophyly_score_for_labels(samples, models[mn], guide_leaves)
            guide_nonmono[mn] = nn
            guide_mono[mn] = mm
    else:
        log("IQ-TREE guide topology not found; RF will proceed without guide penalties")
        for mn in model_names:
            guide_nonmono[mn] = 0.0
            guide_mono[mn] = 1.0

    # SPEEDEMON-style scoring on PCA space
    log("Running SPEEDEMON-style model scoring")
    pc = pca_scores(X_imp, n_pc=10)
    s_rows = []
    for mn in model_names:
        bic, sep, k = speede_score(pc, models[mn])
        s_rows.append((mn, k, bic, sep))
    speede_df = pd.DataFrame(s_rows, columns=["Model", "K_groups", "PseudoBIC", "Separation"])
    speede_df = speede_df.sort_values("PseudoBIC").reset_index(drop=True)
    speede_df["DeltaBIC"] = speede_df["PseudoBIC"] - speede_df["PseudoBIC"].min()
    speede_df.to_csv(os.path.join(TAB_DIR, "speedeMON_model_scores.tsv"), sep="\t", index=False)
    best_s = speede_df.iloc[0]["Model"]
    pd.DataFrame({"Sample": samples, "Population": populations, "Delimited_species": models[best_s]}).to_csv(
        os.path.join(TAB_DIR, "speedeMON_best_model_assignments.tsv"), sep="\t", index=False
    )

    if plt is not None:
        plt.figure(figsize=(7.2, 4.8))
        x = np.arange(len(speede_df))
        plt.bar(x, speede_df["PseudoBIC"], color="#4C78A8")
        plt.xticks(x, speede_df["Model"], rotation=30, ha="right")
        plt.ylabel("Pseudo-BIC (lower better)")
        plt.title("SPEEDEMON species delimitation model comparison")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "speedeMON_model_scores.png"), dpi=320)
        plt.close()

    # Random-forest species delimitation model ranking
    log("Running random-forest species delimitation model ranking (IQ-TREE guided)")
    rf_rows = []
    for mn in model_names:
        labs = models[mn]
        classes = {k: i for i, k in enumerate(sorted(set(labs)))}
        y = np.array([classes[z] for z in labs], dtype=int)
        rf = SimpleRF(n_trees=90, max_depth=8, min_leaf=3, mtry=min(28, int(math.sqrt(Xs.shape[1]))), seed=20260226)
        rf.fit(Xs, y)
        oob, counts = rf.oob_proba(Xs)
        pred = np.argmax(oob, axis=1)
        acc = np.mean(pred == y)
        bacc = balanced_accuracy(y, pred)
        ll = -np.mean(np.log(np.clip(oob[np.arange(y.size), y], 1e-9, 1.0)))
        k = len(classes)
        rf_penalized = ll + 0.03 * (k - 1)
        guide_penalty = guide_penalty_weight * guide_nonmono[mn]
        final_score = rf_penalized + guide_penalty
        rf_rows.append(
            (
                mn, k, acc, bacc, ll, rf_penalized, guide_nonmono[mn], guide_mono[mn],
                guide_penalty, final_score, np.mean(counts > 0)
            )
        )

    rf_df = pd.DataFrame(rf_rows, columns=[
        "Model", "K_groups", "OOB_accuracy", "OOB_balanced_accuracy",
        "OOB_logloss", "RF_penalized_score_raw", "Guide_nonmonophyly",
        "Guide_monophyly_score", "Guide_penalty", "Penalized_score", "Frac_with_OOB"
    ])
    rf_df = rf_df.sort_values("Penalized_score").reset_index(drop=True)
    rf_df.to_csv(os.path.join(TAB_DIR, "random_forest_species_delimitation_scores.tsv"), sep="\t", index=False)
    best_rf = rf_df.iloc[0]["Model"]
    pd.DataFrame({"Sample": samples, "Population": populations, "Delimited_species": models[best_rf]}).to_csv(
        os.path.join(TAB_DIR, "random_forest_best_model_assignments.tsv"), sep="\t", index=False
    )

    if plt is not None:
        plt.figure(figsize=(7.2, 4.8))
        x = np.arange(len(rf_df))
        plt.bar(x, rf_df["Penalized_score"], color="#F58518")
        plt.xticks(x, rf_df["Model"], rotation=30, ha="right")
        plt.ylabel("RF + IQ-TREE guide score (lower better)")
        plt.title("Random forest species delimitation model comparison (IQ-TREE guided)")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "random_forest_model_scores.png"), dpi=320)
        plt.close()

    # DelimitR-style model choice with summary stats + RF
    log("Running DelimitR-style model-choice analysis")
    feat_obs = {}
    train_X = []
    train_y = []
    sims_per_model = 220
    for mi, mn in enumerate(model_names):
        freq, nind, gr = grouped_freqs(X_imp, models[mn])
        obs_vec = summarize_from_freq(freq, nind)
        feat_obs[mn] = obs_vec
        for _ in range(sims_per_model):
            fs = simulate_freq(freq, nind, conc=30.0)
            vec = summarize_from_freq(fs, nind)
            train_X.append(vec)
            train_y.append(mi)
    train_X = np.array(train_X, dtype=float)
    train_y = np.array(train_y, dtype=int)

    rf_m = SimpleRF(n_trees=120, max_depth=9, min_leaf=5, mtry=min(5, train_X.shape[1]), seed=20260226)
    rf_m.fit(train_X, train_y)
    oob, _ = rf_m.oob_proba(train_X)
    oob_pred = np.argmax(oob, axis=1)
    model_choice_oob_acc = np.mean(oob_pred == train_y)

    rows = []
    for mn in model_names:
        pr = rf_m.predict_proba(feat_obs[mn].reshape(1, -1))[0]
        for mi, mname in enumerate(model_names):
            rows.append((mn, mname, pr[mi]))
    post_df = pd.DataFrame(rows, columns=["Observed_feature_from_model", "Predicted_model", "Posterior_like"])
    diag = post_df[post_df["Observed_feature_from_model"] == post_df["Predicted_model"]].copy()
    diag = diag.sort_values("Posterior_like", ascending=False).reset_index(drop=True)
    diag["DelimitR_rank"] = np.arange(1, diag.shape[0] + 1)
    diag["RF_model_choice_OOB_accuracy"] = model_choice_oob_acc
    diag.to_csv(os.path.join(TAB_DIR, "delimitR_model_support.tsv"), sep="\t", index=False)
    post_df.to_csv(os.path.join(TAB_DIR, "delimitR_full_posterior_matrix.tsv"), sep="\t", index=False)

    if plt is not None:
        plt.figure(figsize=(7.2, 4.8))
        plt.bar(np.arange(diag.shape[0]), diag["Posterior_like"], color="#54A24B")
        plt.xticks(np.arange(diag.shape[0]), diag["Predicted_model"], rotation=30, ha="right")
        plt.ylim(0, 1.0)
        plt.ylabel("Posterior-like support")
        plt.title("DelimitR model support (RF on summary simulations)")
        plt.tight_layout()
        plt.savefig(os.path.join(FIG_DIR, "delimitR_model_support.png"), dpi=320)
        plt.close()

    with open(os.path.join(TAB_DIR, "analysis_summary.txt"), "w") as fh:
        fh.write("Species delimitation analyses (unlinked SNPs, 1 SNP per locus)\n")
        fh.write(f"Date: {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n")
        fh.write(f"Samples: {X_imp.shape[1]}\n")
        fh.write(f"Unlinked SNPs: {X_imp.shape[0]}\n")
        fh.write(f"SPEEDEMON best model: {best_s}\n")
        fh.write(f"Random forest best model: {best_rf}\n")
        fh.write("Random forest scoring used IQ-TREE guide-tree monophyly penalty\n")
        fh.write(f"DelimitR top model: {diag.iloc[0]['Predicted_model']}\n")
        fh.write(f"DelimitR OOB model-choice accuracy: {model_choice_oob_acc:.4f}\n")

    log("Completed all requested species delimitation analyses")


if __name__ == "__main__":
    main()
