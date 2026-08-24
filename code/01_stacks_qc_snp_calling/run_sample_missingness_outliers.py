#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
DOWN = os.path.join(BASE, "downstream_analysis_2026-02-23")
VCF = os.path.join(BASE, "populations.snps.vcf")
POPMAP = os.path.join(BASE, "popmap.tsv")
TAB = os.path.join(DOWN, "tables")
FIG = os.path.join(DOWN, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)


def is_missing(gt_field: str) -> bool:
    gt = gt_field.split(":", 1)[0]
    gt = gt.replace("|", "/")
    return gt in {"./.", ".", "./."}


# Parse VCF and count missing genotypes per sample
samples = None
missing_counts = None
total_sites = 0

with open(VCF, "r") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            samples = header[9:]
            missing_counts = np.zeros(len(samples), dtype=np.int64)
            continue
        if line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 10:
            continue
        gts = parts[9:]
        total_sites += 1
        for i, g in enumerate(gts):
            if is_missing(g):
                missing_counts[i] += 1

if samples is None or total_sites == 0:
    raise RuntimeError("Failed to read samples/sites from VCF")

missing_prop = missing_counts / float(total_sites)
callrate = 1.0 - missing_prop

# Join population labels
pop = pd.read_csv(POPMAP, sep="\t")
pop.columns = ["Sample", "Population"]
df = pd.DataFrame({
    "Sample": samples,
    "N_sites": total_sites,
    "N_missing": missing_counts,
    "Missingness": missing_prop,
    "Callrate": callrate,
})
df = df.merge(pop, on="Sample", how="left")

# Outlier detection (IQR rule)
q1 = float(df["Missingness"].quantile(0.25))
q3 = float(df["Missingness"].quantile(0.75))
iqr = q3 - q1
upper = q3 + 1.5 * iqr
lower = max(0.0, q1 - 1.5 * iqr)

df["Outlier_IQR"] = (df["Missingness"] > upper) | (df["Missingness"] < lower)

# z-score as secondary diagnostic
mu = float(df["Missingness"].mean())
sd = float(df["Missingness"].std(ddof=1))
if sd > 0:
    z = (df["Missingness"] - mu) / sd
else:
    z = np.zeros(len(df))
df["Missingness_Z"] = z
df["Outlier_Z3"] = np.abs(df["Missingness_Z"]) > 3.0

# Write tables
df = df.sort_values("Missingness", ascending=False).reset_index(drop=True)
full_table = os.path.join(TAB, "sample_missingness_full_dataset.tsv")
outlier_table = os.path.join(TAB, "sample_missingness_outliers.tsv")
summary_file = os.path.join(TAB, "sample_missingness_summary.txt")

df.to_csv(full_table, sep="\t", index=False)
df[df["Outlier_IQR"] | df["Outlier_Z3"]].to_csv(outlier_table, sep="\t", index=False)

with open(summary_file, "w") as fh:
    fh.write("Per-sample missingness summary (full dataset)\n")
    fh.write(f"VCF: {VCF}\n")
    fh.write(f"Samples: {len(df)}\n")
    fh.write(f"Sites: {total_sites}\n")
    fh.write(f"Mean missingness: {df['Missingness'].mean():.6f}\n")
    fh.write(f"Median missingness: {df['Missingness'].median():.6f}\n")
    fh.write(f"Q1: {q1:.6f}\n")
    fh.write(f"Q3: {q3:.6f}\n")
    fh.write(f"IQR: {iqr:.6f}\n")
    fh.write(f"IQR lower bound: {lower:.6f}\n")
    fh.write(f"IQR upper bound: {upper:.6f}\n")
    fh.write(f"IQR outliers: {int(df['Outlier_IQR'].sum())}\n")
    fh.write(f"|Z|>3 outliers: {int(df['Outlier_Z3'].sum())}\n")

# Plot histogram
plt.rcParams.update({"font.size": 11})
fig, ax = plt.subplots(figsize=(9, 6))
vals_pct = df["Missingness"].to_numpy() * 100.0
ax.hist(vals_pct, bins=20, color="#4C78A8", edgecolor="white", alpha=0.95)
ax.axvline(upper * 100.0, color="#D62728", linestyle="--", linewidth=1.8, label=f"IQR upper ({upper*100:.2f}%)")
ax.axvline(df["Missingness"].median() * 100.0, color="#2F2F2F", linestyle="-", linewidth=1.2, label=f"Median ({df['Missingness'].median()*100:.2f}%)")

out = df[df["Outlier_IQR"]]
if len(out):
    for _, r in out.iterrows():
        ax.axvline(r["Missingness"] * 100.0, color="#E45756", alpha=0.5, linewidth=1.0)

ax.set_xlabel("Per-sample missingness (%)")
ax.set_ylabel("Number of samples")
ax.set_title("Per-sample missingness distribution (full SNP dataset)")
ax.legend(frameon=True, fontsize=9)
fig.tight_layout()

png = os.path.join(FIG, "sample_missingness_histogram_full_dataset.png")
pdf = os.path.join(FIG, "sample_missingness_histogram_full_dataset.pdf")
fig.savefig(png, dpi=400)
fig.savefig(pdf)
plt.close(fig)

print("Wrote:")
print(full_table)
print(outlier_table)
print(summary_file)
print(png)
print(pdf)
