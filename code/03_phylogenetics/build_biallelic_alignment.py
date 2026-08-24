#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd

BASE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
OUT = os.path.join(BASE, "downstream_analysis_2026-02-23", "additional_analyses_2026-02-24", "phylogeny_ml_biallelic_2026-02-24")
WORK = os.path.join(OUT, "work")
TAB = os.path.join(OUT, "tables")
os.makedirs(WORK, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

VCF = os.path.join(BASE, "populations.snps.vcf")
POPMAP = os.path.join(BASE, "popmap.tsv")

IUPAC = {
    frozenset(("A", "G")): "R",
    frozenset(("C", "T")): "Y",
    frozenset(("G", "C")): "S",
    frozenset(("A", "T")): "W",
    frozenset(("G", "T")): "K",
    frozenset(("A", "C")): "M",
}


def gt_to_base(gt, ref, alt):
    a = gt.split(":", 1)[0].replace("|", "/")
    if a in ("./.", "."):
        return "N", np.nan
    if a == "0/0":
        return ref, 0.0
    if a == "1/1":
        return alt, 2.0
    if a in ("0/1", "1/0"):
        return IUPAC.get(frozenset((ref, alt)), "N"), 1.0
    return "N", np.nan


pop = pd.read_csv(POPMAP, sep="\t", header=None, names=["Sample", "Population"])
sample_to_pop = dict(zip(pop["Sample"], pop["Population"]))

samples = None
seqs = None
kept = 0
all_sites = 0

with open(VCF, "r") as fh:
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            hdr = line.rstrip("\n").split("\t")
            samples = hdr[9:]
            seqs = {s: [] for s in samples}
            continue

        p = line.rstrip("\n").split("\t")
        ref = p[3].upper()
        alt = p[4].upper()
        if len(ref) != 1 or len(alt) != 1 or "," in alt:
            continue
        if ref not in "ACGT" or alt not in "ACGT":
            continue

        all_sites += 1
        gts = p[9:]

        dos = []
        bases = []
        for gt in gts:
            b, d = gt_to_base(gt, ref, alt)
            bases.append(b)
            dos.append(d)

        dos = np.array(dos, dtype=float)
        miss = np.mean(~np.isfinite(dos))
        if miss > 0.20:
            continue

        p_alt = np.nanmean(dos) / 2.0
        maf = min(p_alt, 1.0 - p_alt)
        if not np.isfinite(maf) or maf < 0.05:
            continue

        kept += 1
        for s, b in zip(samples, bases):
            seqs[s].append(b)

aln_len = kept
if aln_len == 0:
    raise SystemExit("No sites passed filters")

phy = os.path.join(WORK, "biallelic_individuals_relaxed.phy")
fa = os.path.join(WORK, "biallelic_individuals.fasta")
meta = os.path.join(TAB, "ml_biallelic_sample_metadata.tsv")
summary = os.path.join(TAB, "ml_biallelic_alignment_summary.tsv")

with open(phy, "w") as out:
    out.write(f"{len(samples)} {aln_len}\n")
    for s in samples:
        out.write(f"{s} {''.join(seqs[s])}\n")

with open(fa, "w") as out:
    for s in samples:
        out.write(f">{s}\n{''.join(seqs[s])}\n")

meta_df = pd.DataFrame({
    "Sample": samples,
    "Population": [sample_to_pop.get(s, "Unknown") for s in samples],
})
meta_df.to_csv(meta, sep="\t", index=False)

pd.DataFrame([{
    "N_samples": len(samples),
    "Sites_total_biallelic_snp": all_sites,
    "Sites_retained_after_missing_maf_filter": kept,
    "Missingness_filter_max": 0.20,
    "MAF_filter_min": 0.05,
}]).to_csv(summary, sep="\t", index=False)

print(f"Wrote {phy}")
print(f"Retained sites: {kept} / {all_sites}")
