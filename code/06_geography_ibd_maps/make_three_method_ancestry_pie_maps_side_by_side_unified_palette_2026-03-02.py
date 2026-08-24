#!/tmp/eems_py/bin/python
import os
import csv
import math
import numpy as np
from PIL import Image, ImageDraw, ImageFont

BASE_EXCL = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
DOWN_EXCL = os.path.join(BASE_EXCL, "downstream_analysis_2026-02-23")
FIG_DIR = os.path.join(DOWN_EXCL, "figures")
TAB_DIR = os.path.join(DOWN_EXCL, "tables")

BASE_COMPLETE = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_2026-02-23_complete"
MAP_BASE_IMG = os.path.join(
    BASE_COMPLETE,
    "downstream_analysis_2026-02-23",
    "additional_analyses_2026-02-24",
    "figures",
    "Map_sample_locations_by_population_publication.png",
)

COORD_FILE = "/Users/gspellman/Trowbridgii_analyses/Sample_geographic_coordinants.txt"
POPMAP_FILE = os.path.join(BASE_EXCL, "popmap.tsv")

ADMIX_Q = os.path.join(DOWN_EXCL, "tables", "admixture_nmf_Qmatrix_bestK.tsv")
FAST_BEST = os.path.join(DOWN_EXCL, "faststructure_analysis_2026-02-27", "tables", "faststructure_10rep_best_k_selection.tsv")
FAST_FAM = os.path.join(DOWN_EXCL, "faststructure_analysis_2026-02-27", "work", "runs", "unlinked.fam")
FAST_RUN_DIR = os.path.join(DOWN_EXCL, "faststructure_analysis_2026-02-27", "work", "runs")
POPCL_Q = os.path.join(DOWN_EXCL, "popcluster_native_2026-02-27", "tables", "popcluster_native_Qmatrix_bestK.tsv")

OUT_PNG = os.path.join(FIG_DIR, "Map_Admixture_fastStructure_PopCluster_pies_side_by_side_unified_palette_publication.png")
OUT_PDF = os.path.join(FIG_DIR, "Map_Admixture_fastStructure_PopCluster_pies_side_by_side_unified_palette_publication.pdf")
OUT_SUMMARY = os.path.join(TAB_DIR, "map_three_method_ancestry_pies_unified_palette_summary.txt")

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)


def load_font(size=34, bold=False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for p in candidates:
        if p and os.path.exists(p):
            try:
                return ImageFont.truetype(p, size=size)
            except Exception:
                pass
    return ImageFont.load_default()


def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))


def detect_map_bbox(img_rgb):
    arr = np.array(img_rgb)
    dark = (arr[:, :, 0] < 70) & (arr[:, :, 1] < 70) & (arr[:, :, 2] < 70)
    row_counts = dark.sum(axis=1)
    col_counts = dark.sum(axis=0)

    rows = np.where(row_counts > row_counts.max() * 0.6)[0]
    cols = np.where(col_counts > col_counts.max() * 0.6)[0]
    if len(rows) < 2 or len(cols) < 2:
        return (906, 114, 1845, 2195)

    return (int(cols.min()), int(rows.min()), int(cols.max()), int(rows.max()))


def map_extent(meta_rows, pad_x=0.8, pad_y=0.6):
    lons = [r["Longitude"] for r in meta_rows]
    lats = [r["Latitude"] for r in meta_rows]
    lon_min = float(min(lons) - pad_x)
    lon_max = float(max(lons) + pad_x)
    lat_min = float(min(lats) - pad_y)
    lat_max = float(max(lats) + pad_y)
    return lon_min, lon_max, lat_min, lat_max


def lonlat_to_px(lon, lat, bbox, extent):
    x0, y0, x1, y1 = bbox
    lon_min, lon_max, lat_min, lat_max = extent
    x = x0 + (lon - lon_min) * (x1 - x0) / (lon_max - lon_min)
    y = y1 - (lat - lat_min) * (y1 - y0) / (lat_max - lat_min)
    return float(x), float(y)


def unified_cluster_palette_8():
    # Single palette used across all methods (up to 8 clusters).
    return [
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
        "#9467bd", "#8c564b", "#e377c2", "#17becf",
    ]


def draw_pie(draw, cx, cy, r, fracs, colors):
    start = 0.0
    box = [cx - r, cy - r, cx + r, cy + r]
    for frac, col in zip(fracs, colors):
        if frac <= 0:
            continue
        end = start + float(frac) * 360.0
        draw.pieslice(box, start=start, end=end, fill=col, outline=(255, 255, 255), width=2)
        start = end
    draw.ellipse(box, outline=(255, 255, 255), width=2)


def pie_radius_px(meta_rows, bbox, extent):
    pts = [lonlat_to_px(r["Longitude"], r["Latitude"], bbox, extent) for r in meta_rows]
    pts = np.array(pts)
    if len(pts) <= 1:
        return 18
    dd = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2))
    dd[dd == 0] = np.nan
    nn = np.nanmin(dd, axis=1)
    med = float(np.nanmedian(nn))
    # Original scale ~0.38*median; this doubles that size.
    r = int(round(med * 0.76))
    return max(16, min(34, r))


def read_popmap(path):
    out = []
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            t = s.split("\t")
            if len(t) >= 2:
                out.append((t[0], t[1]))
    return out


def read_coords(path):
    out = {}
    with open(path, "r", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            try:
                out[r["Sample"]] = (float(r["Latitude"]), float(r["Longitude"]))
            except Exception:
                pass
    return out


def read_q_tsv(path):
    with open(path, "r", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        rows = [dict(r) for r in rd]
    cluster_cols = [c for c in (rows[0].keys() if rows else []) if c.startswith("Cluster")]
    q_by_sample = {}
    for r in rows:
        sid = r.get("Sample")
        if not sid:
            continue
        vals = []
        for c in cluster_cols:
            try:
                vals.append(float(r[c]))
            except Exception:
                vals.append(0.0)
        q_by_sample[sid] = vals
    return cluster_cols, q_by_sample


def read_faststructure_best(path):
    with open(path, "r", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        rows = [r for r in rd]
    k = None
    rep = None
    for r in rows:
        if r.get("criterion") == "selected_bestK":
            k = int(float(r.get("value", "0")))
        if r.get("criterion") == "selected_best_rep":
            rep = int(float(r.get("value", "0")))
    if k is None or rep is None:
        raise RuntimeError("Could not parse fastStructure best K/rep")
    return k, rep


def read_faststructure_q(best_k, best_rep):
    fam_samples = []
    with open(FAST_FAM, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            fam_samples.append(s.split()[1])

    meanq_path = os.path.join(FAST_RUN_DIR, f"fs10_rep{best_rep}.{best_k}.meanQ")
    if not os.path.exists(meanq_path):
        raise FileNotFoundError(meanq_path)

    mat = []
    with open(meanq_path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            mat.append([float(x) for x in s.split()])

    if len(mat) != len(fam_samples):
        raise RuntimeError("fastStructure meanQ rows do not match FAM sample count")

    cluster_cols = [f"Cluster{i+1}" for i in range(len(mat[0]))]
    q_by_sample = {sid: vals for sid, vals in zip(fam_samples, mat)}
    return cluster_cols, q_by_sample, meanq_path


def clean_base_map(base_img):
    panel = base_img.copy().convert("RGB")
    d = ImageDraw.Draw(panel)
    # Remove old embedded legends/title artifacts.
    legend_box = (888, 92, 1338, 478)
    donor_box = (1346, 92, 1796, 478)
    panel.paste(panel.crop(donor_box), (legend_box[0], legend_box[1]))
    d.rectangle([560, 18, 2200, 108], fill=(255, 255, 255))
    return panel


def make_panel(base_img, meta_rows, q_by_sample, n_clusters, title, bbox, extent, palette8):
    panel = clean_base_map(base_img)
    draw = ImageDraw.Draw(panel)
    draw.text((620, 32), title, fill=(20, 20, 20), font=load_font(38, bold=False))

    colors = [hex_to_rgb(palette8[i]) for i in range(n_clusters)]
    r = pie_radius_px(meta_rows, bbox, extent)

    for row in meta_rows:
        sid = row["Sample"]
        vals = q_by_sample.get(sid)
        if vals is None:
            continue
        vals = np.array(vals[:n_clusters], dtype=float)
        s = float(vals.sum())
        if s <= 0:
            continue
        vals = vals / s
        x, y = lonlat_to_px(row["Longitude"], row["Latitude"], bbox, extent)
        draw_pie(draw, x, y, r, vals, colors)

    return panel, r


def draw_single_legend(draw, x, y, title, palette8):
    tfont = load_font(30, bold=True)
    lfont = load_font(25, bold=False)
    draw.text((x, y), title, fill=(20, 20, 20), font=tfont)
    yy = y + 48
    for i, h in enumerate(palette8, start=1):
        col = hex_to_rgb(h)
        draw.rectangle([x, yy + 2, x + 34, yy + 26], fill=col, outline=(40, 40, 40), width=1)
        draw.text((x + 46, yy), f"Cluster{i}", fill=(20, 20, 20), font=lfont)
        yy += 32


def main():
    req = [MAP_BASE_IMG, COORD_FILE, POPMAP_FILE, ADMIX_Q, FAST_BEST, FAST_FAM, POPCL_Q]
    miss = [p for p in req if not os.path.exists(p)]
    if miss:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(miss))

    base = Image.open(MAP_BASE_IMG).convert("RGB")
    bbox = detect_map_bbox(base)

    popmap = read_popmap(POPMAP_FILE)
    coords = read_coords(COORD_FILE)

    meta_rows = []
    for sid, pop in popmap:
        if sid in coords:
            lat, lon = coords[sid]
            meta_rows.append({"Sample": sid, "Population": pop, "Latitude": lat, "Longitude": lon})
    meta_rows = sorted(meta_rows, key=lambda r: (r["Population"], r["Sample"]))

    extent = map_extent(meta_rows)

    # Load method Q matrices.
    admix_cols, admix_q = read_q_tsv(ADMIX_Q)
    k_admix = len(admix_cols)

    k_fast, rep_fast = read_faststructure_best(FAST_BEST)
    fast_cols, fast_q, meanq_path = read_faststructure_q(k_fast, rep_fast)

    popcl_cols, popcl_q = read_q_tsv(POPCL_Q)
    k_popcl = len(popcl_cols)

    # Shared palette up to 8 clusters.
    palette8 = unified_cluster_palette_8()

    p_admix, r_admix = make_panel(base, meta_rows, admix_q, min(8, k_admix), f"Admixture (K={k_admix})", bbox, extent, palette8)
    p_fast, r_fast = make_panel(base, meta_rows, fast_q, min(8, k_fast), f"fastStructure (K={k_fast})", bbox, extent, palette8)
    p_popc, r_popc = make_panel(base, meta_rows, popcl_q, min(8, k_popcl), f"PopCluster (K={k_popcl})", bbox, extent, palette8)

    target_h = 2050
    scale = target_h / float(p_admix.height)
    target_w = int(round(p_admix.width * scale))
    p_admix = p_admix.resize((target_w, target_h), Image.Resampling.LANCZOS)
    p_fast = p_fast.resize((target_w, target_h), Image.Resampling.LANCZOS)
    p_popc = p_popc.resize((target_w, target_h), Image.Resampling.LANCZOS)

    leg_w = 330
    side_pad = 18
    top_pad = 80
    bot_pad = 28
    gap = 18

    out_w = side_pad * 2 + leg_w + target_w * 3 + gap * 2
    out_h = top_pad + target_h + bot_pad
    canvas = Image.new("RGB", (out_w, out_h), (255, 255, 255))

    x_leg = side_pad
    x1 = x_leg + leg_w
    x2 = x1 + target_w + gap
    x3 = x2 + target_w + gap
    y0 = top_pad

    canvas.paste(p_admix, (x1, y0))
    canvas.paste(p_fast, (x2, y0))
    canvas.paste(p_popc, (x3, y0))

    d = ImageDraw.Draw(canvas)
    d.text((int(out_w * 0.23), 12), "Geographic ancestry proportions: Admixture, fastStructure, PopCluster", fill=(20, 20, 20), font=load_font(44, bold=True))

    legend_h = 48 + 8 * 32 + 8
    y_leg = y0 + int((target_h - legend_h) / 2)
    draw_single_legend(d, x_leg + 12, y_leg, "Clusters (shared)", palette8)

    canvas.save(OUT_PNG, dpi=(300, 300))
    canvas.save(OUT_PDF, resolution=300)

    with open(OUT_SUMMARY, "w") as f:
        f.write("Three-method side-by-side ancestry pie maps with unified 8-color palette\n")
        f.write(f"Base map image: {MAP_BASE_IMG}\n")
        f.write(f"Detected map bbox: {bbox}\n")
        f.write(f"Geographic extent used: {extent}\n")
        f.write(f"Admixture K={k_admix}\n")
        f.write(f"fastStructure K={k_fast}, rep={rep_fast}, meanQ={meanq_path}\n")
        f.write(f"PopCluster K={k_popcl}\n")
        f.write(f"Pie radius (Admixture): {r_admix} px\n")
        f.write(f"Pie radius (fastStructure): {r_fast} px\n")
        f.write(f"Pie radius (PopCluster): {r_popc} px\n")
        f.write("Shared palette:\n")
        for i, h in enumerate(palette8, start=1):
            f.write(f"  Cluster{i}\t{h}\n")
        f.write(f"Output PNG: {OUT_PNG}\n")
        f.write(f"Output PDF: {OUT_PDF}\n")


if __name__ == "__main__":
    main()
