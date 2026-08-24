#!/usr/bin/env python3
import os
import csv
import math
import struct
from statistics import median
from PIL import Image, ImageDraw, ImageFont, ImageFilter, ImageEnhance

BASE_EXCL = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_excl_MVZ216210_2026-02-26_clean"
DOWN = os.path.join(BASE_EXCL, "downstream_analysis_2026-02-23")
FIG_DIR = os.path.join(DOWN, "figures")
TAB_DIR = os.path.join(DOWN, "tables")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

COORD_FILE = "/Users/gspellman/Trowbridgii_analyses/Sample_geographic_coordinants.txt"
POPMAP_FILE = os.path.join(BASE_EXCL, "popmap.tsv")
ADMIX_Q = os.path.join(DOWN, "tables", "admixture_nmf_Qmatrix_bestK.tsv")
FAST_BEST = os.path.join(DOWN, "faststructure_analysis_2026-02-27", "tables", "faststructure_10rep_best_k_selection.tsv")
FAST_FAM = os.path.join(DOWN, "faststructure_analysis_2026-02-27", "work", "runs", "unlinked.fam")
FAST_RUNS = os.path.join(DOWN, "faststructure_analysis_2026-02-27", "work", "runs")
POPCL_Q = os.path.join(DOWN, "popcluster_native_2026-02-27", "tables", "popcluster_native_Qmatrix_bestK.tsv")
POPCL_BEST = os.path.join(DOWN, "popcluster_native_2026-02-27", "tables", "popcluster_native_best_k.tsv")
SPEEDEMON_ASSIGN = os.path.join(
    DOWN,
    "species_delimitation_speede_rf_delimitr_2026-02-26",
    "tables",
    "speedeMON_best_model_assignments.tsv",
)

COAST_SHP = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/physical/ne_10m_coastline.shp")
STATE_SHP = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/cultural/ne_10m_admin_1_states_provinces_lakes.shp")
STATE_DBF = os.path.expanduser("~/.local/share/cartopy/shapefiles/natural_earth/cultural/ne_10m_admin_1_states_provinces_lakes.dbf")
TOPO_TIF = "/Users/gspellman/Trowbridgii_analyses/Stacks analysis of ddRAD data/stacks_refmap_sorex_2026-02-23_complete/downstream_analysis_2026-02-23/additional_analyses_2026-02-24/work/topography/HYP_50M_SR_W.tif"

OUT_PNG = os.path.join(FIG_DIR, "Map_Admixture_fastStructure_PopCluster_pies_side_by_side_unified_palette_publication.png")
OUT_PDF = os.path.join(FIG_DIR, "Map_Admixture_fastStructure_PopCluster_pies_side_by_side_unified_palette_publication.pdf")
OUT_TXT = os.path.join(TAB_DIR, "map_three_method_ancestry_pies_unified_palette_summary.txt")


def load_font(size=28, bold=False):
    cands = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
    ]
    for c in cands:
        if c and os.path.exists(c):
            try:
                return ImageFont.truetype(c, size=size)
            except Exception:
                pass
    return ImageFont.load_default()


def hex_to_rgb(h):
    h = h.lstrip("#")
    return (int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))


def unified_palette_8():
    return [
        "#1f78b4", "#e31a1c", "#33a02c", "#ff7f00",
        "#6a3d9a", "#17becf", "#a65628", "#fb9a99",
    ]


def read_tsv_dicts(path):
    with open(path, "r", newline="") as f:
        rd = csv.DictReader(f, delimiter="\t")
        return [dict(r) for r in rd]


def read_popmap(path):
    out = []
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            parts = s.split("\t")
            if len(parts) >= 2:
                out.append({"Sample": parts[0], "Population": parts[1]})
    return out


def read_ws(path):
    out = []
    with open(path, "r") as f:
        for ln in f:
            s = ln.strip()
            if s:
                out.append(s.split())
    return out


def read_fast_best(path):
    rows = read_tsv_dicts(path)
    bestk, bestrep = None, None
    for r in rows:
        if r.get("criterion") == "selected_bestK":
            bestk = int(float(r.get("value", "0")))
        elif r.get("criterion") == "selected_best_rep":
            bestrep = int(float(r.get("value", "0")))
    if bestk is None or bestrep is None:
        raise RuntimeError("Could not parse fastStructure bestK/bestRep")
    return bestk, bestrep


def read_popcluster_best_k(path):
    rows = read_tsv_dicts(path)
    if not rows:
        return None
    r0 = rows[0]
    for key in ("best_k", "best_k_DLK2", "best_k_FSTFIS"):
        if key in r0 and str(r0[key]).strip():
            try:
                return int(float(r0[key]))
            except Exception:
                continue
    return None


def read_speede_assignments(path):
    out = {}
    rows = read_tsv_dicts(path)
    for r in rows:
        sid = str(r.get("Sample", "")).strip()
        grp = str(r.get("Delimited_species", "")).strip()
        if sid and grp:
            out[sid] = grp
    return out


def read_shp_parts(path):
    parts_all = []
    with open(path, "rb") as f:
        f.read(100)
        while True:
            rh = f.read(8)
            if len(rh) < 8:
                break
            _rnum, rlen_words = struct.unpack(">2i", rh)
            rec = f.read(rlen_words * 2)
            if len(rec) < 4:
                break
            stype = struct.unpack("<i", rec[0:4])[0]
            if stype not in (3, 5):
                parts_all.append([])
                continue
            off = 4 + 32
            nparts, npts = struct.unpack("<2i", rec[off:off + 8])
            off += 8
            idx = struct.unpack("<" + "i" * nparts, rec[off:off + 4 * nparts]) if nparts else ()
            off += 4 * nparts
            pts = struct.unpack("<" + "d" * (2 * npts), rec[off:off + 16 * npts]) if npts else ()
            xy = [(pts[i], pts[i + 1]) for i in range(0, len(pts), 2)]
            rec_parts = []
            for p in range(nparts):
                s = idx[p]
                e = idx[p + 1] if p + 1 < nparts else npts
                rec_parts.append(xy[s:e])
            parts_all.append(rec_parts)
    return parts_all


def read_dbf(path):
    with open(path, "rb") as f:
        head = f.read(32)
        nrec = struct.unpack("<I", head[4:8])[0]
        rlen = struct.unpack("<H", head[10:12])[0]
        fields = []
        while True:
            b = f.read(1)
            if b == b"\r":
                break
            chunk = b + f.read(31)
            name = chunk[0:11].split(b"\x00", 1)[0].decode("latin1")
            ftype = chr(chunk[11])
            flen = chunk[16]
            fields.append((name, ftype, flen))
        rows = []
        for _ in range(nrec):
            rec = f.read(rlen)
            if not rec or rec[0:1] == b"*":
                continue
            pos = 1
            row = {}
            for name, _ftype, flen in fields:
                raw = rec[pos:pos + flen]
                pos += flen
                row[name] = raw.decode("latin1", errors="ignore").strip()
            rows.append(row)
    return rows


def lonlat_to_xy(lon, lat, extent, frame):
    lon_min, lon_max, lat_min, lat_max = extent
    x0, y0, x1, y1 = frame
    x = x0 + (lon - lon_min) * (x1 - x0) / (lon_max - lon_min)
    y = y1 - (lat - lat_min) * (y1 - y0) / (lat_max - lat_min)
    return x, y


def draw_poly(draw, pts, extent, frame, color, width=1):
    out = [lonlat_to_xy(lon, lat, extent, frame) for lon, lat in pts]
    if len(out) >= 2:
        draw.line(out, fill=color, width=width)


def nice_ticks(vmin, vmax, step):
    s = math.floor(vmin / step) * step
    out = []
    x = s
    while x <= vmax + 1e-8:
        if x >= vmin - 1e-8:
            out.append(round(x, 6))
        x += step
    return out


def fmt_lon(v):
    return f"{abs(v):.1f}°W" if v < 0 else f"{v:.1f}°E"


def fmt_lat(v):
    return f"{abs(v):.1f}°S" if v < 0 else f"{v:.1f}°N"


def draw_pie(draw, cx, cy, r, fracs, cols):
    box = [cx - r, cy - r, cx + r, cy + r]
    start = 0.0
    for frac, col in zip(fracs, cols):
        if frac <= 0:
            continue
        end = start + float(frac) * 360.0
        draw.pieslice(box, start, end, fill=col, outline=(255, 255, 255), width=2)
        start = end
    draw.ellipse(box, outline=(25, 25, 25), width=2)


def terrain_crop(extent, out_size):
    lon_min, lon_max, lat_min, lat_max = extent
    topo = Image.open(TOPO_TIF).convert("RGB")
    tw, th = topo.size
    x0 = int((lon_min + 180.0) / 360.0 * tw)
    x1 = int((lon_max + 180.0) / 360.0 * tw)
    y0 = int((90.0 - lat_max) / 180.0 * th)
    y1 = int((90.0 - lat_min) / 180.0 * th)
    x0, x1 = max(0, min(tw - 1, x0)), max(1, min(tw, x1))
    y0, y1 = max(0, min(th - 1, y0)), max(1, min(th, y1))
    crop = topo.crop((x0, y0, x1, y1)).resize(out_size, Image.Resampling.LANCZOS)

    gray = crop.convert("L")
    gray = ImageEnhance.Contrast(gray).enhance(1.65)
    gray = ImageEnhance.Sharpness(gray).enhance(2.8)
    gray = gray.filter(ImageFilter.UnsharpMask(radius=3.4, percent=320, threshold=1))
    gray = gray.filter(ImageFilter.DETAIL)
    gray = gray.filter(ImageFilter.DETAIL)
    gray = ImageEnhance.Brightness(gray).enhance(1.07)
    return Image.merge("RGB", (gray, gray, gray))


def extract_cluster_cols(rows):
    if not rows:
        return []
    return [c for c in rows[0].keys() if c.startswith("Cluster")]


def qmap_by_sample(rows, cols):
    out = {}
    for r in rows:
        sid = r.get("Sample")
        if not sid:
            continue
        vals = []
        for c in cols:
            try:
                vals.append(float(r.get(c, 0.0)))
            except Exception:
                vals.append(0.0)
        out[sid] = vals
    return out


def main():
    req = [
        COORD_FILE, POPMAP_FILE, ADMIX_Q, FAST_BEST, FAST_FAM, POPCL_Q, POPCL_BEST, SPEEDEMON_ASSIGN,
        TOPO_TIF, COAST_SHP, STATE_SHP, STATE_DBF
    ]
    miss = [p for p in req if not os.path.exists(p)]
    if miss:
        raise FileNotFoundError("Missing required files:\n" + "\n".join(miss))

    pop_rows = read_popmap(POPMAP_FILE)
    coord_rows = read_tsv_dicts(COORD_FILE)
    coord_map = {}
    for r in coord_rows:
        sid = str(r.get("Sample", "")).strip()
        if not sid:
            continue
        try:
            coord_map[sid] = (float(r["Latitude"]), float(r["Longitude"]))
        except Exception:
            continue

    meta = []
    for r in pop_rows:
        sid = r["Sample"]
        if sid in coord_map:
            lat, lon = coord_map[sid]
            meta.append({"Sample": sid, "Population": r["Population"], "Latitude": lat, "Longitude": lon})

    lons = [r["Longitude"] for r in meta]
    lats = [r["Latitude"] for r in meta]
    extent = (min(lons) - 1.35, max(lons) + 1.35, min(lats) - 0.65, max(lats) + 0.65)

    admix_rows = read_tsv_dicts(ADMIX_Q)
    admix_cols = extract_cluster_cols(admix_rows)

    best_k, best_rep = read_fast_best(FAST_BEST)
    fs_samples = [r[1] for r in read_ws(FAST_FAM) if len(r) > 1]
    fs_qfile = os.path.join(FAST_RUNS, f"fs10_rep{best_rep}.{best_k}.meanQ")
    fs_mat = read_ws(fs_qfile)
    fs_rows = []
    for sid, vals in zip(fs_samples, fs_mat):
        row = {"Sample": sid}
        for i, v in enumerate(vals):
            row[f"Cluster{i + 1}"] = float(v)
        fs_rows.append(row)
    fs_cols = extract_cluster_cols(fs_rows)

    popcl_rows = read_tsv_dicts(POPCL_Q)
    popcl_cols = extract_cluster_cols(popcl_rows)
    popcl_best_k = read_popcluster_best_k(POPCL_BEST)

    speede_map = read_speede_assignments(SPEEDEMON_ASSIGN)
    speede_species = []
    for mr in meta:
        sp = speede_map.get(mr["Sample"])
        if sp and sp not in speede_species:
            speede_species.append(sp)
    speede_cols = [f"Species{i+1}" for i in range(len(speede_species))]
    speede_rows = []
    for mr in meta:
        sid = mr["Sample"]
        sp = speede_map.get(sid)
        row = {"Sample": sid}
        for i, sp_name in enumerate(speede_species):
            row[speede_cols[i]] = 1.0 if sp == sp_name else 0.0
        speede_rows.append(row)

    methods = [
        (f"Admixture (K={len(admix_cols)})", admix_rows, admix_cols),
        (f"fastStructure (K={best_k})", fs_rows, fs_cols),
        (f"PopCluster (K={popcl_best_k if popcl_best_k is not None else len(popcl_cols)})", popcl_rows, popcl_cols),
        (f"SPEEDEMON species assignment (K={len(speede_species)})", speede_rows, speede_cols),
    ]

    kmax = max(len(admix_cols), len(fs_cols), len(popcl_cols), len(speede_cols))
    if kmax < 1:
        raise RuntimeError("No cluster columns found in ancestry matrices")
    palette = unified_palette_8()[:max(8, kmax)]
    palette_rgb = [hex_to_rgb(c) for c in palette]

    coast_parts = read_shp_parts(COAST_SHP)
    state_parts = read_shp_parts(STATE_SHP)
    state_attrs = read_dbf(STATE_DBF)

    # Layout tuned for tight alignment and a left-anchored legend.
    W, H = 6700, 3200
    top_pad, bot_pad = 140, 170
    map_h = H - top_pad - bot_pad
    map_w = 1470
    gap = 12
    legend_w = 430
    left_pad = 36
    x1 = left_pad + legend_w + 22
    x2 = x1 + map_w + gap
    x3 = x2 + map_w + gap
    x4 = x3 + map_w + gap

    frames = [
        (x1, top_pad, x1 + map_w, top_pad + map_h),
        (x2, top_pad, x2 + map_w, top_pad + map_h),
        (x3, top_pad, x3 + map_w, top_pad + map_h),
        (x4, top_pad, x4 + map_w, top_pad + map_h),
    ]

    img = Image.new("RGB", (W, H), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    # Render one high-res grayscale terrain tile and reuse across panels.
    tile = terrain_crop(extent, (map_w, map_h))

    lon_ticks = nice_ticks(extent[0], extent[1], 1.5)
    lat_ticks = nice_ticks(extent[2], extent[3], 2.0)

    for i, (title, rows, cols) in enumerate(methods):
        frame = frames[i]
        fx0, fy0, fx1, fy1 = frame
        img.paste(tile, (fx0, fy0))

        # Grid and ticks on outside edges only to keep panel alignment tight.
        for lon in lon_ticks:
            x, _ = lonlat_to_xy(lon, extent[2], extent, frame)
            draw.line([(x, fy0), (x, fy1)], fill=(175, 175, 175), width=1)
        for lat in lat_ticks:
            _, y = lonlat_to_xy(extent[0], lat, extent, frame)
            draw.line([(fx0, y), (fx1, y)], fill=(175, 175, 175), width=1)

        if i == 0:
            for lat in lat_ticks:
                _, y = lonlat_to_xy(extent[0], lat, extent, frame)
                draw.text((fx0 - 130, y - 22), fmt_lat(lat), fill=(20, 20, 20), font=load_font(36))
            for lon in lon_ticks:
                x, _ = lonlat_to_xy(lon, extent[2], extent, frame)
                draw.text((x - 62, fy1 + 18), fmt_lon(lon), fill=(20, 20, 20), font=load_font(33))

        # Clipped linework for coast/state boundaries.
        b_layer = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        bdraw = ImageDraw.Draw(b_layer)
        for rec_parts in coast_parts:
            for part in rec_parts:
                draw_poly(bdraw, part, extent, frame, color=(45, 45, 45, 255), width=3)
        for rec_parts, at in zip(state_parts, state_attrs):
            adm0 = (at.get("adm0_name") or at.get("admin") or "").strip()
            if adm0 not in ("United States of America", "United States"):
                continue
            for part in rec_parts:
                draw_poly(bdraw, part, extent, frame, color=(90, 90, 90, 255), width=2)

        clip = Image.new("L", (W, H), 0)
        cdraw = ImageDraw.Draw(clip)
        cdraw.rectangle([fx0, fy0, fx1, fy1], fill=255)
        clipped = Image.new("RGBA", (W, H), (0, 0, 0, 0))
        clipped.paste(b_layer, (0, 0), clip)
        tmp = img.convert("RGBA")
        tmp.alpha_composite(clipped)
        img = tmp.convert("RGB")
        draw = ImageDraw.Draw(img)

        draw.rectangle([fx0, fy0, fx1, fy1], outline=(48, 48, 48), width=2)

        qmap = qmap_by_sample(rows, cols)
        pts = []
        for mr in meta:
            x, y = lonlat_to_xy(mr["Longitude"], mr["Latitude"], extent, frame)
            pts.append((x, y))
        if len(pts) > 1:
            nearest = []
            for i1, (xa, ya) in enumerate(pts):
                dmin = None
                for i2, (xb, yb) in enumerate(pts):
                    if i1 == i2:
                        continue
                    d = math.hypot(xa - xb, ya - yb)
                    if dmin is None or d < dmin:
                        dmin = d
                if dmin is not None:
                    nearest.append(dmin)
            med = float(median(nearest)) if nearest else 26.0
            radius = max(22, min(42, int(round(med * 0.72))))
        else:
            radius = 28

        for mr in meta:
            sid = mr["Sample"]
            vals = qmap.get(sid)
            if not vals:
                continue
            vals = [float(v) for v in vals]
            s = sum(vals)
            if s <= 0:
                continue
            fracs = [v / s for v in vals]
            colors = palette_rgb[:len(fracs)]
            x, y = lonlat_to_xy(mr["Longitude"], mr["Latitude"], extent, frame)
            draw_pie(draw, x, y, radius, fracs, colors)

        draw.text((fx0 + 8, fy0 - 62), title, fill=(15, 15, 15), font=load_font(46, bold=True))

        if i == len(methods) - 1:
            for lat in lat_ticks:
                _, y = lonlat_to_xy(extent[0], lat, extent, frame)
                draw.text((fx1 + 20, y - 22), fmt_lat(lat), fill=(20, 20, 20), font=load_font(36))

    # Shared legends tightly left of panel 1.
    lx0 = left_pad + 10
    ly0 = top_pad + int(map_h * 0.12)
    draw.text((lx0, ly0 - 56), "Clusters", fill=(15, 15, 15), font=load_font(44, bold=True))
    y = ly0 + 8
    for idx in range(1, 9):
        c = palette[idx - 1]
        draw.rectangle([lx0, y, lx0 + 44, y + 44], fill=hex_to_rgb(c), outline=(20, 20, 20), width=2)
        draw.text((lx0 + 58, y - 1), f"Cluster{idx}", fill=(20, 20, 20), font=load_font(34))
        y += 56

    y += 24
    draw.text((lx0, y), "SPEEDEMON groups", fill=(15, 15, 15), font=load_font(42, bold=True))
    y += 58
    for i, sp in enumerate(speede_species):
        c = palette[i]
        draw.rectangle([lx0, y, lx0 + 44, y + 44], fill=hex_to_rgb(c), outline=(20, 20, 20), width=2)
        draw.text((lx0 + 58, y - 1), sp, fill=(20, 20, 20), font=load_font(30))
        y += 52

    img.save(OUT_PNG, dpi=(400, 400))
    img.save(OUT_PDF, resolution=400)

    with open(OUT_TXT, "w") as f:
        f.write("Four-method ancestry/species-assignment maps (high-resolution grayscale topography renderer).\n")
        f.write("Panels: Admixture, fastStructure, PopCluster, SPEEDEMON species assignment\n")
        f.write("Layout: tightly aligned side-by-side, shared legends left of panel 1\n")
        f.write(f"Extent: lon {extent[0]:.3f} to {extent[1]:.3f}; lat {extent[2]:.3f} to {extent[3]:.3f}\n")
        f.write(f"fastStructure best: K={best_k}, rep={best_rep}, file={fs_qfile}\n")
        f.write(f"Admixture clusters: {len(admix_cols)}\n")
        f.write(f"PopCluster clusters: {len(popcl_cols)}; best_k file value={popcl_best_k}\n")
        f.write(f"SPEEDEMON assignment groups: {len(speede_species)} from {SPEEDEMON_ASSIGN}\n")
        f.write(f"Output PNG: {OUT_PNG}\n")
        f.write(f"Output PDF: {OUT_PDF}\n")


if __name__ == "__main__":
    main()
