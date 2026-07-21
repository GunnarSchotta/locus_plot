#!/usr/bin/env python3
"""
locus_plot.py  —  publication-quality genome browser track figures

Produces IGV-style stacked track plots from BigWig and BED files,
configured via an INI file.

Usage:
    python locus_plot.py --region chr17:18530000-18590000 \\
                         --config tracks.ini \\
                         --out figure.pdf \\
                         [--width 8] [--height-per-unit auto] [--dpi 200]

--width also scales font sizes, line widths, and track heights up or down
together (down to a minimum readable size below ~6.5in); see FONT_SIZES /
set_scale() / auto_hpu() below.

Track types (set with  type = ...):
    bigwig   — filled area signal from a BigWig file
    bed      — rectangular feature annotations from a BED file
    genes    — gene/transcript structures from a BED12 file
    ticks    — vertical tick marks for point features (BED, uses col 2 as pos)
    dynseq   — per-base score bigWig rendered as scaled/colored ACGT letters
               (sequence looked up from a reference FASTA), falling back to
               a plain filled signal view when zoomed out too far for
               individual bases to be legible

Global options in a [_global] section:
    highlight  = 18558000-18567000   # grey shaded column (region coords)
    scalebar   = 5000                # scale bar length in bp (drawn top-right)
    title      =                     # optional figure title

Track options (all optional except file and type):
    file        = /path/to/file.bw or .bed
    type        = bigwig | bed | genes | ticks
    color       = #2166ac            # single color
    height      = 1.5                # relative track height (default 1)
    group       = H3K9me3            # group label shown further left, spanning
                                     # all consecutive tracks with the same group
    label       = ES                 # individual track label (default = section name)
    ylim        = -1,5               # BigWig y-axis limits; omit for auto
    color_col   = 13                 # 1-based BED column for per-feature colors
    color_map   = ES_only:#2166ac,XEN_only:#d6604d,shared:#7B4F9E,negative:#cccccc
    name_col    = 5                  # 1-based column for feature name labels
    strand_col  = 7                  # 1-based column for strand arrows (bed)
    show_names  = true               # whether to draw name labels (default true)
    fasta       = /path/to/genome.fa # reference sequence for dynseq tracks
    a_color     = #109648            # per-base letter color overrides (dynseq)
    c_color     = #255C99
    g_color     = #F7B32B
    t_color     = #D62839
"""

import argparse
import configparser
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
from matplotlib.transforms import Bbox, Affine2D
import pyBigWig
import pyfaidx
import os

# ── Constants ────────────────────────────────────────────────────────────────
NBINS           = 2000
DEFAULT_BW_COLOR  = "#333333"
DEFAULT_BED_COLOR = "#888888"
DEFAULT_GENE_COLOR = "#2244bb"
HIGHLIGHT_COLOR = "#e8e8e8"
COORD_AXIS_HEIGHT = 0.20   # height units for the coordinate axis
NUCLEOTIDE_COLORS = {"A": "#109648", "C": "#255C99", "G": "#F7B32B", "T": "#D62839"}
MIN_LETTER_WIDTH_IN = 0.05   # below this per-base width, letters aren't legible --
                             # dynseq falls back to a filled signal view, like a
                             # genome browser does when zoomed out
AXES_LEFT, AXES_RIGHT = 0.25, 0.97   # must match subplots_adjust() margins in main()

# ── Size scaling ─────────────────────────────────────────────────────────────
# Font sizes / line widths below are tuned for a figure at REF_WIDTH inches
# (the historical default). Wider or narrower figures scale everything up or
# down together via --width -- e.g. a poster-size figure gets bigger text,
# a publication column gets smaller text. Each font role has a floor below
# which it stops shrinking: pushing --width narrower than that trades away
# layout (labels may overlap) rather than producing illegible type.
REF_WIDTH = 8.0

FONT_SIZES = {
    "tiny":  5.5,   # y-range readout, bed/gene secondary annotations
    "small": 6.5,   # coord-axis tick labels, gene name labels
    "label": 7.0,   # track labels, scale bar label
    "group": 7.5,   # group labels
    "base":  8.0,   # default body text (rcParams)
    "title": 9.0,   # chromosome label, figure title
}
FONT_FLOORS = {
    "tiny":  4.5,
    "small": 5.0,
    "label": 5.5,
    "group": 5.5,
    "base":  6.0,
    "title": 7.0,
}
LINE_WIDTHS = {
    "hair":  0.4,   # zero line, strand arrows, spines
    "thin":  0.5,   # coord-axis tick width, gene-body arrows
    "med":   0.8,   # gene body line
    "thick": 1.0,   # point-feature ticks
    "bar":   1.5,   # scale bar
}
LINE_WIDTH_FLOOR = 0.3
ARROW_SCALES = {"strand": 4, "gene": 6}
ARROW_SCALE_FLOOR = 2.5

DEFAULT_HPU = 0.8   # default inches per height unit, at REF_WIDTH
HPU_FLOOR_RATIO = 0.6   # track height never auto-shrinks below 60% of default

_SCALE = 1.0  # set once in main() from --width, via set_scale()


def set_scale(width):
    """Derive the global font/line-width scale factor from figure width.

    Scales up as well as down around REF_WIDTH, so a wider figure gets
    bigger, more legible text instead of staying frozen at the reference
    size. Warns when a role's floor kicks in on the way down, since that's
    the point where shrinking further starts costing legibility.
    """
    global _SCALE
    _SCALE = width / REF_WIDTH
    clamped = [r for r in FONT_SIZES if FONT_SIZES[r] * _SCALE < FONT_FLOORS[r]]
    if _SCALE < HPU_FLOOR_RATIO:
        clamped = clamped + ["track height"]
    if clamped:
        print(f"WARNING: --width {width:.2f}in scales {', '.join(clamped)} "
              f"below their floor; holding at minimum readable size "
              f"-- labels may overlap.", file=sys.stderr)
    return _SCALE


def auto_hpu():
    """Inches per height unit, scaled with width like everything else (up or
    down), down to HPU_FLOOR_RATIO of the default on the narrow end -- so a
    narrow publication figure also gets shorter, instead of staying tall and
    square. Overridden outright by an explicit --height-per-unit."""
    return DEFAULT_HPU * max(_SCALE, HPU_FLOOR_RATIO)


def FS(role):
    return max(FONT_SIZES[role] * _SCALE, FONT_FLOORS[role])


def LW(role):
    return max(LINE_WIDTHS[role] * _SCALE, LINE_WIDTH_FLOOR)


def AS(role):
    return max(ARROW_SCALES[role] * _SCALE, ARROW_SCALE_FLOOR)


LABEL_LINE_HEIGHT = 1.4   # rough line-height multiplier over font point size

def annotation_font_role(track):
    """Font role of the per-feature name labels a track draws inside its own
    panel (gene names, bed feature names) -- distinct from the track's own
    row label to its left. Returns None if the track draws no such labels,
    in which case its configured height needs no adjustment."""
    ttype = track.get("type", "bigwig").lower()
    if ttype == "genes":
        return "small"
    if ttype == "bed":
        show_names = track.get("show_names", "true").lower() != "false"
        if show_names and "name_col" in track:
            return "tiny"
    return None


def track_height_units(track, hpu):
    """Configured height, plus extra room (converted from the scaled label
    font size into height units) if the track draws names inside its panel
    -- otherwise the configured height is used as-is."""
    h = track["height"]
    role = annotation_font_role(track)
    if role is not None:
        pad_inches = (FS(role) / 72.0) * LABEL_LINE_HEIGHT
        h += pad_inches / hpu
    return h


# ── Parsing helpers ──────────────────────────────────────────────────────────

def parse_region(s):
    try:
        chrom, coords = s.split(":")
        a, b = coords.replace(",", "").split("-")
        return chrom, int(a), int(b)
    except Exception:
        sys.exit(f"Cannot parse region '{s}'. Expected chrN:start-end")


def parse_color_map(s):
    d = {}
    for pair in s.split(","):
        if ":" in pair:
            k, v = pair.split(":", 1)
            d[k.strip()] = v.strip()
    return d


def load_config(path):
    cfg = configparser.ConfigParser()
    cfg.read(path)
    glb = dict(cfg["_global"]) if "_global" in cfg else {}
    tracks = []
    for sec in cfg.sections():
        if sec.startswith("_"):
            continue
        t = dict(cfg[sec])
        t.setdefault("label",  sec)
        t.setdefault("type",   "bigwig")
        t.setdefault("height", "1")
        t["height"] = float(t["height"])
        tracks.append(t)
    return glb, tracks


# ── Data fetchers ────────────────────────────────────────────────────────────

def fetch_bigwig(path, chrom, start, end, nbins=NBINS):
    bw = pyBigWig.open(path)
    raw = bw.stats(chrom, start, end, type="mean", nBins=nbins)
    bw.close()
    return np.array([v if v is not None else np.nan for v in raw], dtype=float)


def fetch_bigwig_values(path, chrom, start, end):
    """Per-base-pair values (one per position), unlike fetch_bigwig's binned means."""
    bw = pyBigWig.open(path)
    vals = np.array(bw.values(chrom, start, end), dtype=float)
    bw.close()
    return vals


def fetch_sequence(path, chrom, start, end):
    fasta = pyfaidx.Fasta(path)
    seq = str(fasta[chrom][start:end]).upper()
    fasta.close()
    return seq


def parse_bed_region(path, chrom, start, end):
    rows = []
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(("#", "track", "browser")) or not line.strip():
                    continue
                c = line.rstrip("\n").split("\t")
                if c[0] != chrom:
                    continue
                if int(c[2]) < start or int(c[1]) > end:
                    continue
                rows.append(c)
    except FileNotFoundError:
        print(f"WARNING: file not found: {path}", file=sys.stderr)
    return rows


# ── Track renderers ──────────────────────────────────────────────────────────

def _draw_filled_signal(ax, x, v_fill, color, neg_color, ylim=None):
    """Filled positive/negative area plot with zero line, y-limits, and the
    IGV-style '[min – max]' range readout. Shared by the bigwig track and by
    dynseq's zoomed-out fallback view."""
    v_pos = np.where(v_fill > 0, v_fill, 0.0)
    v_neg = np.where(v_fill < 0, v_fill, 0.0)
    ax.fill_between(x, 0, v_pos, color=color,     alpha=0.9, lw=0)
    ax.fill_between(x, 0, v_neg, color=neg_color, alpha=0.9, lw=0)
    ax.axhline(0, color="black", lw=LW("hair"))

    if ylim is not None:
        ymin, ymax = ylim
    else:
        ymax = float(np.nanmax(v_fill)) if np.any(~np.isnan(v_fill)) else 1.0
        ymin = min(0.0, float(np.nanmin(v_fill)))
    ax.set_ylim(ymin, ymax * 1.05)
    ax.set_yticks([])

    range_str = f"[{ymin:.0f} – {ymax:.1f}]"
    ax.text(0.005, 0.97, range_str,
            transform=ax.transAxes, fontsize=FS("tiny"), va="top", ha="left",
            color="grey", clip_on=True)


def draw_bigwig(ax, track, chrom, start, end):
    color     = track.get("color", DEFAULT_BW_COLOR)
    neg_color = track.get("neg_color", "#d4d4d4")
    vals  = fetch_bigwig(track["file"], chrom, start, end)
    x     = np.linspace(start, end, len(vals))
    v_fill = np.where(np.isnan(vals), 0, vals)
    ylim = tuple(float(v) for v in track["ylim"].split(",")) if "ylim" in track else None
    _draw_filled_signal(ax, x, v_fill, color, neg_color, ylim)


def draw_bed(ax, track, chrom, start, end):
    color_default = track.get("color", DEFAULT_BED_COLOR)
    color_col  = int(track["color_col"])  - 1 if "color_col"  in track else None
    color_map  = parse_color_map(track["color_map"]) if "color_map" in track else {}
    name_col   = int(track["name_col"])   - 1 if "name_col"   in track else None
    strand_col = int(track["strand_col"]) - 1 if "strand_col" in track else None
    show_names = track.get("show_names", "true").lower() != "false"
    feat_h = 0.55
    y_mid  = 0.5

    rows = parse_bed_region(track["file"], chrom, start, end)
    span = end - start
    ax.set_ylim(0, 1)
    ax.set_yticks([])

    # Collect eligible feature-name labels and draw them after the main loop
    # (below), so densely-packed features (e.g. ERV annotation) can have
    # adjacent, colliding labels merged into one "name1/name2" label instead
    # of drawing illegible overlapping text -- same fix as draw_genes.
    label_items = []

    for cols in rows:
        fs = max(int(cols[1]), start)
        fe = min(int(cols[2]), end)
        if fe <= fs:
            continue
        color = color_default
        if color_col is not None and len(cols) > color_col:
            color = color_map.get(cols[color_col], color_default)

        rect = mpatches.FancyBboxPatch(
            (fs, y_mid - feat_h / 2), fe - fs, feat_h,
            boxstyle="square,pad=0",
            facecolor=color, edgecolor="none", zorder=2)
        ax.add_patch(rect)

        # Strand arrows
        strand = None
        if strand_col is not None and len(cols) > strand_col:
            strand = cols[strand_col]
        if strand in ("+", "-"):
            n = max(1, int((fe - fs) / (span / 25)))
            for i in range(1, n + 1):
                xp = fs + i * (fe - fs) / (n + 1)
                dx = 0.004 * span * (1 if strand == "+" else -1)
                ax.annotate("", xy=(xp + dx, y_mid), xytext=(xp, y_mid),
                            arrowprops=dict(arrowstyle="-|>", color="white",
                                            lw=LW("hair"), mutation_scale=AS("strand")), zorder=3)

        # Feature name
        if show_names and name_col is not None and len(cols) > name_col:
            name = cols[name_col]
            if (fe - fs) > 0.005 * span:
                label_items.append((name, (fs, fe)))

    for name, (name_s, name_e) in _merge_overlapping_labels(label_items, span, FS("tiny")):
        ax.text((name_s + name_e) / 2, y_mid - feat_h / 2 - 0.06, name,
                ha="center", va="top", fontsize=FS("tiny"),
                color="#444444", clip_on=True, style="italic")


_GENE_LABEL_FONT = FontProperties(family="sans-serif", style="italic")


def _label_width_bp(text, fontsize_pt, bp_per_inch):
    """Rendered width of a gene-name label, in bp at the current --width/zoom,
    via TextPath's exact font metrics (coordinates in points at the given
    size) -- no canvas draw needed, so this can run before layout/subplots_adjust."""
    if not text:
        return 0.0
    width_pt = TextPath((0, 0), text, size=fontsize_pt, prop=_GENE_LABEL_FONT).get_extents().width
    return (width_pt / 72.0) * bp_per_inch


def _merge_overlapping_labels(items, span, fontsize_pt):
    """Merge adjacent (name, (start, end)) labels whose rendered text would
    collide on screen into a single combined label (e.g. "Sprr2a1/Sprr2a2")
    instead of drawing two overlapping/garbled labels -- neither the genes
    track nor the bed track stack labels onto a second row, so two features
    close enough together previously collided illegibly. Shared by
    draw_genes (gene names) and draw_bed (feature names, e.g. dense ERV
    annotation). Returns a list of (name, (start, end)), sorted by position."""
    items = sorted(items, key=lambda kv: (kv[1][0] + kv[1][1]) / 2)
    if len(items) < 2:
        return items

    fig_width_in = _SCALE * REF_WIDTH
    ax_width_in = fig_width_in * (AXES_RIGHT - AXES_LEFT)
    bp_per_inch = span / ax_width_in if ax_width_in > 0 else 0.0
    pad_bp = (4.0 / 72.0) * bp_per_inch  # ~4pt breathing room between labels

    merged = [items[0]]
    for name, (s, e) in items[1:]:
        prev_name, (prev_s, prev_e) = merged[-1]
        prev_center = (prev_s + prev_e) / 2
        cur_center = (s + e) / 2
        gap_needed = (_label_width_bp(prev_name, fontsize_pt, bp_per_inch) / 2
                      + _label_width_bp(name, fontsize_pt, bp_per_inch) / 2
                      + pad_bp)
        if cur_center - prev_center < gap_needed:
            merged[-1] = (f"{prev_name}/{name}", (min(prev_s, s), max(prev_e, e)))
        else:
            merged.append((name, (s, e)))
    return merged


def draw_genes(ax, track, chrom, start, end):
    color      = track.get("color", DEFAULT_GENE_COLOR)
    name_sep   = track.get("name_sep", "|")
    name_field_raw = track.get("name_field", None)
    name_field = int(name_field_raw) if name_field_raw is not None else None
    rows  = parse_bed_region(track["file"], chrom, start, end)
    span  = end - start
    exon_h, y_mid = 0.45, 0.55

    ax.set_ylim(0, 1)
    ax.set_yticks([])

    # Multiple transcripts of the same gene often overlap almost entirely in
    # view; label each unique gene name once (centered on the union of its
    # transcripts) instead of once per transcript, which produced illegible
    # stacked/overlapping text.
    name_spans = {}

    for cols in rows:
        tx_s = max(int(cols[1]), start)
        tx_e = min(int(cols[2]), end)
        if tx_e <= tx_s:
            continue
        strand = cols[5] if len(cols) > 5 else "+"
        name   = cols[3] if len(cols) > 3 else ""
        if name_field is not None and name_sep in name:
            parts = name.split(name_sep)
            name = parts[name_field] if name_field < len(parts) else name

        if name:
            prev = name_spans.get(name)
            name_spans[name] = (tx_s, tx_e) if prev is None else \
                (min(prev[0], tx_s), max(prev[1], tx_e))

        # Gene body line
        ax.plot([tx_s, tx_e], [y_mid, y_mid], color=color, lw=LW("med"), zorder=1)

        # Exon blocks (BED12: cols 10–11)
        if len(cols) >= 12 and cols[9].strip():
            sizes  = [int(x) for x in cols[10].rstrip(",").split(",") if x]
            starts = [int(x) for x in cols[11].rstrip(",").split(",") if x]
            for sz, st in zip(sizes, starts):
                ex_s = max(int(cols[1]) + st, start)
                ex_e = min(int(cols[1]) + st + sz, end)
                if ex_e > ex_s:
                    ax.add_patch(mpatches.Rectangle(
                        (ex_s, y_mid - exon_h / 2), ex_e - ex_s, exon_h,
                        facecolor=color, edgecolor="none", zorder=2))
        else:
            ax.add_patch(mpatches.Rectangle(
                (tx_s, y_mid - exon_h / 2), tx_e - tx_s, exon_h,
                facecolor=color, edgecolor="none", zorder=2))

        # Directional arrows along gene body
        n = max(1, int((tx_e - tx_s) / (span / 20)))
        for i in range(1, n + 1):
            xp = tx_s + i * (tx_e - tx_s) / (n + 1)
            dx = 0.0035 * span * (1 if strand == "+" else -1)
            ax.annotate("", xy=(xp + dx, y_mid), xytext=(xp, y_mid),
                        arrowprops=dict(arrowstyle="-|>", color=color,
                                        lw=LW("thin"), mutation_scale=AS("gene")), zorder=3)

    # Gene name (italic, below — clip_on=False lets names extend into gap)
    # One label per unique name, centered on the union of its transcripts;
    # adjacent names whose labels would collide on screen are merged into a
    # single "name1/name2" label (see _merge_overlapping_labels).
    for name, (name_s, name_e) in _merge_overlapping_labels(list(name_spans.items()), span, FS("small")):
        ax.text((name_s + name_e) / 2, y_mid - exon_h / 2 - 0.08, name,
                ha="center", va="top", fontsize=FS("small"),
                fontstyle="italic", color=color, clip_on=False)


_LETTER_FONT = FontProperties(family="monospace", weight="bold")
_LETTER_PATHS = {}


def _letter_path(letter):
    """TextPath + bounding box for a glyph, normalized once and cached."""
    if letter not in _LETTER_PATHS:
        tp = TextPath((0, 0), letter, size=1, prop=_LETTER_FONT)
        _LETTER_PATHS[letter] = (tp, tp.get_extents())
    return _LETTER_PATHS[letter]


def draw_letter(ax, letter, x_center, y_base, height, width, color, flip=False):
    """Draw one glyph scaled to fill a (width x height) box anchored at
    (x_center, y_base) -- height independent of width, unlike fontsize-based
    text, since dynseq letter height encodes a score. flip=True mirrors the
    glyph below y_base (dynseq's convention for negative scores)."""
    if height <= 0:
        return
    tp, bbox = _letter_path(letter)
    bbox_w = bbox.width or 1.0
    bbox_h = bbox.height or 1.0
    t = (Affine2D()
         .translate(-bbox.x0 - bbox_w / 2, -bbox.y0)
         .scale(width / bbox_w, (height / bbox_h) * (-1 if flip else 1))
         .translate(x_center, y_base))
    ax.add_patch(mpatches.PathPatch(tp.transformed(t), facecolor=color, edgecolor="none", lw=0, zorder=2))


def draw_dynseq(ax, track, chrom, start, end):
    color     = track.get("color", DEFAULT_BW_COLOR)
    neg_color = track.get("neg_color", "#d4d4d4")
    ylim = tuple(float(v) for v in track["ylim"].split(",")) if "ylim" in track else None

    span = end - start
    fig_width_in = _SCALE * REF_WIDTH
    px_per_bp_in = fig_width_in * (AXES_RIGHT - AXES_LEFT) / span

    if px_per_bp_in < MIN_LETTER_WIDTH_IN:
        # Too zoomed out for legible letters -- fall back to a plain filled
        # signal view, same as a genome browser does at low zoom.
        vals = fetch_bigwig(track["file"], chrom, start, end)
        x = np.linspace(start, end, len(vals))
        v_fill = np.where(np.isnan(vals), 0, vals)
        _draw_filled_signal(ax, x, v_fill, color, neg_color, ylim)
        return

    vals = fetch_bigwig_values(track["file"], chrom, start, end)
    seq  = fetch_sequence(track["fasta"], chrom, start, end)
    letter_colors = {b: track.get(f"{b.lower()}_color", NUCLEOTIDE_COLORS[b]) for b in "ACGT"}

    if ylim is not None:
        ymin, ymax = ylim
    else:
        finite = vals[np.isfinite(vals)]
        ymax = float(finite.max()) if finite.size and finite.max() > 0 else 1.0
        ymin = float(finite.min()) if finite.size and finite.min() < 0 else 0.0

    ax.axhline(0, color="black", lw=LW("hair"))
    ax.set_ylim(ymin * 1.05 if ymin < 0 else 0, ymax * 1.05)
    ax.set_yticks([])

    for i, base in enumerate(seq):
        score = vals[i] if i < len(vals) else np.nan
        if not np.isfinite(score) or score == 0 or base not in letter_colors:
            continue
        draw_letter(ax, base, start + i + 0.5, 0, abs(score), 0.9,
                    letter_colors[base], flip=score < 0)


def draw_ticks(ax, track, chrom, start, end):
    color = track.get("color", "#cc6600")
    rows  = parse_bed_region(track["file"], chrom, start, end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for cols in rows:
        xpos = (int(cols[1]) + int(cols[2])) / 2
        ax.axvline(xpos, color=color, lw=LW("thick"), ymin=0.1, ymax=0.9, zorder=2)


# ── Decoration helpers ───────────────────────────────────────────────────────

def apply_common_style(ax, start, end):
    ax.set_xlim(start, end)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.spines[["top", "right", "bottom", "left"]].set_visible(False)


def add_highlight(axes_list, start, end, region_str):
    if not region_str:
        return
    try:
        a, b = region_str.replace(",", "").split("-")
        hs, he = int(a), int(b)
    except Exception:
        return
    for ax in axes_list:
        ax.axvspan(hs, he, color=HIGHLIGHT_COLOR, zorder=0, lw=0)


def add_scalebar(ax, start, end, bar_bp, position="bottom"):
    """Draw a scale bar (black line + label) in the coord axis strip."""
    if not bar_bp:
        return
    bar_bp = int(bar_bp)
    x_end   = end - 0.02 * (end - start)
    x_start = x_end - bar_bp
    y_pos   = 0.45 if position == "top" else 0.88
    y_text  = y_pos + (0.25 if position == "top" else 0.06)
    ax.plot([x_start, x_end], [y_pos, y_pos],
            transform=ax.get_xaxis_transform(),
            color="black", lw=LW("bar"), clip_on=False, solid_capstyle="butt")
    label = f"{bar_bp // 1000} kb" if bar_bp >= 1000 else f"{bar_bp} bp"
    ax.text((x_start + x_end) / 2, y_text, label,
            transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=FS("label"), color="black")


def add_coord_axis(ax, chrom, start, end, position="bottom"):
    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    span = end - start
    for interval in [5e6, 2e6, 1e6, 5e5, 2e5, 1e5, 5e4, 2e4, 1e4, 5e3, 2e3, 1e3]:
        if span / interval >= 3:
            break
    ticks = np.arange(int(np.ceil(start / interval)) * interval,
                      end + 1, interval)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        [f"{t/1e6:.2f} Mb" if t >= 1e6 else f"{int(t/1e3)} kb" for t in ticks],
        fontsize=FS("small"))
    ax.xaxis.set_tick_params(length=max(3 * _SCALE, 2), width=LW("thin"))
    if position == "top":
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.spines[["bottom", "right", "left"]].set_visible(False)
        ax.spines["top"].set_linewidth(LW("hair"))
    else:
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(LW("hair"))
    ax.set_xlabel(chrom, fontsize=FS("title"), fontweight="bold", labelpad=2)


def add_group_labels(fig, axes, tracks, left_x=0.01):
    """Draw group labels in figure-space, spanning all tracks in each group."""
    # Flush layout so ax.get_position() is accurate
    fig.canvas.draw()
    groups = {}
    for ax, t in zip(axes, tracks):
        g = t.get("group", "")
        if not g:
            continue
        if g not in groups:
            groups[g] = []
        groups[g].append(ax)

    for gname, gaxes in groups.items():
        pos_top = gaxes[0].get_position()
        pos_bot = gaxes[-1].get_position()
        y_top = pos_top.y0 + pos_top.height
        y_bot = pos_bot.y0
        y_mid = (y_top + y_bot) / 2
        fig.text(left_x, y_mid, gname,
                 ha="center", va="center", fontsize=FS("group"), fontweight="bold",
                 rotation=90, transform=fig.transFigure)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Publication-quality genome browser track figures")
    ap.add_argument("--region",  required=True,
                    help="Genomic region, e.g. chr17:18530000-18590000")
    ap.add_argument("--config",  required=True,
                    help="Track config INI file")
    ap.add_argument("--out",     required=True,
                    help="Output file (.pdf, .svg, .png)")
    ap.add_argument("--width",   type=float, default=8,
                    help="Figure width in inches (default 8). Also scales "
                         "font sizes and line widths up or down with figure "
                         "size -- larger for a poster, smaller for a "
                         "publication column; below ~6.5in text holds at a "
                         "minimum readable size and labels may start to "
                         "overlap instead of shrinking further (a warning "
                         "is printed when this happens)")
    ap.add_argument("--height-per-unit", type=float, default=None, dest="hpu",
                    help="Inches per height unit (default: auto, scales "
                         "with --width like text does, down to 60%% of the "
                         "0.8in default on the narrow end; pass a value to "
                         "fix it explicitly)")
    ap.add_argument("--dpi",     type=int,   default=200)
    args = ap.parse_args()

    set_scale(args.width)

    matplotlib.rcParams.update({
        "font.family":      "sans-serif",
        "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":        FS("base"),
        "pdf.fonttype":     42,
        "svg.fonttype":     "none",
    })

    chrom, start, end = parse_region(args.region)
    glb, tracks = load_config(args.config)

    # Figure height: coord axis on top, then tracks. hpu scales down with
    # --width by default (see auto_hpu) so a narrower figure also gets
    # shorter instead of staying tall and square; pass --height-per-unit
    # explicitly to fix it regardless of width.
    hpu = args.hpu if args.hpu is not None else auto_hpu()
    heights = [COORD_AXIS_HEIGHT] + [track_height_units(t, hpu) for t in tracks]
    fig_h = sum(heights) * hpu
    fig_h = max(fig_h, 2.0 * max(_SCALE, HPU_FLOOR_RATIO))

    fig, axes = plt.subplots(
        len(tracks) + 1, 1,
        figsize=(args.width, fig_h),
        gridspec_kw={"height_ratios": heights},
    )

    # Coordinate axis at top
    add_coord_axis(axes[0], chrom, start, end, position="top")
    # Scale bar in the first signal track (top-right corner)
    if len(axes) > 1:
        add_scalebar(axes[1], start, end, glb.get("scalebar", ""))

    # Render tracks (axes[1:] maps to tracks)
    dispatchers = {
        "bigwig": draw_bigwig,
        "bed":    draw_bed,
        "genes":  draw_genes,
        "ticks":  draw_ticks,
        "dynseq": draw_dynseq,
    }
    for ax, track in zip(axes[1:], tracks):
        ttype = track.get("type", "bigwig").lower()
        fn = dispatchers.get(ttype)
        if fn is None:
            print(f"Unknown track type '{ttype}', skipping.", file=sys.stderr)
            ax.set_visible(False)
            continue
        fn(ax, track, chrom, start, end)
        apply_common_style(ax, start, end)

        # Individual track label — placed at visual center of track
        label = track.get("label", "")
        italic = track.get("italic_label", "false").lower() == "true"
        ax.text(-0.02, 0.5, label,
                transform=ax.transAxes,
                fontsize=FS("label"), ha="right", va="center",
                fontstyle="italic" if italic else "normal",
                clip_on=False)

    # Highlight column across all tracks
    add_highlight(list(axes), start, end, glb.get("highlight", ""))

    # Title goes on the coord axis at the top
    if glb.get("title"):
        axes[0].set_title(glb["title"], fontsize=FS("title"), fontweight="bold", pad=4)

    plt.subplots_adjust(hspace=0.02, left=AXES_LEFT, right=AXES_RIGHT,
                        top=0.97, bottom=0.05)

    # Trim excess top/bottom whitespace automatically, but lock the saved
    # width to exactly --width -- a plain bbox_inches="tight" would instead
    # grow the canvas to fit any overflowing labels, silently ignoring
    # --width right when a narrow figure needs it most. Content that doesn't
    # fit at the requested width is clipped instead.
    fig.canvas.draw()
    tight = fig.get_tightbbox(fig.canvas.get_renderer())
    final_bbox = Bbox.from_bounds(0, tight.y0, args.width, tight.height)
    fig.savefig(args.out, dpi=args.dpi, bbox_inches=final_bbox)
    plt.close(fig)
    print(f"Saved → {args.out}")


if __name__ == "__main__":
    main()
