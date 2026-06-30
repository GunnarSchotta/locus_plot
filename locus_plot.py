#!/usr/bin/env python3
"""
locus_plot.py  —  publication-quality genome browser track figures

Produces IGV-style stacked track plots from BigWig and BED files,
configured via an INI file.

Usage:
    python locus_plot.py --region chr17:18530000-18590000 \\
                         --config tracks.ini \\
                         --out figure.pdf \\
                         [--width 8] [--height-per-unit 0.8] [--dpi 200]

Track types (set with  type = ...):
    bigwig   — filled area signal from a BigWig file
    bed      — rectangular feature annotations from a BED file
    genes    — gene/transcript structures from a BED12 file
    ticks    — vertical tick marks for point features (BED, uses col 2 as pos)

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
"""

import argparse
import configparser
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pyBigWig
import os

# ── Constants ────────────────────────────────────────────────────────────────
NBINS           = 2000
DEFAULT_BW_COLOR  = "#333333"
DEFAULT_BED_COLOR = "#888888"
DEFAULT_GENE_COLOR = "#2244bb"
HIGHLIGHT_COLOR = "#e8e8e8"
COORD_AXIS_HEIGHT = 0.20   # height units for the coordinate axis


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

def draw_bigwig(ax, track, chrom, start, end):
    color     = track.get("color", DEFAULT_BW_COLOR)
    neg_color = track.get("neg_color", "#d4d4d4")
    vals  = fetch_bigwig(track["file"], chrom, start, end)
    x     = np.linspace(start, end, len(vals))
    v_fill = np.where(np.isnan(vals), 0, vals)

    v_pos = np.where(v_fill > 0, v_fill, 0.0)
    v_neg = np.where(v_fill < 0, v_fill, 0.0)
    ax.fill_between(x, 0, v_pos, color=color,     alpha=0.9, lw=0)
    ax.fill_between(x, 0, v_neg, color=neg_color, alpha=0.9, lw=0)
    ax.axhline(0, color="black", lw=0.4)

    if "ylim" in track:
        ymin, ymax = [float(v) for v in track["ylim"].split(",")]
    else:
        ymax = float(np.nanmax(v_fill)) if np.any(~np.isnan(vals)) else 1.0
        ymin = min(0.0, float(np.nanmin(v_fill)))
    ax.set_ylim(ymin, ymax * 1.05)
    ax.set_yticks([])

    # Y-range text: "[0 – max]" in top-left corner of track (IGV style)
    range_str = f"[{ymin:.0f} – {ymax:.1f}]"
    ax.text(0.005, 0.97, range_str,
            transform=ax.transAxes, fontsize=5.5, va="top", ha="left",
            color="grey", clip_on=True)


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
                                            lw=0.4, mutation_scale=4), zorder=3)

        # Feature name
        if show_names and name_col is not None and len(cols) > name_col:
            name = cols[name_col]
            if (fe - fs) > 0.005 * span:
                ax.text((fs + fe) / 2, y_mid - feat_h / 2 - 0.06, name,
                        ha="center", va="top", fontsize=5.5,
                        color="#444444", clip_on=True, style="italic")


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

        # Gene body line
        ax.plot([tx_s, tx_e], [y_mid, y_mid], color=color, lw=0.8, zorder=1)

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
                                        lw=0.5, mutation_scale=6), zorder=3)

        # Gene name (italic, below — clip_on=False lets names extend into gap)
        if name:
            ax.text((tx_s + tx_e) / 2, y_mid - exon_h / 2 - 0.08, name,
                    ha="center", va="top", fontsize=6.5,
                    fontstyle="italic", color=color, clip_on=False)


def draw_ticks(ax, track, chrom, start, end):
    color = track.get("color", "#cc6600")
    rows  = parse_bed_region(track["file"], chrom, start, end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for cols in rows:
        xpos = (int(cols[1]) + int(cols[2])) / 2
        ax.axvline(xpos, color=color, lw=1.0, ymin=0.1, ymax=0.9, zorder=2)


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
            color="black", lw=1.5, clip_on=False, solid_capstyle="butt")
    label = f"{bar_bp // 1000} kb" if bar_bp >= 1000 else f"{bar_bp} bp"
    ax.text((x_start + x_end) / 2, y_text, label,
            transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=7, color="black")


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
        fontsize=6.5)
    ax.xaxis.set_tick_params(length=3, width=0.5)
    if position == "top":
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")
        ax.spines[["bottom", "right", "left"]].set_visible(False)
        ax.spines["top"].set_linewidth(0.4)
    else:
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.4)
    ax.set_xlabel(chrom, fontsize=9, fontweight="bold", labelpad=2)


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
                 ha="center", va="center", fontsize=7.5, fontweight="bold",
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
                    help="Figure width in inches (default 8)")
    ap.add_argument("--height-per-unit", type=float, default=0.8, dest="hpu",
                    help="Inches per height unit (default 0.8)")
    ap.add_argument("--dpi",     type=int,   default=200)
    args = ap.parse_args()

    matplotlib.rcParams.update({
        "font.family":      "sans-serif",
        "font.sans-serif":  ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":        8,
        "pdf.fonttype":     42,
        "svg.fonttype":     "none",
    })

    chrom, start, end = parse_region(args.region)
    glb, tracks = load_config(args.config)

    # Figure height: coord axis on top, then tracks
    heights = [COORD_AXIS_HEIGHT] + [t["height"] for t in tracks]
    fig_h = sum(heights) * args.hpu
    fig_h = max(fig_h, 2.0)

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
                fontsize=7, ha="right", va="center",
                fontstyle="italic" if italic else "normal",
                clip_on=False)

    # Highlight column across all tracks
    add_highlight(list(axes), start, end, glb.get("highlight", ""))

    # Title goes on the coord axis at the top
    if glb.get("title"):
        axes[0].set_title(glb["title"], fontsize=9, fontweight="bold", pad=4)

    plt.subplots_adjust(hspace=0.02, left=0.25, right=0.97,
                        top=0.97, bottom=0.05)

    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {args.out}")


if __name__ == "__main__":
    main()
