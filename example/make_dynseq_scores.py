#!/usr/bin/env python3
"""Generate toy per-base importance scores for the locus_plot dynseq example.

dynseq tracks are normally driven by real per-base model output (e.g.
DeepLIFT / in-silico-mutagenesis contribution scores from a BPNet/ChromBPNet-
style model). This repo doesn't ship a trained model, so this script instead
writes a SYNTHETIC score track -- a smooth bump centered on a real, annotated
CTCF site (chr8:127,729,800-127,730,050, see data/ctcf_sites.bed) plus a
little seeded noise -- purely to demonstrate what the dynseq track type looks
like. Swap `file` in tracks_dynseq.ini for your own per-base score BigWig to
visualize real scores.

Requires pyBigWig (pip install pyBigWig  or  conda install -c bioconda pybigwig).

Usage (from repo root):
    python example/make_dynseq_scores.py
"""
import os
import numpy as np
import pyBigWig

OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUTDIR, exist_ok=True)

CHROM = "chr8"
START, END = 127_729_875, 127_729_975   # 100 bp window inside the real CTCF_1 site
CHROM_SIZE = 145_138_636                # hg38 chr8 length

rng = np.random.default_rng(0)
x = np.arange(START, END)
center = (START + END) / 2
bump = 4.0 * np.exp(-((x - center) ** 2) / (2 * 15 ** 2))
noise = rng.normal(0, 0.3, size=len(x))
scores = (bump + noise).astype(np.float32)

out_path = os.path.join(OUTDIR, "ctcf_motif_scores.bw")
bw = pyBigWig.open(out_path, "w")
bw.addHeader([(CHROM, CHROM_SIZE)])
bw.addEntries(CHROM, int(START), values=scores, span=1, step=1)
bw.close()
print(f"Wrote {len(scores)} per-base toy scores -> {out_path}")
