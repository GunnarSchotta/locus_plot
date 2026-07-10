#!/usr/bin/env python3
"""Download the hg38 chr8 reference sequence for the locus_plot dynseq example.

Fetches chr8.fa.gz from the UCSC goldenPath archive and decompresses it into
example/data/chr8.fa. dynseq tracks need the actual reference genome (or at
least the relevant chromosome) to look up the base letter at each position,
the same way a real analysis would use the lab's hg38 FASTA.

Usage (from repo root):
    python example/fetch_fasta.py
"""
import gzip
import os
import shutil
import urllib.request

URL    = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr8.fa.gz"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
GZ     = os.path.join(OUTDIR, "chr8.fa.gz")
FA     = os.path.join(OUTDIR, "chr8.fa")

os.makedirs(OUTDIR, exist_ok=True)

if os.path.exists(FA) and os.path.getsize(FA) > 100_000_000:
    print(f"[skip] {FA} already present ({os.path.getsize(FA)//1_000_000} MB)")
else:
    print(f"Downloading {URL}")

    def _progress(block_num, block_size, total_size):
        done = block_num * block_size
        pct = min(100, done * 100 // total_size) if total_size > 0 else 0
        print(f"\r  {pct:3d}%  {done/1e6:.1f} / {total_size/1e6:.1f} MB", end="", flush=True)

    urllib.request.urlretrieve(URL, GZ, reporthook=_progress)
    print()

    print(f"Decompressing -> {FA}")
    with gzip.open(GZ, "rb") as src, open(FA, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(GZ)

print(f"Done: {FA} ({os.path.getsize(FA)//1_000_000} MB)")
print("A .fai index will be created automatically the first time locus_plot.py reads it.")
