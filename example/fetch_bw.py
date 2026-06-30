#!/usr/bin/env python3
"""Download ENCODE H3K27ac BigWig files for the locus_plot example.

Queries the ENCODE REST API so results always reflect the latest
processed files for GM12878 and K562 cell lines (hg38, fold-change
over control).

Usage (from repo root):
    python example/fetch_bw.py
"""
import json, sys, urllib.request, urllib.parse, os

API    = "https://www.encodeproject.org"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUTDIR, exist_ok=True)

TARGETS = [
    ("GM12878", "GM12878_H3K27ac_fc.bw"),
    ("K562",    "K562_H3K27ac_fc.bw"),
]


def query_encode(biosample):
    params = {
        "type": "File",
        "biosample_ontology.term_name": biosample,
        "assay_title": "Histone ChIP-seq",
        "target.label": "H3K27ac",
        "file_format": "bigWig",
        "output_type": "fold change over control",
        "assembly": "GRCh38",
        "status": "released",
        "limit": 1,
        "format": "json",
    }
    query = "&".join(f"{k}={urllib.parse.quote(str(v))}" for k, v in params.items())
    url = f"{API}/search/?{query}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        data = json.load(r)
    hits = data.get("@graph", [])
    if not hits:
        raise RuntimeError(f"No ENCODE files found for biosample '{biosample}'")
    return API + hits[0]["href"]


def download(url, dest):
    print(f"  → {dest}")
    def _progress(blocks, block_size, total):
        done = blocks * block_size
        pct  = min(100, done * 100 // total) if total > 0 else 0
        print(f"\r  {pct:3d}%  {done / 1e6:.1f} / {total / 1e6:.1f} MB",
              end="", flush=True)
    urllib.request.urlretrieve(url, dest, reporthook=_progress)
    print()


for biosample, fname in TARGETS:
    dest = os.path.join(OUTDIR, fname)
    if os.path.exists(dest) and os.path.getsize(dest) > 1_000_000:
        print(f"[skip] {fname} already present  "
              f"({os.path.getsize(dest) // 1_000_000} MB)")
        continue
    print(f"\n[fetch] {biosample} H3K27ac  fold-change BigWig  hg38")
    try:
        url = query_encode(biosample)
        print(f"  {url}")
        download(url, dest)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

print(f"\nDone.  Files are in {OUTDIR}/")
print("\nGenerate the figure with:")
print("  python locus_plot.py --region chr8:127700000-128050000 "
      "--config example/tracks.ini --out example/output.pdf")
