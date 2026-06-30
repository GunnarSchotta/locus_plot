#!/usr/bin/env python3
"""
Fetch ENCODE H3K27ac peak calls for GM12878 and K562, filter to the
example region, classify each peak as GM12878_only / K562_only / shared,
and write example/data/peaks.bed.

Requires pyBigWig (pip install pyBigWig  or  conda install -c bioconda pybigwig).

Usage (from repo root):
    python example/make_peaks.py
"""
import json, sys, os, tempfile, urllib.request, urllib.parse
import pyBigWig

API    = "https://www.encodeproject.org"
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(OUTDIR, exist_ok=True)

CHROM, REG_START, REG_END = "chr8", 127_700_000, 128_050_000
MIN_SIGNAL = 2.0   # fold-enrichment threshold — filters weak / noisy peaks

SAMPLES = [
    ("GM12878", "GM12878_peaks.bigBed"),
    ("K562",    "K562_peaks.bigBed"),
]


# ── ENCODE query ──────────────────────────────────────────────────────────────

def query_encode_peaks(biosample):
    params = {
        "type": "File",
        "biosample_ontology.term_name": biosample,
        "assay_title": "Histone ChIP-seq",
        "target.label": "H3K27ac",
        "output_type": "peaks",
        "assembly": "GRCh38",
        "status": "released",
        "limit": 1,
        "format": "json",
    }
    url = API + "/search/?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        data = json.load(r)
    hits = data.get("@graph", [])
    if not hits:
        raise RuntimeError(f"No peak file found for {biosample}")
    href = hits[0]["href"]
    print(f"  Found: {API}{href}")
    return API + href


def download_file(url, dest):
    print(f"  Downloading → {dest}")
    def _prog(b, bs, total):
        done = b * bs
        pct  = min(100, done * 100 // total) if total > 0 else 0
        print(f"\r  {pct:3d}%  {done/1e6:.1f} / {total/1e6:.1f} MB",
              end="", flush=True)
    urllib.request.urlretrieve(url, dest, reporthook=_prog)
    print()


# ── Peak extraction ───────────────────────────────────────────────────────────

def load_peaks_bigbed(path, chrom, start, end, min_signal=MIN_SIGNAL):
    """
    Return list of (start, end) for peaks in the region.
    bigBed extra-string format (narrowPeak):
      name  score  strand  signalValue  pValue  qValue  summit_offset
    """
    bb = pyBigWig.open(path)
    raw = bb.entries(chrom, start, end) or []
    bb.close()
    peaks = []
    for ps, pe, extra in raw:
        fields = extra.split("\t")
        try:
            signal = float(fields[3])   # signalValue
        except (IndexError, ValueError):
            signal = 0.0
        if signal < min_signal:
            continue
        peaks.append((max(ps, start), min(pe, end)))
    return peaks


def overlaps(a_s, a_e, b_list, min_ov=1):
    return any(min(a_e, be) - max(a_s, bs) >= min_ov for bs, be in b_list)


# ── Main ─────────────────────────────────────────────────────────────────────

sample_peaks = {}
for biosample, fname in SAMPLES:
    dest = os.path.join(OUTDIR, fname)
    if not (os.path.exists(dest) and os.path.getsize(dest) > 10_000):
        print(f"\n[fetch] {biosample} H3K27ac peaks (hg38, bigBed)")
        url = query_encode_peaks(biosample)
        download_file(url, dest)
    else:
        print(f"[skip]  {fname} already present  "
              f"({os.path.getsize(dest)//1_000_000} MB)")
    peaks = load_peaks_bigbed(dest, CHROM, REG_START, REG_END)
    sample_peaks[biosample] = peaks
    print(f"  → {len(peaks)} peaks in region (signal ≥ {MIN_SIGNAL})")

gm = sample_peaks["GM12878"]
k5 = sample_peaks["K562"]

out_rows = []
for ps, pe in gm:
    cat  = "shared" if overlaps(ps, pe, k5) else "GM12878_only"
    name = f"GM12878_{(ps+pe)//2//1000}kb"
    out_rows.append((ps, pe, name, cat))

for ps, pe in k5:
    if not overlaps(ps, pe, gm):
        name = f"K562_{(ps+pe)//2//1000}kb"
        out_rows.append((ps, pe, name, "K562_only"))

out_rows.sort(key=lambda r: r[0])

out_path = os.path.join(OUTDIR, "peaks.bed")
with open(out_path, "w") as fh:
    for ps, pe, name, cat in out_rows:
        fh.write(f"{CHROM}\t{ps}\t{pe}\t{name}\t{cat}\n")

gm_only = sum(1 for *_, c in out_rows if c == "GM12878_only")
k5_only = sum(1 for *_, c in out_rows if c == "K562_only")
shared  = sum(1 for *_, c in out_rows if c == "shared")
print(f"\nWrote {len(out_rows)} peaks → {out_path}")
print(f"  GM12878-only: {gm_only}  |  K562-only: {k5_only}  |  shared: {shared}")
