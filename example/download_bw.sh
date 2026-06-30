#!/usr/bin/env bash
# Download H3K27ac fold-change BigWig files from the ENCODE portal.
# Queries the ENCODE REST API so results always reflect the latest reprocessed files.
#
# Requirements: curl, wget (or change to curl -O), python3
# Usage: bash example/download_bw.sh

set -euo pipefail

API="https://www.encodeproject.org"
OUTDIR="$(dirname "$(realpath "$0")")/data"
mkdir -p "$OUTDIR"

query_and_download() {
    local biosample="$1"
    local outfile="$2"
    local encoded_biosample="${biosample// /+}"

    echo "──────────────────────────────────────────────"
    echo "Querying ENCODE for: ${biosample} H3K27ac BigWig (hg38, fold-change over control)"

    local json href
    json=$(curl -sL \
        "${API}/search/?type=File\
&biosample_ontology.term_name=${encoded_biosample}\
&assay_title=Histone+ChIP-seq\
&target.label=H3K27ac\
&file_format=bigWig\
&output_type=fold+change+over+control\
&assembly=GRCh38\
&status=released\
&replication_type=isogenic\
&limit=1\
&format=json" \
        -H "Accept: application/json")

    href=$(python3 - <<'PYEOF'
import json, sys
d = json.loads(sys.stdin.read())
if not d.get("@graph"):
    raise SystemExit("No files found — check biosample name or ENCODE portal availability")
print(d["@graph"][0]["href"])
PYEOF
        <<< "$json")

    echo "  File: ${API}${href}"
    wget -q --show-progress -c -O "${OUTDIR}/${outfile}" "${API}${href}"
    echo "  Saved: ${OUTDIR}/${outfile}"
}

query_and_download "GM12878" "GM12878_H3K27ac_fc.bw"
query_and_download "K562"    "K562_H3K27ac_fc.bw"

echo ""
echo "══════════════════════════════════════════════"
echo "Download complete. Generate the example figure:"
echo ""
echo "  python locus_plot.py \\"
echo "      --region chr8:127700000-128050000 \\"
echo "      --config example/tracks.ini \\"
echo "      --out example/output.pdf"
echo ""
