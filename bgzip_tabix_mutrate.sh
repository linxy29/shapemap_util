#!/usr/bin/env bash
# Re-compress every *mutrate.txt.gz in a folder with bgzip and build a
# tabix index so single-gene rows can be fetched in O(log N).
#
# Usage:
#   bgzip_tabix_mutrate.sh <folder> [parallel_jobs]
#
# Output (placed next to each input):
#   <stem>.mutrate.txt.bgz       bgzipped, sorted by gene then start
#   <stem>.mutrate.txt.bgz.tbi   tabix index (col 1 = gene, col 2/3 = start/end)
#
# Lookup example:
#   tabix <stem>.mutrate.txt.bgz Camk2a

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "usage: $0 <folder> [parallel_jobs]" >&2
    exit 1
fi

DIR=$1
JOBS=${2:-1}

if [[ ! -d "$DIR" ]]; then
    echo "error: not a directory: $DIR" >&2
    exit 1
fi

command -v bgzip >/dev/null || { echo "error: bgzip not on PATH" >&2; exit 1; }
command -v tabix >/dev/null || { echo "error: tabix not on PATH" >&2; exit 1; }

shopt -s nullglob
files=( "$DIR"/*mutrate.txt.gz )
shopt -u nullglob

if [[ ${#files[@]} -eq 0 ]]; then
    echo "no *mutrate.txt.gz files in $DIR" >&2
    exit 0
fi

process_one() {
    local f=$1
    local out="${f%.gz}.bgz"

    if [[ -s "${out}.tbi" ]]; then
        echo "skip (already indexed): $f"
        return 0
    fi

    echo "indexing $f"
    # Sort by gene (col 1) then start (col 2, numeric) before bgzipping —
    # tabix requires this order even if the input usually already is.
    zcat "$f" \
        | LC_ALL=C sort -k1,1 -k2,2n -S 2G \
        | bgzip -c > "${out}.tmp"
    mv "${out}.tmp" "$out"
    tabix -f -s 1 -b 2 -e 3 "$out"
}

export -f process_one

if [[ "$JOBS" -gt 1 ]]; then
    printf '%s\n' "${files[@]}" \
        | xargs -n 1 -P "$JOBS" -I{} bash -c 'process_one "$@"' _ {}
else
    for f in "${files[@]}"; do
        process_one "$f"
    done
fi

echo "done. Test with: tabix <file>.mutrate.txt.bgz <gene>"
