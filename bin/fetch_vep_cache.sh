#!/usr/bin/env bash
# Stream-download and extract an Ensembl indexed VEP cache, verifying Ensembl's published
# CHECKSUMS ('sum' algorithm) on the same stream, so the ~28 GB tarball never touches disk.
# On checksum mismatch the extracted cache is deleted and the script fails.
#
# Usage: fetch_vep_cache.sh <species> <assembly> <release> <outdir>
#   e.g. fetch_vep_cache.sh homo_sapiens GRCh38 116 refs/vep
set -euo pipefail
species=$1; assembly=$2; release=$3; outdir=$4
base="https://ftp.ensembl.org/pub/release-${release}/variation/indexed_vep_cache"
tarball="${species}_vep_${release}_${assembly}.tar.gz"
target="${outdir}/${species}/${release}_${assembly}"
mkdir -p "$outdir"

if [[ -f "${target}/.download_ok" ]]; then
    echo "cache present: ${target}"; exit 0
fi

expected=$(curl -fsSL --retry 5 "${base}/CHECKSUMS" | awk -v f="$tarball" '$3==f {print $1, $2}')
[[ -n "$expected" ]] || { echo "no checksum published for $tarball" >&2; exit 1; }

sumfile=$(mktemp)
curl -fsSL --retry 5 --retry-delay 30 "${base}/${tarball}" \
    | tee >(sum > "$sumfile") \
    | tar -xz -C "$outdir"
sleep 2   # let the process substitution flush
observed=$(awk '{print $1, $2}' "$sumfile")

if [[ "$observed" != "$expected" ]]; then
    echo "checksum mismatch for $tarball: expected '$expected' observed '$observed'" >&2
    rm -rf "$target"
    exit 1
fi
printf 'source=%s/%s\nchecksum_sum=%s\ndownloaded_utc=%s\n' "$base" "$tarball" "$observed" "$(date -u +%FT%TZ)" > "${target}/.download_ok"
echo "cache ok: ${target} (sum ${observed})"
