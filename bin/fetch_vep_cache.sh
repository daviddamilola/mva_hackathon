#!/usr/bin/env bash
# Download and extract an Ensembl indexed VEP cache, verified against Ensembl's published CHECKSUMS
# ('sum' algorithm, 1 KiB blocks; macOS and GNU sum agree).
#
# Resumable: the tarball is fetched with 'curl -C -' and stall detection (--speed-limit/-time),
# retried until its size equals the server's Content-Length. A dropped or stalled connection then
# costs a retry, not the whole 28 GB (a first streaming version stalled at 13 GB with no timeout).
# After verification the tarball is extracted and deleted. Peak disk ~2x the tarball size.
# The '.download_ok' marker is written only after verification; the pipeline checks for it.
#
# Usage: fetch_vep_cache.sh <species> <assembly> <release> <outdir>
#   e.g. fetch_vep_cache.sh homo_sapiens GRCh38 116 refs/vep
set -euo pipefail
species=$1; assembly=$2; release=$3; outdir=$4
base="https://ftp.ensembl.org/pub/release-${release}/variation/indexed_vep_cache"
tarball="${species}_vep_${release}_${assembly}.tar.gz"
target="${outdir}/${species}/${release}_${assembly}"
partial="${outdir}/${tarball}.part"
mkdir -p "$outdir"

if [[ -f "${target}/.download_ok" ]]; then
    echo "cache present: ${target}"; exit 0
fi

expected=$(curl -fsSL --retry 5 "${base}/CHECKSUMS" | awk -v f="$tarball" '$3==f {print $1, $2}')
[[ -n "$expected" ]] || { echo "no checksum published for $tarball" >&2; exit 1; }
size=$(curl -fsSIL --retry 5 "${base}/${tarball}" | awk 'tolower($1)=="content-length:" {v=$2} END {gsub(/\r/,"",v); print v}')
[[ -n "$size" ]] || { echo "could not read Content-Length" >&2; exit 1; }

for attempt in $(seq 1 50); do
    have=$(stat -f %z "$partial" 2>/dev/null || stat -c %s "$partial" 2>/dev/null || echo 0)
    [[ "$have" -ge "$size" ]] && break
    echo "attempt ${attempt}: ${have}/${size} bytes"
    curl -fsSL -C - --speed-limit 50000 --speed-time 120 --connect-timeout 60 \
        -o "$partial" "${base}/${tarball}" || sleep 20
done
have=$(stat -f %z "$partial" 2>/dev/null || stat -c %s "$partial")
[[ "$have" -eq "$size" ]] || { echo "incomplete download: ${have}/${size}" >&2; exit 1; }

observed=$(sum "$partial" | awk '{print $1, $2}')
if [[ "$observed" != "$expected" ]]; then
    echo "checksum mismatch for $tarball: expected '$expected' observed '$observed'; deleting" >&2
    rm -f "$partial"
    exit 1
fi

tar -xzf "$partial" -C "$outdir"
rm -f "$partial"
printf 'source=%s/%s\nchecksum_sum=%s\nbytes=%s\ndownloaded_utc=%s\n' "$base" "$tarball" "$observed" "$size" "$(date -u +%FT%TZ)" > "${target}/.download_ok"
echo "cache ok: ${target} (sum ${observed})"
