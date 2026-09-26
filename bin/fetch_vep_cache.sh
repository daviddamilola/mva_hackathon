#!/usr/bin/env bash
# Download and extract an Ensembl indexed VEP cache, verified against Ensembl's published CHECKSUMS
# ('sum' algorithm, 1 KiB blocks; macOS and GNU sum agree).
#
# Why chunked and parallel: a single long-lived HTTPS connection to ftp.ensembl.org was throttled
# or reset repeatedly (6 stalls in one run), while fresh connections reached ~28 MB/s. The bytes
# still missing from '<tarball>.part' are fetched as byte-range chunks over fresh connections,
# PARALLEL at a time. Each chunk has stall detection and is retried until its size is exact. The
# chunks are appended in order, and the whole file is checked against Ensembl's checksum before
# extraction. Bytes already on disk are kept. The '.download_ok' marker is written only after
# verification; the pipeline checks for it.
#
# Usage: fetch_vep_cache.sh <species> <assembly> <release> <outdir> [parallel=4] [chunk_mib=512]
set -euo pipefail
species=$1; assembly=$2; release=$3; outdir=$4; parallel=${5:-4}; chunk_mib=${6:-512}
base="https://ftp.ensembl.org/pub/release-${release}/variation/indexed_vep_cache"
tarball="${species}_vep_${release}_${assembly}.tar.gz"
url="${base}/${tarball}"
target="${outdir}/${species}/${release}_${assembly}"
partial="${outdir}/${tarball}.part"
chunkdir="${outdir}/.chunks_${release}"
mkdir -p "$outdir" "$chunkdir"
fsize() { stat -f %z "$1" 2>/dev/null || stat -c %s "$1" 2>/dev/null || echo 0; }

if [[ -f "${target}/.download_ok" ]]; then echo "cache present: ${target}"; exit 0; fi

expected=$(curl -fsSL --retry 5 "${base}/CHECKSUMS" | awk -v f="$tarball" '$3==f {print $1, $2}')
[[ -n "$expected" ]] || { echo "no checksum published for $tarball" >&2; exit 1; }
size=$(curl -fsSIL --retry 5 "$url" | awk 'tolower($1)=="content-length:" {v=$2} END {gsub(/\r/,"",v); print v}')
[[ -n "$size" ]] || { echo "could not read Content-Length" >&2; exit 1; }
touch "$partial"
have=$(fsize "$partial")
echo "have ${have}/${size} bytes; fetching the rest in ${chunk_mib} MiB chunks, ${parallel} in parallel"

# chunk list: index start end
step=$((chunk_mib * 1024 * 1024)); i=0; : > "${chunkdir}/list"
for ((s = have; s < size; s += step)); do
    e=$((s + step - 1)); (( e >= size )) && e=$((size - 1))
    printf '%05d %d %d\n' "$i" "$s" "$e" >> "${chunkdir}/list"; i=$((i + 1))
done

fetch_chunk() {
    local idx=$1 s=$2 e=$3 out="${chunkdir}/c_$1" want=$(( $3 - $2 + 1 ))
    for attempt in $(seq 1 20); do
        [[ -f "$out" && $(fsize "$out") -eq $want ]] && return 0
        curl -fsSL -r "${s}-${e}" --speed-limit 100000 --speed-time 60 --connect-timeout 60 -o "$out" "$url" || sleep 5
    done
    echo "chunk ${idx} failed after 20 attempts" >&2; return 1
}
export -f fetch_chunk fsize; export chunkdir url
xargs -P "$parallel" -L 1 bash -c 'fetch_chunk "$0" "$1" "$2" && echo "chunk $0 ok"' < "${chunkdir}/list"

while read -r idx s e; do
    [[ $(fsize "${chunkdir}/c_${idx}") -eq $((e - s + 1)) ]] || { echo "chunk ${idx} incomplete" >&2; exit 1; }
    cat "${chunkdir}/c_${idx}" >> "$partial" && rm "${chunkdir}/c_${idx}"
done < "${chunkdir}/list"
rm -rf "$chunkdir"

have=$(fsize "$partial")
[[ "$have" -eq "$size" ]] || { echo "incomplete download: ${have}/${size}" >&2; exit 1; }
observed=$(sum "$partial" | awk '{print $1, $2}')
if [[ "$observed" != "$expected" ]]; then
    echo "checksum mismatch for $tarball: expected '$expected' observed '$observed'; deleting" >&2
    rm -f "$partial"; exit 1
fi
echo "checksum ok (${observed}); extracting"
tar -xzf "$partial" -C "$outdir"
rm -f "$partial"
printf 'source=%s\nchecksum_sum=%s\nbytes=%s\ndownloaded_utc=%s\n' "$url" "$observed" "$size" "$(date -u +%FT%TZ)" > "${target}/.download_ok"
echo "cache ok: ${target} (sum ${observed})"
