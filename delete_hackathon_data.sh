#!/usr/bin/env bash
# delete_hackathon_data.sh
#
# Remove every copy of the "Rare Disease, Real Kid: MVA Hackathon 2026" subject data from this
# checkout, following the organisers' post-hackathon instructions (email of 28 Oct 2026):
#
#   DELETE  the FASTQ, VCF, BAM and CRAM files and any copy, subset, slice or reformat of them,
#           including index files; intermediate files carrying per-variant genotypes (filtered or
#           annotated variant tables), caches and notebook state holding them; stored prompts/logs
#           containing pasted blocks of variant data; model weights trained on the raw data.
#   KEEP    the ranked candidate list (submission), HPO terms, gene/pathway rankings, mechanism
#           analyses, drug candidates, code, report, pitch, leaderboard entry.
#
# Dry run by default: prints what would happen. `--execute` performs it. Every action is written to
# deletion_log_<timestamp>.txt (paths and sizes only). `--scan-claude-logs` additionally reports
# (never deletes) Claude Code transcripts under ~/.claude that contain variant-like blocks.
#
# What it does, in order:
#   0. refuses to run unless launched from the repo root with no nextflow process alive and the
#      hf-mount-nfs mount at local/ unmounted
#   1. preserves software-only material: work/conda -> .conda_cache (nextflow.config now points there)
#   2. moves the retained submission artefacts out of results/ into submission/:
#        results/report/                      report (md, tex, pdf) and its aggregate figures
#        results/pgx/open_targets/*.track2.*  ranked drug candidates (the Track 2 submission list)
#        results/pgx/open_targets/*.pgx_ot.summary.json   aggregate PGx evidence counts
#        results/ot_literature|ot_tractability|ot_neighbourhood   public Open Targets query results
#   3. deletes:
#        data/                 downloaded VCF + index, clinical phenotype docx, HF download cache
#        work/                 every Nextflow task directory (lane CRAMs, slices, pileups, tables)
#        results/              everything else: CRAM/CRAI, panel/annotated VCFs, variant tables,
#                              PGx calls, QC, MultiQC, execution reports, logs, resource copies
#        results_*.zip         the submission archive (bundles panel VCF subsets)
#        .nextflow/ .nextflow.log* .nextflow.pid   session cache and logs
#   4. verifies: no FASTQ/VCF/BAM/CRAM/index left outside refs/ (public resources) and tests/
#      (synthetic), and nothing on the delete list in git history (any branch)
#   5. prints the outside-the-repo checklist that this script cannot do for you
set -euo pipefail

EXECUTE=0; SCAN_CLAUDE=0
for arg in "$@"; do
    case "$arg" in
        --execute) EXECUTE=1 ;;
        --scan-claude-logs) SCAN_CLAUDE=1 ;;
        -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
        *) echo "unknown option: $arg" >&2; exit 2 ;;
    esac
done

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
[ -f main.nf ] && [ -d .git ] || { echo "run from the mva_hackathon repo root" >&2; exit 1; }
if pgrep -f "nextflow run" >/dev/null 2>&1; then echo "a nextflow run is alive; stop it first" >&2; exit 1; fi
if mount | grep -q " ${ROOT}/local "; then echo "local/ is still mounted (hf-mount-nfs); unmount it first" >&2; exit 1; fi

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="${ROOT}/deletion_log_${STAMP}.txt"
MODE=$([ "$EXECUTE" = 1 ] && echo EXECUTE || echo DRY-RUN)
FAILED=()
log() { printf '%s\n' "$*" | tee -a "$LOG"; }
size() { [ -e "$1" ] && du -sh "$1" 2>/dev/null | cut -f1 || echo "-"; }

log "# MVA Hackathon 2026 data deletion ($MODE) $(date -u +%Y-%m-%dT%H:%M:%SZ) in $ROOT"
log "# git: $(git rev-parse --abbrev-ref HEAD) @ $(git rev-parse --short HEAD)"

# ---- 1. preserve software-only conda envs ---------------------------------------------------
log ""; log "## 1. preserve conda envs (software only, no data)"
if [ -d work/conda ]; then
    log "move work/conda ($(size work/conda)) -> .conda_cache"
    if [ "$EXECUTE" = 1 ]; then
        if [ -d .conda_cache ]; then
            for d in work/conda/*; do [ -e ".conda_cache/$(basename "$d")" ] || mv "$d" .conda_cache/; done
        else
            mv work/conda .conda_cache
        fi
    fi
else
    log "work/conda not present"
fi

# ---- 2. retain the submission artefacts -------------------------------------------------------
log ""; log "## 2. retain submission artefacts -> submission/"
keep_move() {   # src dst
    [ -e "$1" ] || { log "  (absent) $1"; return 0; }
    log "  keep  $1 ($(size "$1")) -> $2"
    if [ "$EXECUTE" = 1 ]; then mkdir -p "$(dirname "$2")"; mv "$1" "$2"; fi
}
keep_move results/report submission/report
for f in results/pgx/open_targets/*.track2.* results/pgx/open_targets/*.pgx_ot.summary.json; do
    [ -e "$f" ] && keep_move "$f" "submission/track2_candidates/$(basename "$f")"
done
for d in ot_literature ot_tractability ot_neighbourhood; do keep_move "results/$d" "submission/open_targets/$d"; done
[ -e submission ] || [ "$EXECUTE" = 0 ] || true

# ---- 3. delete ---------------------------------------------------------------------------------
log ""; log "## 3. delete"
TARGETS=(data work results)
for z in ./*.zip; do [ -e "$z" ] && TARGETS+=("${z#./}"); done
TARGETS+=(.nextflow .nextflow.pid)
for l in .nextflow.log .nextflow.log.*; do [ -e "$l" ] && TARGETS+=("$l"); done
for t in "${TARGETS[@]}"; do
    [ -e "$t" ] || continue
    case "$t" in data|work|results)
        # record what was inside, by top level, before it goes
        while IFS= read -r entry; do log "  rm   $t/$entry ($(size "$t/$entry"))"; done < <(ls -A "$t" 2>/dev/null)
        ;;
    *) log "  rm   $t ($(size "$t"))" ;;
    esac
    if [ "$EXECUTE" = 1 ]; then
        # Finder re-creates .DS_Store inside an open folder, which makes the first rm fail on "not empty"
        rm -rf -- "$t" 2>/dev/null || { sleep 2; rm -rf -- "$t"; } || FAILED+=("$t")
    fi
done
if [ "${#FAILED[@]}" -gt 0 ]; then log "FAILED to remove: ${FAILED[*]} (close Finder windows on them and re-run)"; fi

# ---- 4. verify ---------------------------------------------------------------------------------
log ""; log "## 4. verify"
LEFT=$(find . \( -path ./refs -o -path ./tests -o -path ./.git -o -path ./.conda_cache \) -prune -o -type f \
        \( -iname '*.fastq*' -o -iname '*.fq.gz' -o -iname '*.vcf' -o -iname '*.vcf.gz' -o -iname '*.bcf' -o -iname '*.bam' \
           -o -iname '*.cram' -o -iname '*.bai' -o -iname '*.crai' -o -iname '*.tbi' -o -iname '*.csi' -o -iname '*.g.vcf*' \) -print | sort)
if [ -n "$LEFT" ]; then log "sequence/variant files still present (review):"; log "$LEFT"; else log "no FASTQ/VCF/BAM/CRAM/index files remain outside refs/ (public) and tests/ (synthetic)"; fi
GITHIT=$(git log --all --pretty=format: --name-only --diff-filter=A | sort -u \
         | grep -iE '\.(vcf|bcf|bam|bai|cram|crai|fastq|fq|tbi|csi)(\.gz)?$|^results/|^data/|^work/' | grep -v '^tests/stub/synthetic' || true)
if [ -n "$GITHIT" ]; then log "git history contains paths on the delete list (needs git filter-repo):"; log "$GITHIT"; else log "git history (all branches): no FASTQ/VCF/BAM/CRAM/index or data/results/work paths were ever committed"; fi
TRACKED_IDS=$(git ls-files | xargs grep -l "WGS_EX2312012" 2>/dev/null | tr '\n' ' ' || true)
log "tracked files mentioning the sample accession (identifiers only, no genotypes): ${TRACKED_IDS:-none}"
log "retained: $(ls -d submission/* 2>/dev/null | tr '\n' ' ')"
log "disk free: $(df -h "$ROOT" | awk 'NR==2{print $4}')"

# ---- optional: Claude Code transcripts (report only) -------------------------------------------
if [ "$SCAN_CLAUDE" = 1 ]; then
    log ""; log "## Claude Code transcripts with variant-like content (report only; delete by hand if you agree)"
    PROJ="$HOME/.claude/projects/$(echo "$ROOT" | sed 's#/#-#g')"
    if [ -d "$PROJ" ]; then
        for f in "$PROJ"/*.jsonl; do
            [ -e "$f" ] || continue
            n=$(grep -cE 'WGS_EX2312012.*[0-9]+\\t[ACGT]+\\t[ACGT]+\\t|[0-9]+\\t\.\\t[ACGT]+\\t[ACGT]+\\t[0-9.]+\\t(PASS|\.)|\b[01][/|][01]:[0-9]+,[0-9]+' "$f" || true)
            [ "$n" -gt 0 ] && log "  $n variant-like lines: $f ($(size "$f"))"
        done
    else
        log "  no transcript directory at $PROJ"
    fi
fi

# ---- 5. checklist for things outside this checkout --------------------------------------------
log ""; log "## 5. outside this checkout (not done by this script)"
log "  - hf-mount-nfs: make sure the bucket mount at local/ is gone and no other mount point exists"
log "  - ~/.cache/huggingface: hub/xet caches (the pipeline disabled the xet chunk cache; check anyway)"
log "  - Claude Code / notebook transcripts with pasted variant blocks (run with --scan-claude-logs)"
log "  - shell history entries are fine (commands, no variant data); Trash is not involved (rm -rf)"
log "  - Time Machine / iCloud / other backups that captured data/, work/ or results/"
log "  - any copy on another machine or cloud instance; teammates confirm their own"
log "  - after deleting: email cnb@sagebase.org to confirm (deadline 23 Nov 2026)"
log ""
log "$MODE complete. Log: $LOG"
[ "$EXECUTE" = 1 ] || log "(nothing was changed; re-run with --execute to delete)"
