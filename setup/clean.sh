#!/bin/bash
#
# Remove what a previous run left behind.
#
# Three kinds of leftover, and they are not equally safe to delete.
#
# Scratch is always disposable: BAMSurgeon's tmp directories, its
# per-mutation log folders, the unsorted intermediates, the working copies
# made for addindel. None of it is read by anything.
#
# Results are the output of the pipeline and deleting them means running
# it again, which is sometimes exactly what is wanted and sometimes not.
#
# Inputs and reference data are never touched. They took an hour to fetch
# and are not produced by any step.
#
# Nothing is deleted without being listed first unless --force is given.
#
# Usage:
#   ./clean.sh                 scratch only, the safe default
#   ./clean.sh --results       scratch and results, keeping inputs
#   ./clean.sh --results B003  one sample's results
#   ./clean.sh --logs          the SLURM logs in steps/
#   ./clean.sh --all           everything the pipeline produced
#   add --force to skip the confirmation
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

DO_RESULTS=0; DO_LOGS=0; FORCE=0
TARGETS=""
for a in "$@"; do
    case "$a" in
        --results) DO_RESULTS=1 ;;
        --logs)    DO_LOGS=1 ;;
        --all)     DO_RESULTS=1; DO_LOGS=1 ;;
        --force)   FORCE=1 ;;
        --*)       echo "unknown option: $a" >&2; exit 1 ;;
        *)         TARGETS="${TARGETS} $a" ;;
    esac
done
TARGETS="${TARGETS:-${SAMPLES}}"

declare -a DOOMED
add() { [ -e "$1" ] && DOOMED+=("$1"); }

say "scratch"
for SID in ${TARGETS}; do
    d="${OUT_DIR}/${SID}"
    [ -d "${d}" ] || continue
    for p in "${d}"/addsnv.tmp "${d}"/addindel.tmp \
             "${d}"/addsnv_logs_* "${d}"/addindel_logs_* \
             "${d}"/tumour_raw.bam* "${d}"/raw.bam* \
             "${d}"/snv_raw.bam* "${d}"/indel_raw.bam* \
             "${d}"/with_snv.bam* "${d}"/with_indel.bam* \
             "${d}"/"${SID}"_input.bam* \
             "${d}"/t.rg.bam* "${d}"/n.rg.bam* \
             "${d}"/rest.bam "${d}"/merged.bam \
             "${d}"/locus.bam "${d}"/locus_thin.bam \
             "${d}"/target_thin.bam; do
        add "${p}"
    done
done

if [ "${DO_RESULTS}" = "1" ]; then
    say "results"
    for SID in ${TARGETS}; do
        add "${OUT_DIR}/${SID}"
    done
    add "${OUT_DIR}/loh_simulation.tsv"
    add "${OUT_DIR}/loh_design.tsv"
    add "${OUT_DIR}/summary.tsv"
fi

if [ "${DO_LOGS}" = "1" ]; then
    say "logs"
    for f in "${DEMO_ROOT}"/steps/*.log; do
        add "${f}"
    done
fi

if [ "${#DOOMED[@]}" -eq 0 ]; then
    note "nothing to remove"
    exit 0
fi

total=0
say "would remove"
for p in "${DOOMED[@]}"; do
    sz=$(du -sh "${p}" 2>/dev/null | cut -f1)
    b=$(du -sb "${p}" 2>/dev/null | cut -f1)
    total=$((total + ${b:-0}))
    printf '  %8s  %s\n' "${sz}" "${p#${DATA_DIR}/}"
done
note ""
note "$(echo ${total} | awk '{printf "%.1f GB in %d items", $1/1073741824, '"${#DOOMED[@]}"'}')"

if [ "${FORCE}" != "1" ]; then
    echo ""
    read -r -p "  remove these? [y/N] " ans
    case "${ans}" in
        y|Y|yes) ;;
        *) note "nothing removed"; exit 0 ;;
    esac
fi

for p in "${DOOMED[@]}"; do
    rm -rf "${p}"
done

say "done"
note "removed ${#DOOMED[@]} items"
note ""
note "inputs and reference data were not touched:"
note "  ${SAMPLE_DIR}"
note "  ${REF_DIR}"
