#!/bin/bash
#
# Build the conda environments the pipeline needs.
#
# Seven of them, because the tools do not coexist. mhcflurry and pvactools
# want tensorflow; VEP wants a perl stack; OptiType is pinned to an older
# python by razers3; LOHHLA's R dependencies hold samtools at a version
# too old for the rest. Keeping them apart is not tidiness, it is the only
# arrangement conda can resolve.
#
# The yml files pin exact versions — these are the ones the published
# numbers came from. If conda cannot resolve one, channels do drop old
# builds over time; see the note at the end.
#
# Environments already present are left alone. Rebuilding one takes as
# long as building it, and none of them changes.
#
# Usage:
#   ./install.sh                every environment
#   ./install.sh mhc_env        one of them
#   ./install.sh --list         say what would be built and stop
#   ./install.sh --force        rebuild even if present
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

YML_DIR="${HERE}/environments"
FORCE=0; LIST=0
TARGETS=""
for a in "$@"; do
    case "$a" in
        --force) FORCE=1 ;;
        --list)  LIST=1 ;;
        --*)     echo "unknown option: $a" >&2; exit 1 ;;
        *)       TARGETS="${TARGETS} $a" ;;
    esac
done

WHAT_FOR=(
  "${ENV_BAMSURGEON}:placing mutations"
  "${ENV_GATK}:variant calling"
  "${ENV_OPTITYPE}:HLA typing"
  "${ENV_MHCFLURRY}:binding prediction"
  "${ENV_LOHHLA}:HLA loss"
  "${ENV_VEP}:annotation"
  "${ENV_PVACSEQ}:neoantigens"
)

if [ -z "${TARGETS}" ]; then
    for pair in "${WHAT_FOR[@]}"; do
        TARGETS="${TARGETS} ${pair%%:*}"
    done
fi

purpose_of() {
    for pair in "${WHAT_FOR[@]}"; do
        [ "${pair%%:*}" = "$1" ] && { echo "${pair#*:}"; return; }
    done
    echo "unknown"
}

step_header "environments"

if [ -z "${CONDA_SH:-}" ] || [ ! -f "${CONDA_SH}" ]; then
    fail "conda not found"
    note "install miniconda from https://docs.conda.io/en/latest/miniconda.html"
    note "or set CONDA_SH in config.sh to its profile.d/conda.sh"
    exit 1
fi
# shellcheck disable=SC1090
source "${CONDA_SH}"
note "conda at $(dirname "$(dirname "${CONDA_SH}")")"

# mamba resolves these far faster than conda does; twenty minutes against
# an hour is the usual difference on a stack this size
SOLVER="conda"
if command -v mamba >/dev/null 2>&1; then
    SOLVER="mamba"
fi
note "solver ${SOLVER}"

if [ "${LIST}" = "1" ]; then
    say "would build"
    for e in ${TARGETS}; do
        yml="${YML_DIR}/${e}.yml"
        state="missing yml"
        [ -s "${yml}" ] && state="$(wc -l < "${yml}") packages"
        conda env list | awk '{print $1}' | grep -qx "${e}" && state="already present"
        printf "  %-18s %-22s %s\n" "${e}" "$(purpose_of "${e}")" "${state}"
    done
    exit 0
fi

BUILT=0; SKIPPED=0; FAILED=0
declare -a BROKEN

for e in ${TARGETS}; do
    yml="${YML_DIR}/${e}.yml"

    echo ""
    echo "----------------------------------------------------------------------"
    printf " %s — %s\n" "${e}" "$(purpose_of "${e}")"

    if [ ! -s "${yml}" ]; then
        fail "no ${yml}"
        FAILED=$((FAILED+1)); BROKEN+=("${e}")
        continue
    fi

    if conda env list | awk '{print $1}' | grep -qx "${e}"; then
        if [ "${FORCE}" = "0" ]; then
            note "already present — skipping"
            SKIPPED=$((SKIPPED+1)); continue
        fi
        note "removing the existing environment"
        conda env remove -n "${e}" -y > /dev/null 2>&1
    fi

    note "$(date '+%T') building from $(basename "${yml}")"
    T0=$(date +%s)

    if ${SOLVER} env create -f "${yml}" -n "${e}" \
            > "/tmp/install_${e}.log" 2>&1; then
        T1=$(date +%s)
        n=$(conda list -n "${e}" 2>/dev/null | grep -vc "^#")
        printf "  built in %dm%02ds, %s packages\n" \
            $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${n}"
        BUILT=$((BUILT+1))
    else
        fail "failed — see /tmp/install_${e}.log"
        tail -8 "/tmp/install_${e}.log" | sed 's/^/     /'
        FAILED=$((FAILED+1)); BROKEN+=("${e}")
    fi
done

step_header "done"
note "built ${BUILT}   already present ${SKIPPED}   failed ${FAILED}"

if [ "${FAILED}" -gt 0 ]; then
    echo ""
    say "if an environment would not resolve"
    echo ""
    note "The yml files pin exact versions, and conda channels remove old"
    note "builds after a year or two. When that happens the fix is to relax"
    note "the pins:"
    note ""
    note "  sed -E 's/=[0-9][^=]*$//' setup/environments/NAME.yml > relaxed.yml"
    note "  conda env create -f relaxed.yml -n NAME"
    note ""
    note "That gives current versions of the same packages rather than the"
    note "exact stack these results came from. The pipeline will run; the"
    note "numbers may shift slightly."
    note ""
    note "unresolved: ${BROKEN[*]}"
    exit 1
fi

echo ""
note "next: ./setup/get_lohhla.sh, then ./setup/fetch_data.sh"
