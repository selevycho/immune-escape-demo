#!/bin/bash
#
# Run the whole demo.
#
# Eight stages, in order, on five samples. Each stage skips what is
# already done, so an interrupted run can simply be started again and will
# pick up where it stopped.
#
# On a cluster the stages are submitted to SLURM with dependencies, so
# nothing starts before its input exists and the whole thing runs
# unattended. Without SLURM they run in sequence in this shell, which
# works but ties up the terminal for about three hours.
#
# Expect roughly:
#   substitutions      5 min per sample
#   indels             2 min
#   variant calling   20 min
#   HLA typing         3 min
#   loss simulation    1 min
#   LOHHLA             1 min
#   neoantigens        2 min
#   pVACseq           12 min
#
# About 45 minutes of compute per sample; in parallel on a cluster, under
# an hour in total.
#
# Usage:
#   ./run_demo.sh                 everything
#   ./run_demo.sh --local         sequentially in this shell, no SLURM
#   ./run_demo.sh --from 3        start at stage 3
#   ./run_demo.sh --only 2        one stage
#   ./run_demo.sh --dry-run       say what would happen
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/config.sh"

LOCAL=0; DRY=0; FROM=1; ONLY=""
while [ $# -gt 0 ]; do
    case "$1" in
        --local)   LOCAL=1 ;;
        --dry-run) DRY=1 ;;
        --from)    FROM="$2"; shift ;;
        --only)    ONLY="$2"; shift ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
    shift
done

# stage number : script : what it does
STAGES=(
  "1:01a_addsnv.sh:place substitutions"
  "2:01b_addindel.sh:place indels"
  "3:02_mutect2.sh:call variants"
  "4:03_optitype.sh:type HLA"
  "5:04_loh.sh:simulate HLA loss"
  "6:06_lohhla.sh:test for HLA loss"
  "7:05_neoantigens.sh:predict binding with mhcflurry"
  "8:07_pvacseq.sh:predict binding with NetMHCpan"
)

step_header "IMMUNE ESCAPE CONVERGENCE — demo"
note "samples   $(echo ${SAMPLES} | tr '\n' ' ')"
note "data      ${DATA_DIR}"
note "results   ${OUT_DIR}"
note "mode      $([ "${LOCAL}" = "1" ] && echo "sequential" || echo "SLURM")"

# ---------------------------------------------------------------- checks
# a dry run is asked what would happen, not whether it could happen now,
# so the dependency check is reported rather than enforced
if [ "${DRY}" = "1" ]; then
    say "dependencies"
    "${DEMO_ROOT}/setup/check_deps.sh" --quiet 2>&1 | sed 's/^/  /' || true
fi

if [ "${DRY}" = "0" ] && \
   ! "${DEMO_ROOT}/setup/check_deps.sh" --quiet > /tmp/demo_check.$$ 2>&1; then
    cat /tmp/demo_check.$$
    rm -f /tmp/demo_check.$$
    echo ""
    fail "something required is missing — see above"
    note "run setup/check_deps.sh for the full picture"
    exit 1
fi
rm -f /tmp/demo_check.$$
[ "${DRY}" = "0" ] && note "dependencies ok"

if [ "${DRY}" = "1" ]; then
    say "would run"
    for entry in "${STAGES[@]}"; do
        n="${entry%%:*}"
        rest="${entry#*:}"
        script="${rest%%:*}"
        what="${rest#*:}"
        [ -n "${ONLY}" ] && [ "${ONLY}" != "${n}" ] && continue
        [ "${n}" -lt "${FROM}" ] && continue
        printf "  %s  %-22s %s\n" "${n}" "${script}" "${what}"
    done
    exit 0
fi

mkdir -p "${OUT_DIR}"
T_START=$(date +%s)

# ------------------------------------------------------------ submission
if [ "${LOCAL}" = "0" ] && command -v sbatch >/dev/null; then
    say "submitting"
    PREV_IDS=""
    for entry in "${STAGES[@]}"; do
        n="${entry%%:*}"
        rest="${entry#*:}"
        script="${rest%%:*}"
        what="${rest#*:}"
        [ -n "${ONLY}" ] && [ "${ONLY}" != "${n}" ] && continue
        [ "${n}" -lt "${FROM}" ] && continue

        # one job per sample, so a slow sample does not hold up the rest;
        # the next stage waits for all of them
        IDS=""
        for SID in ${SAMPLES}; do
            DEP=""
            [ -n "${PREV_IDS}" ] && DEP="--dependency=afterok:${PREV_IDS}"
            id=$(sbatch --parsable ${DEP} \
                 --export=ALL,DATA_DIR="${DATA_DIR}" \
                 "${DEMO_ROOT}/steps/${script}" "${SID}" 2>/dev/null)
            [ -n "${id}" ] && IDS="${IDS}${IDS:+:}${id}"
        done

        if [ -z "${IDS}" ]; then
            fail "stage ${n} (${script}) could not be submitted"
            exit 1
        fi
        printf "  %s  %-22s %s\n" "${n}" "${what}" \
               "$(echo ${IDS} | tr ':' ' ' | wc -w) jobs"
        PREV_IDS="${IDS}"
    done

    echo ""
    note "submitted. The stages are chained: each waits for the one"
    note "before it, so nothing runs against a file that is not there."
    note ""
    note "watch progress with"
    note "  squeue -u \$USER"
    note ""
    note "when the queue empties"
    note "  ${DEMO_ROOT}/steps/08_report.sh"
    exit 0
fi

# ------------------------------------------------------------ sequential
say "running here"
note "no SLURM, or --local given; this will take a few hours"

for entry in "${STAGES[@]}"; do
    n="${entry%%:*}"
    rest="${entry#*:}"
    script="${rest%%:*}"
    what="${rest#*:}"
    [ -n "${ONLY}" ] && [ "${ONLY}" != "${n}" ] && continue
    [ "${n}" -lt "${FROM}" ] && continue

    echo ""
    echo "######################################################################"
    printf "# stage %s of %s — %s\n" "${n}" "${#STAGES[@]}" "${what}"
    echo "######################################################################"

    T0=$(date +%s)
    if ! bash "${DEMO_ROOT}/steps/${script}"; then
        T1=$(date +%s)
        echo ""
        fail "stage ${n} (${script}) failed after $(( (T1-T0)/60 ))m"
        note ""
        note "Fix what it reports and start again from here:"
        note "  ./run_demo.sh --from ${n}"
        note ""
        note "Stages skip what is already done, so nothing is repeated."
        exit 1
    fi
    T1=$(date +%s)
    note "stage ${n} took $(( (T1-T0)/60 ))m$(( (T1-T0)%60 ))s"
done

# ---------------------------------------------------------------- report
echo ""
bash "${DEMO_ROOT}/steps/08_report.sh"

T_END=$(date +%s)
step_header "DONE"
printf "  total %dh%02dm\n" \
    $(( (T_END-T_START)/3600 )) $(( ((T_END-T_START)%3600)/60 ))
note "results in ${OUT_DIR}"
