#!/bin/bash
#SBATCH --job-name=demo_neo
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
#
# Step 5 — neoantigen prediction with mhcflurry.
#
# Ported from the cohort pipeline's step 5.
#
# Takes the mutations that were placed into each sample, translates them
# into mutant peptides, and predicts binding against that sample's own HLA
# genotype as called in step 3.
#
# mhcflurry lives in its own environment because tensorflow does not
# coexist with the rest of the stack.
#
# This step reads the truth file rather than the called variants: the
# question is what the placed mutations would present, not what the caller
# managed to recover. Comparing the two is the report's job.
#
# Usage:
#   sbatch 05_neoantigens.sh                every sample in config.sh
#   sbatch 05_neoantigens.sh B003 B007      named samples
#   LENGTHS=9,10 sbatch 05_neoantigens.sh   9-mers and 10-mers only
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
LENGTHS="${LENGTHS:-8,9,10,11}"
ENGINE="${DEMO_ROOT}/lib/make_neoantigens.py"

step_header "STEP 5 — neoantigen prediction"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "lengths   ${LENGTHS}"
note "started   $(date '+%F %T') on $(hostname)"

if [ ! -s "${ENGINE}" ]; then
    fail "missing ${ENGINE}"
    exit 1
fi
if [ ! -s "${GENCODE_PROT}" ]; then
    fail "GENCODE protein FASTA missing — run setup/fetch_data.sh"
    exit 1
fi

use_env "${ENV_MHCFLURRY}" || exit 1
command -v mhcflurry-predict >/dev/null || {
    fail "mhcflurry-predict not available in ${ENV_MHCFLURRY}"
    exit 1
}

OK=0; SKIP=0; WAIT=0; FAIL=0
T_ALL=$(date +%s)

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    TRUTH="${SAMPLE_DIR}/${SID}_truth.tsv"
    HLA="${OUT}/optitype/${SID}_result.tsv"
    NEO="${OUT}/neoantigens"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${TRUTH}" ]; then
        fail "truth file missing"
        FAIL=$((FAIL+1)); continue
    fi
    if [ ! -s "${HLA}" ]; then
        note "HLA type not ready — step 3 has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${NEO}/neoantigens_per_mutation.tsv" ]; then
        n=$(awk 'NR>1' "${NEO}/neoantigens_per_mutation.tsv" | wc -l)
        s=$(awk -F'\t' 'NR>1 {t+=$6} END {print t+0}' \
            "${NEO}/neoantigens_per_mutation.tsv")
        note "already predicted: ${n} mutations, ${s} strong binders"
        SKIP=$((SKIP+1)); continue
    fi

    T0=$(date +%s)
    note "$(date '+%T') predicting"

    if ! python "${ENGINE}" \
            "${TRUTH}" "${HLA}" "${GENCODE_PROT}" "${NEO}" "${LENGTHS}" \
            > "${OUT}/neoantigens.log" 2>&1; then
        fail "prediction failed:"
        tail -8 "${OUT}/neoantigens.log" | sed 's/^/     /'
        FAIL=$((FAIL+1)); continue
    fi

    T1=$(date +%s)
    grep -E "missense inputs|peptides tested|strong binders|weak binders|at least one binder" \
         "${OUT}/neoantigens.log" | sed 's/^/  /'
    printf "  done in %dm%02ds\n" $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
    OK=$((OK+1))
done

T_END=$(date +%s)
step_header "STEP 5 done"
printf "  elapsed %dm%02ds   predicted %d   skipped %d   waiting %d   failed %d\n" \
    $(( (T_END-T_ALL)/60 )) $(( (T_END-T_ALL)%60 )) \
    "${OK}" "${SKIP}" "${WAIT}" "${FAIL}"

say "binders across the demo"
printf "  %-8s %10s %10s %8s %8s\n" \
       "sample" "mutations" "peptides" "strong" "weak"
T_MUT=0; T_STRONG=0; T_WEAK=0; N=0
for SID in ${SAMPLES}; do
    F="${OUT_DIR}/${SID}/neoantigens/neoantigens_per_mutation.tsv"
    [ -s "${F}" ] || continue
    read -r nm st wk pep < <(awk -F'\t' '
        NR>1 {n++; s+=$6; w+=$7; p+=$4}
        END {print n+0, s+0, w+0, p+0}' "${F}")
    printf "  %-8s %10s %10s %8s %8s\n" "${SID}" "${nm}" "${pep}" "${st}" "${wk}"
    T_MUT=$((T_MUT+nm)); T_STRONG=$((T_STRONG+st)); T_WEAK=$((T_WEAK+wk))
    N=$((N+1))
done

if [ "${N}" -gt 0 ]; then
    echo ""
    note "samples predicted        ${N}"
    note "missense mutations used  ${T_MUT}"
    note "strong binders           ${T_STRONG}"
    note "weak binders             ${T_WEAK}"
    printf "  strong binders per mutation: %.2f\n" \
        "$(echo "${T_STRONG} ${T_MUT}" | awk '{print ($2>0 ? $1/$2 : 0)}')"
    echo ""
    note "The genotype belongs to the 1000 Genomes backbone, not to the"
    note "TCGA donor whose mutations these are. These are the peptides"
    note "those mutations would present in someone carrying this HLA."
fi

[ "${FAIL}" -gt 0 ] && exit 1
exit 0
