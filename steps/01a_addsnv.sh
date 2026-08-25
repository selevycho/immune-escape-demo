#!/bin/bash
#SBATCH --job-name=demo_addsnv
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#
# Step 1a — inject substitutions.
#
# Ported from the cohort pipeline's step 2, changed only where the demo
# differs: the mutations come from the truth file shipped with the data
# rather than from a lifted MAF, and paths come from config.sh.
#
# Memory: every parallel bwa process loads the full hg38 index, about
# 6 GB. Two threads inside 16 GB is the combination that held across forty
# samples; more threads were killed partway and bwa then returned an empty
# SAM, which addsnv reports as a site that would not realign.
#
# Runtime scales with mutation count, roughly three seconds per SNV.
#
# Usage:
#   sbatch 01a_addsnv.sh              every sample in config.sh
#   sbatch 01a_addsnv.sh B003 B007    named samples
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
THREADS="${SLURM_CPUS_PER_TASK:-2}"

step_header "STEP 1a — substitutions"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "threads   ${THREADS}"
note "started   $(date '+%F %T') on $(hostname)"

OK=0; SKIP=0; FAIL=0

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    NORMAL="${SAMPLE_DIR}/${SID}_normal.bam"
    TRUTH="${SAMPLE_DIR}/${SID}_truth.tsv"
    TUMOUR="${OUT}/${SID}_tumour_snv.bam"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${NORMAL}" ] || [ ! -s "${TRUTH}" ]; then
        fail "input missing — run setup/fetch_data.sh"
        FAIL=$((FAIL+1)); continue
    fi
    if [ -s "${TUMOUR}" ] && [ -s "${TUMOUR}.bai" ]; then
        note "already injected — skipping"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${OUT}"

    # BAMSurgeon variant format for substitutions:
    #   chrom  pos  pos  VAF  altbase
    awk -F'\t' 'NR>1 && $8=="SNP" {print $2"\t"$3"\t"$3"\t"$7"\t"$5}' \
        "${TRUTH}" > "${OUT}/snvs.txt"

    N_SNV=$(wc -l < "${OUT}/snvs.txt")
    VAF_MED=$(awk '{print $4}' "${OUT}/snvs.txt" | sort -n \
              | awk '{a[NR]=$1} END {if(NR>0) printf "%.3f", a[int(NR/2)+1]; else print "NA"}')
    note "${N_SNV} substitutions, median VAF ${VAF_MED}"

    if [ "${N_SNV}" -lt 1 ]; then
        fail "no substitutions in the truth set"
        FAIL=$((FAIL+1)); continue
    fi

    use_env "${ENV_BAMSURGEON}" || exit 1
    PICARD=$(find "${CONDA_PREFIX}" -name "picard.jar" 2>/dev/null | head -1)
    if [ -z "${PICARD}" ]; then
        fail "picard.jar not found in ${ENV_BAMSURGEON}"
        FAIL=$((FAIL+1)); continue
    fi

    T0=$(date +%s)
    note "$(date '+%T') running addsnv.py ..."

    # addsnv.py writes its scratch directory into the working directory
    cd "${OUT}"
    rm -rf addsnv.tmp

    if ! addsnv.py \
            -v "${OUT}/snvs.txt" \
            -f "${NORMAL}" \
            -r "${REF_FASTA}" \
            -o "${OUT}/tumour_raw.bam" \
            -p "${THREADS}" \
            --aligner mem \
            --picardjar "${PICARD}" \
            --force --insane \
            > "${OUT}/addsnv.log" 2>&1; then
        fail "addsnv.py failed:"
        tail -6 "${OUT}/addsnv.log" | sed 's/^/     /'
        rm -rf "${OUT}/addsnv.tmp" "${OUT}/tumour_raw.bam"*
        FAIL=$((FAIL+1)); continue
    fi

    use_env "${ENV_GATK}" || exit 1
    if ! samtools sort -@ "${THREADS}" -m 1G -o "${TUMOUR}" \
            "${OUT}/tumour_raw.bam" 2>> "${OUT}/addsnv.log"; then
        fail "sorting the injected BAM failed"
        rm -f "${TUMOUR}"
        FAIL=$((FAIL+1)); continue
    fi
    samtools index -@ "${THREADS}" "${TUMOUR}"

    # BAMSurgeon leaves a per-mutation log directory behind
    rm -rf "${OUT}/tumour_raw.bam"* "${OUT}/addsnv.tmp" \
           "${OUT}"/addsnv_logs_*

    T1=$(date +%s)
    printf "  injected in %dm%02ds  |  %s\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
           "$(ls -lh ${TUMOUR} | awk '{print $5}')"
    OK=$((OK+1))
done

step_header "STEP 1a done"
note "injected ${OK}   skipped ${SKIP}   failed ${FAIL}"
[ "${FAIL}" -gt 0 ] && exit 1
exit 0
