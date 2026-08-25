#!/bin/bash
#SBATCH --job-name=demo_mutect2
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#
# Step 2 — somatic variant calling.
#
# Ported from the cohort pipeline's step 3. Two things carried over
# unchanged because both were learned the hard way.
#
# Read groups are rewritten first. Tumour and normal are both sliced from
# the same 1000 Genomes individual and carry an identical SM tag; Mutect2
# refuses a pair it cannot distinguish and the error does not say why.
# samtools addreplacerg rewrites the tag as a stream without realigning.
#
# The unfiltered VCF is kept. FilterMutectCalls only annotates, so the
# unfiltered file holds every variant the caller emitted — and "never
# called" and "called then rejected" are different failures. The report
# needs to tell them apart.
#
# One run covers both kinds of mutation, because step 1b leaves both in
# the same BAM. That is also how a real sample arrives.
#
# Usage:
#   sbatch 02_mutect2.sh              every sample in config.sh
#   sbatch 02_mutect2.sh B003 B007    named samples
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
THREADS="${SLURM_CPUS_PER_TASK:-2}"

step_header "STEP 2 — Mutect2"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "threads   ${THREADS}"
note "started   $(date '+%F %T') on $(hostname)"

use_env "${ENV_GATK}" || exit 1

CALLED=0; SKIP=0; WAIT=0; FAIL=0
T_ALL=$(date +%s)

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    NORMAL="${SAMPLE_DIR}/${SID}_normal.bam"
    TUMOUR="${OUT}/${SID}_tumour.bam"
    MDIR="${OUT}/mutect2"
    RAW="${MDIR}/${SID}.unfiltered.vcf.gz"
    VCF="${MDIR}/${SID}.filtered.vcf.gz"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${TUMOUR}" ]; then
        note "tumour BAM not ready — step 1b has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi
    if [ ! -s "${NORMAL}" ]; then
        fail "normal BAM missing"
        FAIL=$((FAIL+1)); continue
    fi
    if [ -s "${VCF}" ]; then
        note "already called — skipping"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${MDIR}"
    T0=$(date +%s)

    # both BAMs come from the same individual, so SM must be rewritten
    note "$(date '+%T') rewriting read groups"
    if ! samtools addreplacerg -@ "${THREADS}" \
            -r "ID:TUMOR\tSM:${SID}_TUMOR\tPL:ILLUMINA\tLB:panel\tPU:sim" \
            -o "${OUT}/t.rg.bam" "${TUMOUR}" 2>/dev/null \
       || ! samtools index -@ "${THREADS}" "${OUT}/t.rg.bam"; then
        fail "could not prepare the tumour BAM"
        rm -f "${OUT}"/t.rg.bam*
        FAIL=$((FAIL+1)); continue
    fi

    if ! samtools addreplacerg -@ "${THREADS}" \
            -r "ID:NORMAL\tSM:${SID}_NORMAL\tPL:ILLUMINA\tLB:panel\tPU:sim" \
            -o "${OUT}/n.rg.bam" "${NORMAL}" 2>/dev/null \
       || ! samtools index -@ "${THREADS}" "${OUT}/n.rg.bam"; then
        fail "could not prepare the normal BAM"
        rm -f "${OUT}"/[tn].rg.bam*
        FAIL=$((FAIL+1)); continue
    fi

    note "$(date '+%T') Mutect2 over the panel"
    if ! gatk --java-options "-Xmx${MEM_GB}g" Mutect2 \
            -R "${REF_FASTA}" \
            -I "${OUT}/t.rg.bam" \
            -I "${OUT}/n.rg.bam" \
            -normal "${SID}_NORMAL" \
            -L "${PANEL_BED}" \
            --native-pair-hmm-threads "${THREADS}" \
            -O "${RAW}" \
            > "${MDIR}/mutect2.log" 2>&1; then
        fail "Mutect2 failed:"
        grep -iE "error|exception" "${MDIR}/mutect2.log" | tail -4 \
            | sed 's/^/     /'
        rm -f "${OUT}"/[tn].rg.bam*
        FAIL=$((FAIL+1)); continue
    fi

    note "$(date '+%T') FilterMutectCalls"
    if ! gatk --java-options "-Xmx$(( MEM_GB / 2 ))g" FilterMutectCalls \
            -R "${REF_FASTA}" \
            -V "${RAW}" \
            -O "${VCF}" \
            > "${MDIR}/filter.log" 2>&1; then
        fail "FilterMutectCalls failed:"
        grep -iE "error|exception" "${MDIR}/filter.log" | tail -4 \
            | sed 's/^/     /'
        rm -f "${OUT}"/[tn].rg.bam*
        FAIL=$((FAIL+1)); continue
    fi

    # the rewritten copies are only needed while the caller runs
    rm -f "${OUT}"/[tn].rg.bam*

    T1=$(date +%s)

    N_ALL=$(zcat "${VCF}" | grep -vc "^#")
    N_PASS=$(zcat "${VCF}" | grep -v "^#" | awk -F'\t' '$7=="PASS"' | wc -l)
    N_IND=$(zcat "${VCF}" | grep -v "^#" \
            | awk -F'\t' '$7=="PASS" && length($4)!=length($5)' | wc -l)
    N_TRUTH=$(awk -F'\t' 'NR>1' "${SAMPLE_DIR}/${SID}_truth.tsv" | wc -l)

    printf "  %dm%02ds  |  %s records, %s PASS (%s indels) against %s placed\n" \
        $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
        "${N_ALL}" "${N_PASS}" "${N_IND}" "${N_TRUTH}"
    CALLED=$((CALLED+1))
done

T_END=$(date +%s)
step_header "STEP 2 done"
printf "  elapsed %dh%02dm   called %d   skipped %d   waiting %d   failed %d\n" \
    $(( (T_END-T_ALL)/3600 )) $(( ((T_END-T_ALL)%3600)/60 )) \
    "${CALLED}" "${SKIP}" "${WAIT}" "${FAIL}"

n_want=$(echo ${SAMPLES} | wc -w)
n_have=$(ls "${OUT_DIR}"/*/mutect2/*.filtered.vcf.gz 2>/dev/null | wc -l)
note "VCFs: ${n_have} / ${n_want}"

[ "${FAIL}" -gt 0 ] && exit 1
exit 0
