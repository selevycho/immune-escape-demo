#!/bin/bash
#SBATCH --job-name=demo_loh
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.log
#
# Step 4 — simulate loss of heterozygosity at one HLA locus.
#
# Ported from the cohort pipeline's step 6a.
#
# What is simulated, and what is not. Thinning is locus-wide: a random
# fraction of the reads across the whole window is dropped. It reduces
# total coverage there without creating allelic imbalance, and LOHHLA
# tests for allelic imbalance. An allele-specific version was attempted
# and abandoned — on panel data at this depth the two alleles cannot be
# separated cleanly, and reads split 152 against 20 where an even split
# was needed. So what follows measures sensitivity to coverage loss, and
# the report says so rather than claiming otherwise.
#
# The locus to hit is the first heterozygous one in the order A, B, C. A
# homozygous locus has no second allele to lose. The untouched loci are
# within-sample controls: a method that reports loss there is reporting
# something that was not done.
#
# Both BAMs are copied into a private directory per sample. LOHHLA writes
# intermediates beside its input and derives working directory names from
# the input filename, so two samples cannot share a directory. The normal
# copy is never modified — it is the reference the tumour is compared
# against, and thinning it would invert the question.
#
# Usage:
#   sbatch 04_loh.sh                    every sample in config.sh
#   sbatch 04_loh.sh B003 B007          named samples
#   LOH_KEEP=0.5 sbatch 04_loh.sh       keep half instead of 35%
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"
KEEP="${LOH_KEEP:-0.35}"

# hg38 windows, the same ones used for typing
A_S=29932000; A_E=29956000
B_S=31343000; B_E=31367000
C_S=31258000; C_E=31282000

step_header "STEP 4 — simulating HLA loss"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "keeping   ${KEEP} of the target locus, other loci untouched"
note "started   $(date '+%F %T') on $(hostname)"

use_env "${ENV_GATK}" || exit 1

DESIGN="${OUT_DIR}/loh_design.tsv"
[ -s "${DESIGN}" ] || printf \
  "sample\ttarget_locus\ttarget_alleles\tcontrols\tkeep\tA_before\tA_after\tB_before\tB_after\tC_before\tC_after\n" \
  > "${DESIGN}"

OK=0; SKIP=0; HOMO=0; WAIT=0; FAIL=0

depth_of() {   # $1 bam  $2 start  $3 end
    samtools depth -a -r "chr6:${2}-${3}" "$1" 2>/dev/null \
      | awk '{s+=$3;n++} END {printf "%.1f", (n?s/n:0)}'
}

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    NORMAL="${SAMPLE_DIR}/${SID}_normal.bam"
    TUMOUR="${OUT}/${SID}_tumour.bam"
    HLA="${OUT}/optitype/${SID}_result.tsv"
    LDIR="${OUT}/loh"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${TUMOUR}" ] || [ ! -s "${HLA}" ]; then
        note "tumour BAM or HLA type not ready"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${LDIR}/${SID}_tumour_LOH.bam.bai" ]; then
        note "already built — skipping"
        SKIP=$((SKIP+1)); continue
    fi

    read -r A1 A2 B1 B2 C1 C2 < <(awk 'NR==2 {print $2,$3,$4,$5,$6,$7}' "${HLA}")
    printf "   A: %-10s %-10s  B: %-10s %-10s  C: %-10s %-10s\n" \
           "${A1}" "${A2}" "${B1}" "${B2}" "${C1}" "${C2}"

    # first heterozygous locus wins, in the order A, B, C
    TARGET=""; T_S=""; T_E=""; T_ALLELES=""; CONTROLS=""
    if   [ "${A1}" != "${A2}" ] && [ "${A1#*\*}" != "${A1}" ] && [ "${A2#*\*}" != "${A2}" ]; then
        TARGET=A; T_S=${A_S}; T_E=${A_E}; T_ALLELES="${A1}/${A2}"; CONTROLS="B,C"
    elif [ "${B1}" != "${B2}" ] && [ "${B1#*\*}" != "${B1}" ] && [ "${B2#*\*}" != "${B2}" ]; then
        TARGET=B; T_S=${B_S}; T_E=${B_E}; T_ALLELES="${B1}/${B2}"; CONTROLS="A,C"
    elif [ "${C1}" != "${C2}" ] && [ "${C1#*\*}" != "${C1}" ] && [ "${C2#*\*}" != "${C2}" ]; then
        TARGET=C; T_S=${C_S}; T_E=${C_E}; T_ALLELES="${C1}/${C2}"; CONTROLS="A,B"
    else
        note "homozygous at every typed locus — nothing to lose"
        printf "%s\thomozygous_all\n" "${SID}" >> "${OUT_DIR}/loh_skipped.log"
        HOMO=$((HOMO+1)); continue
    fi

    note "target HLA-${TARGET} (${T_ALLELES})   controls ${CONTROLS}"
    mkdir -p "${LDIR}"
    T0=$(date +%s)

    # ---- private copies; nothing upstream is modified ----
    if [ ! -s "${LDIR}/${SID}_normal.bam.bai" ]; then
        note "copying the normal BAM, untouched"
        cp "${NORMAL}" "${LDIR}/${SID}_normal.bam"
        cp "${NORMAL}.bai" "${LDIR}/${SID}_normal.bam.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${LDIR}/${SID}_normal.bam"
    fi
    if [ ! -s "${LDIR}/${SID}_tumour_original.bam.bai" ]; then
        note "copying the tumour BAM as a pre-loss reference"
        cp "${TUMOUR}" "${LDIR}/${SID}_tumour_original.bam"
        cp "${TUMOUR}.bai" "${LDIR}/${SID}_tumour_original.bam.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${LDIR}/${SID}_tumour_original.bam"
    fi

    SRC="${LDIR}/${SID}_tumour_original.bam"
    DA0=$(depth_of "${SRC}" ${A_S} ${A_E})
    DB0=$(depth_of "${SRC}" ${B_S} ${B_E})
    DC0=$(depth_of "${SRC}" ${C_S} ${C_E})

    # ---- everything outside the target, plus a thinned target ----
    note "$(date '+%T') thinning HLA-${TARGET} to ${KEEP}"
    samtools view -h "${SRC}" \
      | awk -v s="${T_S}" -v e="${T_E}" \
            'BEGIN{OFS="\t"} /^@/ {print; next}
             !($3=="chr6" && $4>=s && $4<=e) {print}' \
      | samtools view -b -o "${LDIR}/rest.bam" - 2>/dev/null

    samtools view -b -s "${KEEP}" "${SRC}" "chr6:${T_S}-${T_E}" \
        > "${LDIR}/target_thin.bam" 2>/dev/null

    samtools merge -f -@ "${THREADS}" -o "${LDIR}/merged.bam" \
        "${LDIR}/rest.bam" "${LDIR}/target_thin.bam" 2>/dev/null
    samtools sort -@ "${THREADS}" -m 1G \
        -o "${LDIR}/${SID}_tumour_LOH.bam" "${LDIR}/merged.bam" 2>/dev/null
    samtools index -@ "${THREADS}" "${LDIR}/${SID}_tumour_LOH.bam"
    rm -f "${LDIR}/rest.bam" "${LDIR}/target_thin.bam" "${LDIR}/merged.bam"

    if [ ! -s "${LDIR}/${SID}_tumour_LOH.bam" ]; then
        fail "could not build the LOH BAM"
        FAIL=$((FAIL+1)); continue
    fi

    LOHBAM="${LDIR}/${SID}_tumour_LOH.bam"
    DA1=$(depth_of "${LOHBAM}" ${A_S} ${A_E})
    DB1=$(depth_of "${LOHBAM}" ${B_S} ${B_E})
    DC1=$(depth_of "${LOHBAM}" ${C_S} ${C_E})

    # ---- what LOHHLA needs beside the BAMs ----
    # all six alleles, named the way LOHHLA writes them
    awk 'NR==2 {
      for (i = 2; i <= 7; i++) {
        if ($i ~ /\*/) {
          split($i, x, "*"); split(x[2], p, ":")
          printf "hla_%s_%s_%s\n", tolower(x[1]), p[1], p[2]
        }
      }
    }' "${HLA}" > "${LDIR}/${SID}_alleles.txt"

    # purity and ploidy are exact by construction here: nothing was
    # diluted and no whole-genome copy number was altered
    printf "Ploidy\ttumorPurity\ttumorPloidy\t\n%s_tumour_LOH\t2\t1.0\t2\t\n" \
        "${SID}" > "${LDIR}/${SID}_solutions.txt"

    printf "%s\tHLA-%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${SID}" "${TARGET}" "${T_ALLELES}" "${CONTROLS}" "${KEEP}" \
        "${DA0}" "${DA1}" "${DB0}" "${DB1}" "${DC0}" "${DC1}" >> "${DESIGN}"

    T1=$(date +%s)
    printf "  built in %dm%02ds  |  %s\n" \
        $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
        "$(ls -lh ${LOHBAM} | awk '{print $5}')"
    for L in A B C; do
        eval "b=\$D${L}0; a=\$D${L}1"
        printf "   HLA-%s  %7s -> %7s x%s\n" "${L}" "${b}" "${a}" \
            "$([ "${TARGET}" = "${L}" ] && echo '   <- target' || echo '   control')"
    done
    OK=$((OK+1))
done

step_header "STEP 4 done"
note "built ${OK}   skipped ${SKIP}   homozygous ${HOMO}   waiting ${WAIT}   failed ${FAIL}"

if [ -s "${DESIGN}" ]; then
    say "design"
    column -t -s $'\t' "${DESIGN}" 2>/dev/null | sed 's/^/  /' || cat "${DESIGN}"
fi

[ "${FAIL}" -gt 0 ] && exit 1
exit 0
