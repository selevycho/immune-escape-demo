#!/bin/bash
#SBATCH --job-name=demo_addindel
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#
# Step 1b — inject indels into the substituted BAM.
#
# Ported from the cohort pipeline's step 8. Two things learned there are
# preserved exactly.
#
# addindel.py needs its input under a name it can index beside, so the
# substituted BAM is copied to a working name together with its index
# rather than being passed in place.
#
# addindel.py never emits a coordinate-sorted BAM. Indexing one fails with
# a message about the index rather than about the order, so the sort is
# not optional.
#
# Indels inside the MHC are dropped up front. Realignment there returns
# reads with zero mapping quality because hg38 carries alternate MHC
# haplotypes, and a position that could never have worked should not be
# counted as a failure.
#
# A sample with no indels still produces a tumour BAM: the substituted one
# is carried forward under the final name, so every downstream step finds
# what it expects.
#
# Usage:
#   sbatch 01b_addindel.sh              every sample in config.sh
#   sbatch 01b_addindel.sh B003 B007    named samples
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
THREADS="${SLURM_CPUS_PER_TASK:-2}"

MHC_S=29600000
MHC_E=33100000

step_header "STEP 1b — indels"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "threads   ${THREADS}"
note "started   $(date '+%F %T') on $(hostname)"

OK=0; SKIP=0; FAIL=0; NONE=0; WAIT=0

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    TRUTH="${SAMPLE_DIR}/${SID}_truth.tsv"
    SRC="${OUT}/${SID}_tumour_snv.bam"
    FINAL="${OUT}/${SID}_tumour.bam"
    WORK="${OUT}/${SID}_input.bam"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${SRC}" ]; then
        note "substituted BAM not ready — step 1a has not reached this sample"
        WAIT=$((WAIT+1)); continue
    fi
    if [ -s "${FINAL}.bai" ]; then
        note "already built — skipping"
        SKIP=$((SKIP+1)); continue
    fi

    # BAMSurgeon variant format for indels:
    #   DEL: chrom  start  end  VAF  DEL
    #   INS: chrom  start  start  VAF  INS  sequence
    # For a deletion the MAF gives the deleted bases as the reference
    # allele, so the end coordinate follows from its length.
    # Columns by name, not by number — see the note in 01a_addsnv.sh.
    awk -F'\t' -v ms="${MHC_S}" -v me="${MHC_E}" '
        NR==1 { for (i=1; i<=NF; i++) c[$i]=i; next }
        {
            chrom = $c["chrom"]; pos = $c["pos"]
            ref = $c["ref"]; alt = $c["alt"]; vaf = $c["vaf"]
            if (chrom=="chr6" && pos>=ms && pos<me) next
        }
        $c["type"]=="DEL" {
            end = pos + length(ref) - 1
            if (end <= pos) end = pos + 1
            printf "%s\t%s\t%d\t%s\tDEL\n", chrom, pos, end, vaf
        }
        $c["type"]=="INS" {
            seq = (alt == "-" || alt == "") ? "A" : alt
            printf "%s\t%s\t%s\t%s\tINS\t%s\n", chrom, pos, pos, vaf, seq
        }
    ' "${TRUTH}" > "${OUT}/indels.txt"

    N_TOTAL=$(awk -F'\t' '
        NR==1 { for (i=1; i<=NF; i++) c[$i]=i; next }
        $c["type"]=="DEL" || $c["type"]=="INS"' "${TRUTH}" | wc -l)
    N_IND=$(wc -l < "${OUT}/indels.txt")
    N_MHC=$(( N_TOTAL - N_IND ))
    note "${N_TOTAL} indels in the truth set, ${N_MHC} dropped from the MHC, ${N_IND} to inject"

    if [ "${N_IND}" -lt 1 ]; then
        note "nothing to inject — carrying the substituted BAM forward"
        cp "${SRC}" "${FINAL}"
        cp "${SRC}.bai" "${FINAL}.bai" 2>/dev/null
        NONE=$((NONE+1)); continue
    fi

    T0=$(date +%s)

    # ---- a working copy under a name addindel can index beside ----
    use_env "${ENV_GATK}" || exit 1
    if [ ! -s "${WORK}.bai" ]; then
        note "copying the substituted BAM ..."
        cp "${SRC}" "${WORK}"
        cp "${SRC}.bai" "${WORK}.bai" 2>/dev/null \
            || samtools index -@ "${THREADS}" "${WORK}"
    fi

    use_env "${ENV_BAMSURGEON}" || exit 1
    PICARD=$(find "${CONDA_PREFIX}" -name "picard.jar" 2>/dev/null | head -1)
    if [ -z "${PICARD}" ]; then
        fail "picard.jar not found"
        FAIL=$((FAIL+1)); continue
    fi

    cd "${OUT}"
    rm -rf addindel.tmp
    note "$(date '+%T') addindel.py on ${N_IND} indels ..."

    if ! addindel.py \
            -v "${OUT}/indels.txt" \
            -f "${WORK}" -r "${REF_FASTA}" \
            -o "${OUT}/raw.bam" \
            -p "${THREADS}" --aligner mem \
            --picardjar "${PICARD}" --force --insane \
            > "${OUT}/addindel.log" 2>&1; then
        fail "addindel.py failed:"
        tail -6 "${OUT}/addindel.log" | sed 's/^/     /'
        rm -rf "${OUT}/raw.bam"* "${OUT}/addindel.tmp"
        FAIL=$((FAIL+1)); continue
    fi

    # addindel never emits a coordinate-sorted BAM
    use_env "${ENV_GATK}" || exit 1
    samtools sort -@ "${THREADS}" -m 1G -o "${FINAL}" "${OUT}/raw.bam" 2>/dev/null
    samtools index -@ "${THREADS}" "${FINAL}"
    rm -rf "${OUT}/raw.bam"* "${OUT}/addindel.tmp" "${OUT}"/addindel_logs_*

    if [ ! -s "${FINAL}" ]; then
        fail "sorting failed"
        FAIL=$((FAIL+1)); continue
    fi

    T1=$(date +%s)
    printf "  injected in %dm%02ds  |  %s\n" \
           $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) \
           "$(ls -lh ${FINAL} | awk '{print $5}')"

    # ---- verify, reading the pileup the way indels are actually written ----
    # mpileup marks a deletion in the column before the deleted base, as
    # -1G, and shows the base itself as *. Looking only at the deletion's
    # own coordinate finds nothing.
    HIT=0
    while IFS=$'\t' read -r chrom start end vaf typ seq; do
        if [ "${typ}" = "DEL" ]; then
            win="${chrom}:$((start-1))-$((start+1))"
        else
            win="${chrom}:${start}-${start}"
        fi
        n=$(samtools mpileup -B -q 0 -Q 0 -r "${win}" -f "${REF_FASTA}" \
                "${FINAL}" 2>/dev/null \
            | awk '{print $5}' | grep -o '[+-][0-9]\+\|\*' | wc -l)
        [ "${n}" -gt 0 ] && HIT=$((HIT+1))
    done < "${OUT}/indels.txt"

    printf "  verified: %d / %d landed (%d%%)\n" \
           "${HIT}" "${N_IND}" $(( 100 * HIT / N_IND ))

    rm -f "${WORK}" "${WORK}.bai"
    OK=$((OK+1))
done

step_header "STEP 1b done"
note "injected ${OK}   skipped ${SKIP}   no indels ${NONE}   waiting ${WAIT}   failed ${FAIL}"
[ "${FAIL}" -gt 0 ] && exit 1
exit 0
