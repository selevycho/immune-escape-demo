#!/bin/bash
#SBATCH --job-name=demo_lohhla
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.log
#
# Step 6 — LOHHLA, one locus per invocation.
#
# Ported from the cohort pipeline's step 6b.
#
# LOHHLA processes every locus in the allele file within a single run and
# aborts the whole run if any one of them fails — a locus whose two
# alleles differ at only a handful of covered positions gives its t-test
# constant data, and R halts, discarding the loci that had worked. Running
# each locus separately means a failure costs that locus only.
#
# The allele FASTA keeps every IMGT subtype of each called allele rather
# than one representative. Reducing it was tried and made things worse:
# the reference run that succeeded used 834 sequences for six alleles, and
# cutting that to six dropped usable mismatch positions from 95 to 3.
#
# LOHHLA here is a modified copy — novoalign replaced by bwa mem -a, GATK
# 3 jar calls by samtools equivalents, hg19 HLA coordinates by hg38 ones
# with chr prefixes, invalid multi-letter optparse short flags removed,
# and a name-sort inserted before FASTQ conversion so read pairs survive.
# The statistics are untouched. See soft/lohhla/PATCHES.md.
#
# It insists on bare BAM filenames rather than paths, and does not create
# the per-BAM working directories it expects. Both are handled below.
#
# Expect results for some samples and not others. At panel depth many loci
# have too few positions distinguishing the two alleles for the test to
# run at all; across the full cohort 23 of 40 produced any output. A blank
# result is the method's limit, not a failure of this script.
#
# Usage:
#   sbatch 06_lohhla.sh              every sample in config.sh
#   sbatch 06_lohhla.sh B003 B007    named samples
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
MIN_COV="${MIN_COV:-5}"

HLA_FASTA_ALL="${LOHHLA_DIR}/data/hla_all_lohhla.fasta"
HLA_EXON="${LOHHLA_DIR}/data/hla.dat"
DESIGN="${OUT_DIR}/loh_design.tsv"

step_header "STEP 6 — LOHHLA"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "minCov    ${MIN_COV}"
note "started   $(date '+%F %T') on $(hostname)"

for f in "${LOHHLA_DIR}/LOHHLAscript.R" "${HLA_FASTA_ALL}"; do
    [ -s "${f}" ] || { fail "missing ${f}"; exit 1; }
done
if [ ! -s "${HLA_EXON}" ]; then
    fail "missing ${HLA_EXON}"
    note "hla.dat is part of the IMGT release; setup/fetch_data.sh should"
    note "have unpacked it beside the FASTA"
    exit 1
fi

use_env "${ENV_LOHHLA}" || exit 1

# the patched script needs samtools view -N, which arrived in 1.12;
# lohhla_env is pinned to an older samtools by its R dependencies
GATK_BIN="$(conda info --base)/envs/${ENV_GATK}/bin"
[ -d "${GATK_BIN}" ] && export PATH="${GATK_BIN}:${PATH}"
note "samtools  $(samtools --version | head -1)"

N_RUN=0; N_OK=0; N_FAIL=0; N_WAIT=0; N_SKIP=0

for SID in ${TARGETS}; do
    DIR="${OUT_DIR}/${SID}/loh"
    NORMAL="${SID}_normal.bam"
    TUMOUR="${SID}_tumour_LOH.bam"
    ALLELES="${DIR}/${SID}_alleles.txt"
    SOLUTIONS="${DIR}/${SID}_solutions.txt"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${DIR}/${TUMOUR}" ] || [ ! -s "${ALLELES}" ]; then
        note "LOH BAM or allele file missing — run step 4 first"
        N_WAIT=$((N_WAIT+1)); continue
    fi

    TARGET_LOCUS=$(awk -F'\t' -v s="${SID}" 'NR>1 && $1==s {print $2}' \
                   "${DESIGN}" 2>/dev/null | head -1)
    note "locus with simulated loss: ${TARGET_LOCUS:-unknown}"

    cd "${DIR}"

    for LOCUS in a b c; do
        # only heterozygous loci are worth running; LOHHLA refuses the rest
        grep "^hla_${LOCUS}_" "${ALLELES}" | sort -u > "${DIR}/al_${LOCUS}.txt"
        N_AL=$(wc -l < "${DIR}/al_${LOCUS}.txt")
        if [ "${N_AL}" -lt 2 ]; then
            note "HLA-${LOCUS^^}: homozygous or untyped — skipping"
            rm -f "${DIR}/al_${LOCUS}.txt"
            N_SKIP=$((N_SKIP+1)); continue
        fi

        OUTDIR="${DIR}/out_${LOCUS}"
        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            note "HLA-${LOCUS^^}: already done"
            N_OK=$((N_OK+1)); continue
        fi

        # every IMGT subtype of these two alleles — see the note above
        SUB="${DIR}/sub_${LOCUS}.fasta"
        rm -f "${SUB}"*
        python3 - "${DIR}/al_${LOCUS}.txt" "${HLA_FASTA_ALL}" "${SUB}" << 'PYEOF'
import sys
keep = tuple({l.strip() for l in open(sys.argv[1]) if l.strip()})
out, w, n = [], False, 0
for line in open(sys.argv[2]):
    if line.startswith(">"):
        w = any(line[1:].strip().startswith(k) for k in keep)
        if w:
            n += 1
    if w:
        out.append(line)
open(sys.argv[3], "w").writelines(out)
print("       FASTA: %d sequences" % n)
PYEOF

        if [ "$(grep -c '^>' ${SUB})" -lt 2 ]; then
            note "HLA-${LOCUS^^}: fewer than two sequences matched — skipping"
            N_FAIL=$((N_FAIL+1)); continue
        fi

        rm -rf "${OUTDIR}"
        mkdir -p "${OUTDIR}/${NORMAL%.bam}" "${OUTDIR}/${TUMOUR%.bam}"

        note "HLA-${LOCUS^^}: $(date '+%T') running"
        N_RUN=$((N_RUN+1))
        T0=$(date +%s)

        Rscript "${LOHHLA_DIR}/LOHHLAscript.R" \
            --patientId "${SID}" \
            --outputDir "${OUTDIR}" \
            --normalBAMfile "${NORMAL}" \
            --tumorBAMfile "${TUMOUR}" \
            --hlaPath "${DIR}/al_${LOCUS}.txt" \
            --HLAfastaLoc "${SUB}" \
            --HLAexonLoc "${HLA_EXON}" \
            --CopyNumLoc "${SOLUTIONS}" \
            --mappingStep TRUE \
            --fishingStep FALSE \
            --coverageStep TRUE \
            --plottingStep TRUE \
            --cleanUp FALSE \
            --minCoverageFilter "${MIN_COV}" \
            --numMisMatch 1 \
            --ignoreWarnings TRUE \
            --novoDir "" --gatkDir "" \
            > "${DIR}/lohhla_${LOCUS}.log" 2> "${DIR}/lohhla_${LOCUS}.err"
        RC=$?
        T1=$(date +%s)

        PRED=$(ls "${OUTDIR}"/*HLAlossPrediction*.txt 2>/dev/null | head -1)
        if [ -n "${PRED}" ] && [ "$(wc -l < ${PRED})" -gt 1 ]; then
            printf "       done in %dm%02ds\n" \
                $(( (T1-T0)/60 )) $(( (T1-T0)%60 ))
            awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
              { printf "       %s vs %s   CN %.3f / %.3f   P=%.4g   sites=%s\n",
                $h["HLA_A_type1"], $h["HLA_A_type2"],
                $h["HLA_type1copyNum_withBAFBin"],
                $h["HLA_type2copyNum_withBAFBin"],
                $h["PVal"], $h["numMisMatchSitesCov"] }' "${PRED}" 2>/dev/null
            N_OK=$((N_OK+1))
        else
            note "       no prediction (exit ${RC}): $(tail -2 ${DIR}/lohhla_${LOCUS}.err | head -1)"
            printf "%s\thla_%s\trc%s\n" "${SID}" "${LOCUS}" "${RC}" \
                >> "${OUT_DIR}/lohhla_failures.log"
            N_FAIL=$((N_FAIL+1))
        fi
    done
done

step_header "STEP 6 done"
note "locus runs ${N_RUN}   succeeded ${N_OK}   no result ${N_FAIL}   homozygous ${N_SKIP}   waiting ${N_WAIT}"

say "predictions"
printf " %-8s %-7s %-20s %-20s %8s %8s %10s %6s %s\n" \
       "sample" "locus" "allele 1" "allele 2" "CN1" "CN2" "PVal" "sites" "role"
FOUND=0
for d in "${OUT_DIR}"/*/loh/out_*/; do
    P=$(ls "${d}"*HLAlossPrediction*.txt 2>/dev/null | head -1)
    [ -n "${P}" ] && [ "$(wc -l < ${P})" -gt 1 ] || continue
    s=$(echo "$d" | sed "s|${OUT_DIR}/||; s|/loh/out_.*||")
    loc=$(echo "$d" | sed 's|.*/out_||; s|/||')
    tgt=$(awk -F'\t' -v x="$s" -v l="HLA-${loc^^}" \
          'NR>1 && $1==x && $2==l {print "yes"}' "${DESIGN}" 2>/dev/null | head -1)
    awk -F'\t' -v s="$s" -v loc="HLA-${loc^^}" -v tgt="${tgt:-control}" '
      NR==1 {for (i=1;i<=NF;i++) h[$i]=i; next}
      { printf " %-8s %-7s %-20s %-20s %8.3f %8.3f %10.4g %6s %s\n",
        s, loc, $h["HLA_A_type1"], $h["HLA_A_type2"],
        $h["HLA_type1copyNum_withBAFBin"], $h["HLA_type2copyNum_withBAFBin"],
        $h["PVal"], $h["numMisMatchSitesCov"],
        (tgt == "yes" ? "TARGET" : "control") }' "${P}" 2>/dev/null
    FOUND=$((FOUND+1))
done

if [ "${FOUND}" -eq 0 ]; then
    echo ""
    note "No locus produced a prediction. At this depth that is the"
    note "expected outcome for many samples: the two alleles of a locus"
    note "must differ at enough covered positions for the test to have"
    note "anything to compare, and panel data often does not provide it."
fi

exit 0
