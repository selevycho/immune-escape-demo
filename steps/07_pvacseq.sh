#!/bin/bash
#SBATCH --job-name=demo_pvac
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.log
#
# Step 7 — neoantigen prediction with pVACseq and NetMHCpan.
#
# Ported from the cohort pipeline's step 7.
#
# This is the second, independent route to the answer step 5 already gave.
# The mhcflurry path took mutations from the truth file with the protein
# consequence already annotated and cut peptides straight out of GENCODE
# sequences. Here the Mutect2 VCF goes through VEP, which recomputes the
# consequence from the genome, and pVACseq builds the peptides itself
# before handing them to NetMHCpan-4.1.
#
# Running both on the same sample is the point. They share no code and
# start from different files — one from what was placed, the other from
# what was called — so agreement is evidence that neither carries a
# systematic error, and disagreement localises where one of them does.
#
# Three details that cost time to find:
#
#   Sample order in the VCF. Mutect2 writes the normal in column 10 and
#   the tumour in column 11. Passing them the wrong way round produces a
#   valid-looking run with an empty result, since the "tumour" then has no
#   somatic variants at all.
#
#   Duplicate HLA alleles. A homozygous locus appears twice in the
#   OptiType output; pVACseq wants each allele once.
#
#   VEP needs its plugins on an explicit --dir_plugins path, and VEP and
#   pVACseq live in separate environments because VEP's perl stack and
#   pvactools' tensorflow cannot be resolved together.
#
# Without netMHCpan this step is skipped rather than failed. It is
# licensed per user by DTU and often absent; everything else still runs.
#
# Usage:
#   sbatch 07_pvacseq.sh              every sample in config.sh
#   sbatch 07_pvacseq.sh B003 B007    named samples
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
CPUS="${SLURM_CPUS_PER_TASK:-4}"
ALGORITHMS="${ALGORITHMS:-NetMHCpan}"
EPITOPE_LENGTHS="${EPITOPE_LENGTHS:-8,9,10,11}"
VEP_PLUGINS="${VEP_PLUGINS:-${HOME}/vep_plugins}"

step_header "STEP 7 — pVACseq with ${ALGORITHMS}"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "started   $(date '+%F %T') on $(hostname)"

if [ ! -d "${VEP_CACHE}/homo_sapiens" ]; then
    fail "VEP cache missing — run setup/fetch_data.sh without --no-vep"
    exit 1
fi
if [ ! -x "${NETMHC_DIR}/netMHCpan" ]; then
    note ""
    note "netMHCpan is not present at ${NETMHC_DIR}"
    note ""
    note "It cannot be redistributed: DTU licenses it per user. Request"
    note "one at https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/"
    note "unpack the archive there, and edit its netMHCpan script to set"
    note "NMHOME to that path."
    note ""
    note "This step is skipped. Step 5 already predicted binding with"
    note "mhcflurry, so the demo still reports neoantigens — it simply"
    note "cannot cross-check them against a second engine."
    exit 0
fi
if [ ! -s "${VEP_PLUGINS}/Wildtype.pm" ]; then
    fail "VEP plugins missing at ${VEP_PLUGINS}"
    note "pVACseq needs the Wildtype and Frameshift plugins; set"
    note "VEP_PLUGINS to where they are installed"
    exit 1
fi

OK=0; SKIP=0; WAIT=0; FAIL=0
T_ALL=$(date +%s)

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    VCF="${OUT}/mutect2/${SID}.filtered.vcf.gz"
    HLA_TSV="${OUT}/optitype/${SID}_result.tsv"
    PDIR="${OUT}/pvacseq"
    ANNOT="${PDIR}/${SID}.vep.vcf"
    RESULT_DIR="${PDIR}/out"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${VCF}" ] || [ ! -s "${HLA_TSV}" ]; then
        note "VCF or HLA type not ready"
        WAIT=$((WAIT+1)); continue
    fi

    FINAL=$(ls "${RESULT_DIR}"/MHC_Class_I/*.all_epitopes.tsv 2>/dev/null | head -1)
    if [ -n "${FINAL}" ] && [ "$(wc -l < ${FINAL})" -gt 1 ]; then
        note "already done: $(( $(wc -l < ${FINAL}) - 1 )) epitopes"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${PDIR}"
    T0=$(date +%s)

    # HLA alleles, each once, in pVACseq notation
    ALLELES=$(awk 'NR==2 {
        n = 0
        for (i = 2; i <= 7; i++)
            if ($i ~ /\*/ && !seen[$i]++) {
                if (n++) printf ","
                printf "HLA-%s", $i
            }
    }' "${HLA_TSV}")
    note "HLA: ${ALLELES}"

    if [ -z "${ALLELES}" ]; then
        fail "no usable HLA alleles"
        FAIL=$((FAIL+1)); continue
    fi

    # ---- 1. VEP ----
    if [ ! -s "${ANNOT}" ]; then
        use_env "${ENV_VEP}" || exit 1
        note "[1/2] $(date '+%T') VEP"
        vep \
            --input_file "${VCF}" \
            --output_file "${ANNOT}" \
            --format vcf --vcf \
            --symbol --terms SO --tsl --biotype --hgvs \
            --fasta "${REF_FASTA}" \
            --offline --cache --dir_cache "${VEP_CACHE}" \
            --plugin Frameshift --plugin Wildtype \
            --dir_plugins "${VEP_PLUGINS}" \
            --pick --transcript_version \
            --fork "${CPUS}" --force_overwrite \
            > "${PDIR}/vep.log" 2>&1
        RC=$?
        if [ "${RC}" -ne 0 ] || [ ! -s "${ANNOT}" ]; then
            fail "VEP failed (exit ${RC}):"
            tail -6 "${PDIR}/vep.log" | sed 's/^/     /'
            FAIL=$((FAIL+1)); continue
        fi
        note "      $(grep -vc '^#' ${ANNOT}) records annotated"
    else
        note "[1/2] VEP output already present"
    fi

    # ---- 2. pVACseq ----
    # Mutect2 column order: 10 = normal, 11 = tumour
    NORMAL_NAME=$(grep '^#CHROM' "${ANNOT}" | cut -f10)
    TUMOUR_NAME=$(grep '^#CHROM' "${ANNOT}" | cut -f11)
    note "[2/2] $(date '+%T') pVACseq: tumour ${TUMOUR_NAME}, normal ${NORMAL_NAME}"

    use_env "${ENV_PVACSEQ}" || exit 1
    export PATH="${NETMHC_DIR}:${PATH}"
    rm -rf "${RESULT_DIR}"

    pvacseq run \
        "${ANNOT}" \
        "${TUMOUR_NAME}" \
        "${ALLELES}" \
        ${ALGORITHMS} \
        "${RESULT_DIR}" \
        -e1 "${EPITOPE_LENGTHS}" \
        --normal-sample-name "${NORMAL_NAME}" \
        --pass-only \
        --normal-cov 5 --tdna-cov 5 \
        --normal-vaf 0.02 --tdna-vaf 0.05 \
        --maximum-transcript-support-level 5 \
        --n-threads "${CPUS}" \
        > "${PDIR}/pvacseq.log" 2>&1
    RC=$?
    T1=$(date +%s)

    ALL=$(ls "${RESULT_DIR}"/MHC_Class_I/*.all_epitopes.tsv 2>/dev/null | head -1)
    FILT=$(ls "${RESULT_DIR}"/MHC_Class_I/*.filtered.tsv 2>/dev/null | head -1)

    if [ -z "${ALL}" ] || [ "$(wc -l < ${ALL})" -le 1 ]; then
        fail "no epitopes produced (exit ${RC})"
        grep -iE "error|exception|empty" "${PDIR}/pvacseq.log" | tail -4 \
            | sed 's/^/     /'
        FAIL=$((FAIL+1)); continue
    fi

    N_ALL=$(( $(wc -l < ${ALL}) - 1 ))
    N_FILT=0
    [ -n "${FILT}" ] && N_FILT=$(( $(wc -l < ${FILT}) - 1 ))
    printf "      %dm%02ds   %d epitopes, %d after pVACseq's own filter\n" \
        $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${N_ALL}" "${N_FILT}"

    # Binders are counted from all_epitopes at the same percentile
    # thresholds step 5 used, so the two routes are directly comparable.
    #
    # pVACseq's filtered.tsv is not used for this. Its chain — binding,
    # coverage, transcript support, then one best epitope per mutation —
    # shortlists vaccine candidates and reduced 49 690 epitopes to three
    # on the first sample it was tried on. That is the right answer to a
    # different question.
    awk -F'\t' 'NR==1 {
        for (i=1;i<=NF;i++) if ($i == "Best MT Percentile") col=i
        next
    }
    col && $col != "NA" {
        if ($col+0 < 0.5) sb++
        else if ($col+0 < 2) wb++
    }
    END { printf "      strong (<0.5%%): %d   weak (0.5-2%%): %d\n", sb+0, wb+0 }' \
        "${ALL}"

    AGG=$(ls "${RESULT_DIR}"/MHC_Class_I/*.aggregated.tsv 2>/dev/null | head -1)
    [ -n "${AGG}" ] && note "      mutations with an epitope: $(( $(wc -l < ${AGG}) - 1 ))"
    OK=$((OK+1))
done

T_END=$(date +%s)
step_header "STEP 7 done"
printf "  elapsed %dh%02dm   done %d   skipped %d   waiting %d   failed %d\n" \
    $(( (T_END-T_ALL)/3600 )) $(( ((T_END-T_ALL)%3600)/60 )) \
    "${OK}" "${SKIP}" "${WAIT}" "${FAIL}"

say "the two routes side by side"
printf "  %-8s %12s %12s %12s %12s\n" \
       "sample" "mhcflurry" "NetMHCpan" "mf strong" "pv strong"
for SID in ${SAMPLES}; do
    MF="${OUT_DIR}/${SID}/neoantigens/neoantigens_per_mutation.tsv"
    A=$(ls "${OUT_DIR}/${SID}"/pvacseq/out/MHC_Class_I/*.all_epitopes.tsv \
        2>/dev/null | head -1)
    [ -s "${MF}" ] || [ -n "${A}" ] || continue

    mf_pep=0; mf_str=0
    if [ -s "${MF}" ]; then
        read -r mf_pep mf_str < <(awk -F'\t' 'NR>1 {p+=$4; s+=$6}
                                   END {print p+0, s+0}' "${MF}")
    fi
    pv_pep=0; pv_str=0
    if [ -n "${A}" ]; then
        pv_pep=$(( $(wc -l < ${A}) - 1 ))
        pv_str=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if ($i=="Best MT Percentile") c=i; next}
                  c && $c != "NA" && $c+0 < 0.5 {n++} END {print n+0}' "${A}")
    fi
    printf "  %-8s %12s %12s %12s %12s\n" \
        "${SID}" "${mf_pep}" "${pv_pep}" "${mf_str}" "${pv_str}"
done

echo ""
note "The counts differ because the two routes count different things:"
note "mhcflurry scores every peptide-allele pair, pVACseq one row per"
note "epitope per transcript. What is comparable is the strong-binder"
note "count and which mutations produce them."

[ "${FAIL}" -gt 0 ] && exit 1
exit 0
