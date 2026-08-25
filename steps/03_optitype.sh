#!/bin/bash
#SBATCH --job-name=demo_optitype
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.log
#
# Step 3 — HLA typing.
#
# Ported from the cohort pipeline's step 4.
#
# The typing runs on the normal BAM, not the tumour. A person has one HLA
# genotype and the injected mutations do not change it, so the untouched
# reads are the honest source. It also means this step does not have to
# wait for the injection.
#
# Reads come from three windows on chromosome 6 rather than from the whole
# file. OptiType aligns every read against the full IMGT allele set, and
# handing it a panel-wide BAM makes that far slower for no gain.
#
# Only properly paired primary alignments are taken. OptiType pairs reads
# by name and a singleton or a secondary alignment either confuses that or
# contributes a fragment of evidence that is not independent.
#
# Threads matter more here than anywhere else in the pipeline. Given fewer
# than it wants, OptiType does not fail — it returns a plausible default
# genotype, and an earlier run produced fifteen identical genotypes in a
# row before anyone noticed. The demo counts distinct genotypes at the end
# for exactly this reason.
#
# Usage:
#   sbatch 03_optitype.sh              every sample in config.sh
#   sbatch 03_optitype.sh B003 B007    named samples
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

TARGETS="${*:-${SAMPLES}}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# hg38 windows for the class I loci, with flanks: OptiType realigns whole
# regions against allele sequences and exon-level intervals leave it with
# fragmentary coverage
HLA_A="chr6:29932000-29956000"
HLA_B="chr6:31343000-31367000"
HLA_C="chr6:31258000-31282000"

step_header "STEP 3 — HLA typing"
note "samples   $(echo ${TARGETS} | tr '\n' ' ')"
note "threads   ${THREADS}"
note "started   $(date '+%F %T') on $(hostname)"

TYPED=0; SKIP=0; FAIL=0

for SID in ${TARGETS}; do
    OUT="${OUT_DIR}/${SID}"
    NORMAL="${SAMPLE_DIR}/${SID}_normal.bam"
    ODIR="${OUT}/optitype"
    RESULT="${ODIR}/${SID}_result.tsv"

    echo ""
    echo "----------------------------------------------------------------------"
    echo " ${SID}"

    if [ ! -s "${NORMAL}" ]; then
        fail "normal BAM missing"
        FAIL=$((FAIL+1)); continue
    fi
    if [ -s "${RESULT}" ]; then
        note "already typed — skipping"
        SKIP=$((SKIP+1)); continue
    fi

    mkdir -p "${ODIR}" "${OUT}/hla_fastq"
    T0=$(date +%s)

    # ---- collect the reads over the three loci ----
    use_env "${ENV_GATK}" || exit 1

    note "$(date '+%T') collecting reads over HLA-A, -B and -C"
    samtools view -b -f 2 -F 0x900 "${NORMAL}" \
        "${HLA_A}" "${HLA_B}" "${HLA_C}" \
        > "${OUT}/hla_fastq/hla.bam" 2>/dev/null

    # OptiType pairs by read name, so the BAM is sorted by name before the
    # fastq is written; coordinate order would separate the mates
    samtools sort -n -@ "${THREADS}" -m 1G \
        -o "${OUT}/hla_fastq/hla.qsort.bam" \
        "${OUT}/hla_fastq/hla.bam" 2>/dev/null

    samtools fastq -@ "${THREADS}" \
        -1 "${OUT}/hla_fastq/${SID}_1.fastq" \
        -2 "${OUT}/hla_fastq/${SID}_2.fastq" \
        -0 /dev/null -s /dev/null -n \
        "${OUT}/hla_fastq/hla.qsort.bam" 2>/dev/null

    rm -f "${OUT}/hla_fastq/hla.bam" "${OUT}/hla_fastq/hla.qsort.bam"

    N1=$(( $(wc -l < "${OUT}/hla_fastq/${SID}_1.fastq") / 4 ))
    N2=$(( $(wc -l < "${OUT}/hla_fastq/${SID}_2.fastq") / 4 ))
    note "${N1} read pairs"

    if [ "${N1}" -lt 100 ]; then
        fail "too few reads to type — ${N1} pairs"
        FAIL=$((FAIL+1)); continue
    fi
    if [ "${N1}" -ne "${N2}" ]; then
        fail "mates do not match: ${N1} against ${N2}"
        FAIL=$((FAIL+1)); continue
    fi

    # ---- run OptiType ----
    use_env "${ENV_OPTITYPE}" || exit 1

    # OptiType needs a config naming the razers3 binary by absolute path;
    # it is written per run so that the path follows the environment
    # rather than being fixed at install time
    RAZERS="$(command -v razers3)"
    if [ -z "${RAZERS}" ]; then
        fail "razers3 not found in ${ENV_OPTITYPE}"
        FAIL=$((FAIL+1)); continue
    fi

    cat > "${ODIR}/config.ini" << EOF
[mapping]
razers3=${RAZERS}
threads=${THREADS}

[ilp]
solver=glpk
threads=1

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
EOF

    note "$(date '+%T') OptiTypePipeline.py"
    if ! OptiTypePipeline.py \
            -i "${OUT}/hla_fastq/${SID}_1.fastq" \
               "${OUT}/hla_fastq/${SID}_2.fastq" \
            --dna \
            -c "${ODIR}/config.ini" \
            -o "${ODIR}" \
            -p "${SID}" \
            > "${ODIR}/optitype.log" 2>&1; then
        fail "OptiTypePipeline.py failed:"
        tail -6 "${ODIR}/optitype.log" | sed 's/^/     /'
        FAIL=$((FAIL+1)); continue
    fi

    # OptiType writes into a timestamped subdirectory when -p is not
    # honoured by the version in use; the result is found either way
    if [ ! -s "${RESULT}" ]; then
        found="$(find "${ODIR}" -name "*_result.tsv" -newermt "-1 hour" \
                 2>/dev/null | head -1)"
        [ -n "${found}" ] && cp "${found}" "${RESULT}"
    fi

    if [ ! -s "${RESULT}" ]; then
        fail "OptiType produced no result file"
        FAIL=$((FAIL+1)); continue
    fi

    rm -f "${OUT}/hla_fastq/${SID}"_[12].fastq

    T1=$(date +%s)

    GT="$(awk -F'\t' 'NR==2 {printf "A %s/%s  B %s/%s  C %s/%s",
                             $2, $3, $4, $5, $6, $7}' "${RESULT}")"
    READS="$(awk -F'\t' 'NR==2 {print $8}' "${RESULT}")"

    printf "  %dm%02ds  |  %s\n" \
        $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${GT}"
    note "  ${READS} reads reached the solver"
    TYPED=$((TYPED+1))
done

step_header "STEP 3 done"
note "typed ${TYPED}   skipped ${SKIP}   failed ${FAIL}"

# ---- the check that matters ----
# Identical genotypes across samples mean OptiType defaulted rather than
# typed. Forty backbones were chosen so that none is reused, so every
# genotype should differ.
say "distinct genotypes"
n_res=0; 
tmp="$(mktemp)"
for SID in ${SAMPLES}; do
    r="${OUT_DIR}/${SID}/optitype/${SID}_result.tsv"
    [ -s "${r}" ] || continue
    n_res=$((n_res+1))
    awk -F'\t' 'NR==2 {print $2"|"$3"|"$4"|"$5"|"$6"|"$7}' "${r}" >> "${tmp}"
done
n_uniq=$(sort -u "${tmp}" | wc -l)
rm -f "${tmp}"

note "${n_uniq} distinct among ${n_res} typed"
if [ "${n_res}" -gt 1 ] && [ "${n_uniq}" -lt "${n_res}" ]; then
    fail "two or more samples share a genotype — OptiType may have"
    fail "defaulted; check that it had the threads it asked for"
fi

[ "${FAIL}" -gt 0 ] && exit 1
exit 0
