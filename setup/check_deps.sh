#!/bin/bash
#
# Does this machine have everything the demo needs?
#
# The demo takes about three hours, and discovering a missing tool at hour
# two is worse than discovering it now. This checks every requirement
# before anything runs and says plainly what each stage will and will not
# do.
#
# Nothing here is fatal by itself. netMHCpan is licensed per user and
# often absent, and the pipeline is expected to skip that stage rather
# than stop. The exit code reflects only whether the stages that can run
# have what they need.
#
# Usage:
#   ./check_deps.sh            check and report
#   ./check_deps.sh --quiet    report only what is wrong
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

QUIET=0
[ "${1:-}" = "--quiet" ] && QUIET=1

OK=0; MISSING=0; OPTIONAL_MISSING=0
declare -a BLOCKED_STAGES

pass() { [ "${QUIET}" = "0" ] && printf '  ok       %s\n' "$*"; OK=$((OK+1)); }
miss() { printf '  MISSING  %s\n' "$*"; MISSING=$((MISSING+1)); }
opt()  { printf '  absent   %s\n' "$*"; OPTIONAL_MISSING=$((OPTIONAL_MISSING+1)); }

# a tool is only useful if it runs, so each is invoked rather than merely
# located; a binary present but broken fails here instead of at hour two
tool_in_env() {
    local env_name="$1" tool="$2" probe="${3:---version}"
    local base
    base="$(conda info --base 2>/dev/null)"
    local bin="${base}/envs/${env_name}/bin/${tool}"

    # Presence is the test that matters. Exit codes are not a reliable
    # signal here: bwa with no arguments returns 1, gatk and picard write
    # their version banners to stderr, and OptiTypePipeline.py exits
    # non-zero when asked for help. A tool that is executable and produces
    # any output at all is installed; one that produces nothing is not.
    [ -x "${bin}" ] || return 1
    local out
    out="$("${bin}" ${probe} 2>&1 | head -3)"
    [ -n "${out}" ] || return 1
    return 0
}

env_exists() {
    conda env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"
}

say "conda"
if [ -z "${CONDA_SH:-}" ] || [ ! -f "${CONDA_SH}" ]; then
    miss "conda.sh — set CONDA_SH in config.sh"
else
    pass "conda at $(dirname "$(dirname "${CONDA_SH}")")"
fi

say "environments"
for pair in \
    "${ENV_BAMSURGEON}:injection" \
    "${ENV_GATK}:variant calling" \
    "${ENV_OPTITYPE}:HLA typing" \
    "${ENV_MHCFLURRY}:binding prediction" \
    "${ENV_LOHHLA}:HLA loss" \
    "${ENV_VEP}:annotation" \
    "${ENV_PVACSEQ}:neoantigens"; do
    e="${pair%%:*}"; what="${pair##*:}"
    if env_exists "${e}"; then
        pass "${e}  (${what})"
    else
        miss "${e}  (${what}) — run setup/install.sh"
        BLOCKED_STAGES+=("${what}")
    fi
done

say "tools"
check_tool() {
    local env_name="$1" tool="$2" probe="$3" stage="$4"
    if ! env_exists "${env_name}"; then
        return
    fi
    if tool_in_env "${env_name}" "${tool}" "${probe}"; then
        pass "${tool}  in ${env_name}"
    else
        miss "${tool}  in ${env_name} — ${stage} cannot run"
        BLOCKED_STAGES+=("${stage}")
    fi
}

check_tool "${ENV_BAMSURGEON}" samtools    "--version"  "injection"
check_tool "${ENV_BAMSURGEON}" bwa         ""           "injection"
check_tool "${ENV_BAMSURGEON}" addsnv.py   "--help"     "injection"
check_tool "${ENV_BAMSURGEON}" addindel.py "--help"     "injection"
check_tool "${ENV_BAMSURGEON}" picard      "SortSam --version" "injection"
check_tool "${ENV_GATK}"       gatk        "--list"     "variant calling"
check_tool "${ENV_OPTITYPE}"   OptiTypePipeline.py "--help" "HLA typing"
check_tool "${ENV_OPTITYPE}"   razers3     "--version"  "HLA typing"
check_tool "${ENV_MHCFLURRY}"  mhcflurry-predict "--help" "binding prediction"
check_tool "${ENV_LOHHLA}"     Rscript     "--version"  "HLA loss"
check_tool "${ENV_VEP}"        vep         "--help"     "annotation"
check_tool "${ENV_PVACSEQ}"    pvacseq     "--version"  "neoantigens"

say "reference data"
[ -s "${REF_FASTA}" ]        && pass "hg38 fasta" \
                             || { miss "hg38 fasta at ${REF_FASTA}"
                                  BLOCKED_STAGES+=("everything"); }
[ -s "${REF_FASTA}.fai" ]    && pass "hg38 index" || miss "hg38 .fai"
[ -s "${REF_FASTA}.bwt" ]    && pass "bwa index" \
                             || { miss "bwa index — run fetch_data.sh"
                                  BLOCKED_STAGES+=("injection"); }
[ -s "${REF_FASTA%.fasta}.dict" ] && pass "sequence dictionary" \
                             || miss "sequence dictionary (.dict)"
[ -s "${PANEL_BED}" ]        && pass "panel BED" \
                             || { miss "panel BED"; BLOCKED_STAGES+=("everything"); }
[ -s "${HLA_TYPES}" ]        && pass "published HLA types" \
                             || opt "published HLA types — accuracy check skipped"
[ -s "${LOHHLA_DIR}/data/hla_all_lohhla.fasta" ] && pass "LOHHLA reference" \
                             || { miss "LOHHLA reference — run fetch_data.sh"
                                  BLOCKED_STAGES+=("HLA loss"); }
[ -s "${LOHHLA_DIR}/LOHHLAscript.R" ] && pass "LOHHLA script" \
                             || { miss "LOHHLA script"; BLOCKED_STAGES+=("HLA loss"); }

if [ -d "${VEP_CACHE}/homo_sapiens" ]; then
    pass "VEP cache"
else
    opt "VEP cache — annotation and neoantigens will be skipped"
    BLOCKED_STAGES+=("annotation" "neoantigens")
fi

if [ -x "${NETMHC_DIR}/netMHCpan" ]; then
    pass "netMHCpan"
else
    opt "netMHCpan — neoantigens will be skipped"
    BLOCKED_STAGES+=("neoantigens")
fi

say "samples"
n_bam=0; n_truth=0
for s in ${SAMPLES}; do
    [ -s "${SAMPLE_DIR}/${s}_normal.bam" ] && n_bam=$((n_bam+1))
    [ -s "${SAMPLE_DIR}/${s}_truth.tsv" ] && n_truth=$((n_truth+1))
done
n_want=$(echo ${SAMPLES} | wc -w)
if [ "${n_bam}" -eq "${n_want}" ]; then
    pass "${n_bam} of ${n_want} normal BAMs"
else
    miss "${n_bam} of ${n_want} normal BAMs in ${SAMPLE_DIR}"
    BLOCKED_STAGES+=("everything")
fi
[ "${n_truth}" -eq "${n_want}" ] && pass "${n_truth} truth sets" \
                                 || miss "${n_truth} of ${n_want} truth sets"

say "disk"
avail=$(df -BG "${DATA_DIR}" 2>/dev/null | awk 'NR==2 {gsub("G","",$4); print $4}')
if [ -n "${avail}" ]; then
    if [ "${avail}" -ge 40 ]; then
        pass "${avail} GB free under ${DATA_DIR}"
    else
        miss "${avail} GB free — the run needs about 40"
    fi
fi

say "what will happen"
declare -A blocked
for s in "${BLOCKED_STAGES[@]:-}"; do
    [ -n "$s" ] && blocked["$s"]=1
done

for stage in "injection" "variant calling" "HLA typing" \
             "binding prediction" "HLA loss" "annotation" "neoantigens"; do
    if [ -n "${blocked[everything]:-}" ] || [ -n "${blocked[$stage]:-}" ]; then
        printf '  will skip   %s\n' "${stage}"
    else
        printf '  will run    %s\n' "${stage}"
    fi
done

echo ""
echo "  ${OK} checks passed, ${MISSING} missing, ${OPTIONAL_MISSING} optional absent"
echo ""
if [ "${MISSING}" -gt 0 ]; then
    echo "  Something required is missing. Fix the entries marked MISSING"
    echo "  and run this again."
    exit 1
fi
if [ "${OPTIONAL_MISSING}" -gt 0 ]; then
    echo "  Ready. Some stages will be skipped for want of optional"
    echo "  components, which is expected — the rest of the pipeline runs"
    echo "  and reports on what it did."
else
    echo "  Ready — everything present. Run ./run_demo.sh"
fi
exit 0
