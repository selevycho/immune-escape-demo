#!/bin/bash
#
# Every path the demo needs, in one place.
#
# The scripts under steps/ read these variables and contain no absolute
# paths of their own. Editing this file is the whole of the configuration.
#
# Data lives outside the repository. The reference genome alone is 3 GB
# and the VEP cache 26 GB, so DATA_DIR should point at a filesystem with
# room — a scratch or work directory rather than a home directory.

DEMO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ------------------------------------------------------------------ data
# everything downloaded or produced goes here; nothing under DEMO_ROOT is
# written to except logs
DATA_DIR="${DATA_DIR:-${HOME}/immune_demo_data}"

REF_DIR="${REF_DIR:-${DATA_DIR}/reference}"
SAMPLE_DIR="${SAMPLE_DIR:-${DATA_DIR}/samples}"
OUT_DIR="${OUT_DIR:-${DATA_DIR}/results}"

REF_FASTA="${REF_FASTA:-${REF_DIR}/Homo_sapiens_assembly38.fasta}"
VEP_CACHE="${VEP_CACHE:-${REF_DIR}/vep_cache}"
HLA_TYPES="${HLA_TYPES:-${REF_DIR}/1000G_HLA_types.txt}"
GENCODE_GTF="${GENCODE_GTF:-${REF_DIR}/gencode.v46.annotation.gtf.gz}"
GENCODE_PROT="${GENCODE_PROT:-${REF_DIR}/gencode.v46.pc_translations.fa.gz}"

# ------------------------------------------------------- in the repository
PANEL_BED="${PANEL_BED:-${DEMO_ROOT}/panel/panel.bed}"
PANEL_TSV="${PANEL_TSV:-${DEMO_ROOT}/panel/panel_annotated.tsv}"
LOHHLA_DIR="${LOHHLA_DIR:-${DEMO_ROOT}/soft/lohhla}"

# ---------------------------------------------------------------- licensed
# request from https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/
# and unpack here; the pVACseq step is skipped with a message when absent
NETMHC_DIR="${NETMHC_DIR:-${DEMO_ROOT}/soft/netMHCpan-4.1}"

# ------------------------------------------------------------------ zenodo
ZENODO_DOI="10.5281/zenodo.22086838"
ZENODO_RECORD="22086838"

# ------------------------------------------------------------ environments
ENV_BAMSURGEON="${ENV_BAMSURGEON:-bamsurgeon_env}"
ENV_GATK="${ENV_GATK:-bio_work}"
ENV_OPTITYPE="${ENV_OPTITYPE:-optitype_stable}"
ENV_MHCFLURRY="${ENV_MHCFLURRY:-mhc_env}"
ENV_LOHHLA="${ENV_LOHHLA:-lohhla_env}"
ENV_VEP="${ENV_VEP:-vep_env}"
ENV_PVACSEQ="${ENV_PVACSEQ:-pvac_env}"

if [ -z "${CONDA_SH:-}" ]; then
    for c in "$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh" \
             "${HOME}/miniconda3/etc/profile.d/conda.sh" \
             "${HOME}/anaconda3/etc/profile.d/conda.sh" \
             "/opt/conda/etc/profile.d/conda.sh"; do
        [ -f "$c" ] && CONDA_SH="$c" && break
    done
fi

# --------------------------------------------------------------- samples
# chosen so that every stage produces something: each carries indels as
# well as substitutions, all three HLA loci type successfully, and two of
# the five yield a LOHHLA call
SAMPLES="${SAMPLES:-B003 B007 B012 B020 O011}"

THREADS="${THREADS:-4}"
MEM_GB="${MEM_GB:-12}"

# --------------------------------------------------------------- helpers
use_env() {
    local env_name="$1"
    if [ -z "${CONDA_SH:-}" ] || [ ! -f "${CONDA_SH}" ]; then
        echo "ERROR: conda.sh not found; set CONDA_SH in config.sh" >&2
        return 1
    fi
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    conda activate "${env_name}" 2>/dev/null || {
        echo "ERROR: environment '${env_name}' does not exist" >&2
        echo "       run setup/install.sh first" >&2
        return 1
    }
}

say() { printf '\n=== %s ===\n' "$*"; }
note() { printf '  %s\n' "$*"; }
fail() { printf '  ERROR: %s\n' "$*" >&2; }

# functions are not inherited by a child shell the way variables are,
# so they have to be exported explicitly
export -f use_env say note fail step_header 2>/dev/null || true

export DEMO_ROOT DATA_DIR REF_DIR SAMPLE_DIR OUT_DIR
export REF_FASTA VEP_CACHE HLA_TYPES PANEL_BED PANEL_TSV
export GENCODE_GTF GENCODE_PROT
export LOHHLA_DIR NETMHC_DIR ZENODO_DOI ZENODO_RECORD
export ENV_BAMSURGEON ENV_GATK ENV_OPTITYPE ENV_MHCFLURRY
export ENV_LOHHLA ENV_VEP ENV_PVACSEQ CONDA_SH SAMPLES THREADS MEM_GB

step_header() {
    echo ""
    echo "======================================================================"
    echo " $*"
    echo "======================================================================"
}
export -f step_header
