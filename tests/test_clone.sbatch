#!/bin/bash
#SBATCH --job-name=test_clone
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=%x_%j.log
#
# Clone the repository the way someone else would, and fetch everything it
# says it needs.
#
# This is the check that matters. Everything up to now ran in a directory
# that had been assembled by hand, with symlinks into another workspace
# and files that were never committed. A clone has only what is on GitHub,
# and a fresh data directory has only what fetch_data.sh can download.
#
# What it cannot check: that conda installs, that netMHCpan can be
# licensed, that a reader understands the README. Those stay with the
# reader. What it does check is every path, every URL and every file the
# repository claims to carry.
#
# netMHCpan is copied in from the existing installation rather than
# downloaded, because it cannot be downloaded. Doing it here as a separate
# step mirrors what the README asks of the reader.
#
# Usage:
#   sbatch test_clone.sbatch                    clone and fetch
#   sbatch test_clone.sbatch --run              and then run the pipeline
set -uo pipefail

RUN_PIPELINE=0
[ "${1:-}" = "--run" ] && RUN_PIPELINE=1

REPO="https://github.com/selevycho/immune-escape-demo.git"

# compute nodes carry a minimal userland: no git, and ws_find is not on
# PATH. Both come from the login environment, so the paths are resolved at
# submission time and passed in.
WSV="${WSV:-}"
if [ -z "${WSV}" ]; then
    echo "ERROR: WSV not set — submit with"
    echo "  sbatch --export=ALL,WSV=\$(ws_find demo_verify) test_clone.sbatch"
    exit 1
fi

CLONE="${WSV}/repo"
DATA="${WSV}/data"

# --export=ALL carries the submitting shell's environment across, and
# config.sh honours anything already set. That is right for a normal run
# and wrong here: every path would point back at the workspace this was
# submitted from, and the test would pass without touching the clone.
# Only WSV and SRC_WS survive.
for v in DATA_DIR REF_DIR SAMPLE_DIR OUT_DIR REF_FASTA VEP_CACHE \
         HLA_TYPES GENCODE_GTF GENCODE_PROT PANEL_BED PANEL_TSV \
         LOHHLA_DIR NETMHC_DIR DEMO_ROOT; do
    unset "${v}"
done

echo "======================================================================"
echo " cloning and fetching, as a stranger would"
echo " workspace ${WSV}"
echo " started   $(date '+%F %T') on $(hostname)"
echo "======================================================================"

# ------------------------------------------------------------------ clone
# The clone happens on the login node: compute nodes here have no git.
# That is a property of this cluster, not of the pipeline, and it does not
# weaken the test — what matters is that the directory contains only what
# GitHub served, which it does.
echo ""
echo "=== clone ==="
if [ ! -d "${CLONE}" ]; then
    echo "  ERROR: no clone at ${CLONE}"
    echo "  clone it on a login node first:"
    echo "    git clone ${REPO} ${CLONE}"
    exit 1
fi
echo "  cloned from ${REPO}"
echo "  $(find "${CLONE}" -type f -not -path '*/.git/*' | wc -l) files"
echo "  $(du -sh --exclude=.git "${CLONE}" | cut -f1)"

# a file the pipeline needs but the repository forgot to commit shows up
# here rather than three stages into a run
echo ""
echo "=== what the clone carries ==="
for f in config.sh run_demo.sh README.md LICENSE \
         setup/install.sh setup/fetch_data.sh setup/check_deps.sh \
         setup/get_lohhla.sh \
         lib/report.py lib/make_neoantigens.py \
         panel/panel.bed steps/01a_addsnv.sh steps/08_report.sh \
         soft/lohhla/lohhla.patch; do
    printf "  %-32s %s\n" "${f}" \
        "$([ -s "${CLONE}/${f}" ] && echo present || echo MISSING)"
done

cd "${CLONE}"
export DATA_DIR="${DATA}"
mkdir -p "${DATA}"

# ----------------------------------------------------------------- lohhla
echo ""
echo "=== LOHHLA ==="
if ./setup/get_lohhla.sh 2>&1 | tail -20; then
    echo "  ok"
else
    echo "  FAILED — the pipeline can run without it, minus stage 5"
fi

# --------------------------------------------------------------- netmhcpan
echo ""
echo "=== netMHCpan ==="
SRC_WS="${SRC_WS:-}"
if [ -n "${SRC_WS}" ] && [ -d "${SRC_WS}/soft/netmhcpan/netMHCpan-4.1" ]; then
    cp -r "${SRC_WS}/soft/netmhcpan/netMHCpan-4.1" "${CLONE}/soft/"
    sed -i "s|setenv       NMHOME  .*|setenv       NMHOME  ${CLONE}/soft/netMHCpan-4.1|" \
        "${CLONE}/soft/netMHCpan-4.1/netMHCpan"
    if "${CLONE}/soft/netMHCpan-4.1/netMHCpan" -version >/dev/null 2>&1; then
        echo "  copied and NMHOME rewritten — this is what the README asks"
        echo "  a reader to do by hand after obtaining a licence"
    else
        echo "  copied but does not run"
    fi
else
    echo "  no local copy to place; stage 7 will be skipped"
fi

# ------------------------------------------------------------------ fetch
echo ""
echo "=== fetching reference data ==="
echo "  about 30 GB, plus an hour for the bwa index"
T0=$(date +%s)
if ./setup/fetch_data.sh > "${WSV}/fetch.log" 2>&1; then
    T1=$(date +%s)
    printf "  done in %dh%02dm\n" \
        $(( (T1-T0)/3600 )) $(( ((T1-T0)%3600)/60 ))
else
    echo "  fetch reported problems:"
    grep -E "ERROR|failed" "${WSV}/fetch.log" | head -10 | sed 's/^/    /'
fi
echo ""
grep -A20 "=== summary ===" "${WSV}/fetch.log" | tail -12

# ------------------------------------------------------------------ check
echo ""
echo "=== check_deps ==="
./setup/check_deps.sh 2>&1 | tail -22

if [ "${RUN_PIPELINE}" = "0" ]; then
    echo ""
    echo "======================================================================"
    echo " fetched. To run the pipeline against this clone:"
    echo ""
    echo "   cd ${CLONE}"
    echo "   export DATA_DIR=${DATA}"
    echo "   ./run_demo.sh"
    echo ""
    echo " or resubmit this with --run"
    echo "======================================================================"
    exit 0
fi

# -------------------------------------------------------------------- run
echo ""
echo "=== running the pipeline ==="
./run_demo.sh --local

echo ""
echo " finished $(date '+%F %T')"
