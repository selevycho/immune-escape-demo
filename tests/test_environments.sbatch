#!/bin/bash
#SBATCH --job-name=test_envs
#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
#
# Do the environment files still build?
#
# The yml files pin exact versions, which is what makes the results
# reproducible and also what will eventually break: conda channels remove
# old builds after a year or two. This rebuilds every one of them from
# scratch and reports which resolve and which do not.
#
# They go to a separate prefix rather than the usual envs directory, so a
# failed build cannot touch the environments the pipeline actually uses.
# Nothing here modifies anything the demo depends on.
#
# Usage:
#   sbatch test_environments.sbatch
#   sbatch test_environments.sbatch optitype_stable    just one
set -uo pipefail

# SLURM copies the submitted script into its own spool directory, so
# BASH_SOURCE points at /var/spool/slurm rather than at the repository.
# The path is passed in at submission time instead.
DEMO="${DEMO_ROOT:-${HOME}/immune_escape_demo}"
HERE="${DEMO}/tests"
YML_DIR="${DEMO}/setup/environments"
PREFIX="${TESTENV:-${HOME}/miniconda3/envs_test}"
WORK="$(mktemp -d)"

ENVS="${*:-bamsurgeon_env bio_work optitype_stable mhc_env lohhla_env vep_env pvac_env}"

echo "======================================================================"
echo " rebuilding environments from scratch"
echo " into ${PREFIX}"
echo " started $(date '+%F %T') on $(hostname)"
echo "======================================================================"
echo ""
echo "  This does not touch the environments the pipeline uses. It answers"
echo "  one question: would someone cloning this repository today be able"
echo "  to build what it asks for."
echo ""

source "$(conda info --base)/etc/profile.d/conda.sh"

SOLVER="conda"
command -v mamba >/dev/null 2>&1 && SOLVER="mamba"
echo "  solver: ${SOLVER}"
echo ""

mkdir -p "${PREFIX}"

OK=0; FAILED=0
declare -a BROKEN

for E in ${ENVS}; do
    SRC="${YML_DIR}/${E}.yml"
    echo "----------------------------------------------------------------------"
    echo " ${E}"

    if [ ! -s "${SRC}" ]; then
        echo "   no ${SRC}"
        FAILED=$((FAILED+1)); BROKEN+=("${E}")
        continue
    fi

    # name and prefix are stripped so the environment lands where this
    # script decides rather than where it was exported from
    grep -v "^name:" "${SRC}" | grep -v "^prefix:" > "${WORK}/${E}.yml"

    n_pkg=$(grep -c "^  - " "${WORK}/${E}.yml")
    echo "   ${n_pkg} packages requested"

    if [ -d "${PREFIX}/${E}" ]; then
        echo "   already built here — removing to test a clean build"
        rm -rf "${PREFIX}/${E}"
    fi

    T0=$(date +%s)
    if ${SOLVER} env create -f "${WORK}/${E}.yml" -p "${PREFIX}/${E}" \
            > "${WORK}/${E}.log" 2>&1; then
        T1=$(date +%s)
        n=$(ls "${PREFIX}/${E}/conda-meta"/*.json 2>/dev/null | wc -l)
        printf "   built in %dm%02ds, %s packages installed\n" \
            $(( (T1-T0)/60 )) $(( (T1-T0)%60 )) "${n}"
        OK=$((OK+1))
    else
        echo "   FAILED:"
        grep -iE "PackagesNotFound|UnsatisfiableError|nothing provides|conflict" \
            "${WORK}/${E}.log" | head -6 | sed 's/^/     /'
        [ -s "${WORK}/${E}.log" ] && tail -4 "${WORK}/${E}.log" | sed 's/^/     /'
        cp "${WORK}/${E}.log" "${HERE}/failed_${E}.log"
        echo "     full log: ${HERE}/failed_${E}.log"
        FAILED=$((FAILED+1)); BROKEN+=("${E}")
    fi
done

echo ""
echo "======================================================================"
echo " do the tools inside actually run"
echo "======================================================================"
echo ""

# a package that installs is not necessarily a package that works, so each
# tool is invoked; exit codes are unreliable here (bwa returns 1 with no
# arguments, gatk writes its banner to stderr) so any output at all counts
probe() {
    local env_name="$1" tool="$2" args="${3:-}"
    local bin="${PREFIX}/${env_name}/bin/${tool}"
    printf "  %-22s %-22s " "${env_name}" "${tool}"
    if [ ! -x "${bin}" ]; then
        echo "not installed"
        return 1
    fi
    if [ -n "$("${bin}" ${args} 2>&1 | head -3)" ]; then
        echo "ok"
        return 0
    fi
    echo "installed but silent"
    return 1
}

TOOLS_OK=0; TOOLS_BAD=0
run_probe() { probe "$@" && TOOLS_OK=$((TOOLS_OK+1)) || TOOLS_BAD=$((TOOLS_BAD+1)); }

run_probe bamsurgeon_env  samtools            "--version"
run_probe bamsurgeon_env  bwa                 ""
run_probe bamsurgeon_env  addsnv.py           "--help"
run_probe bamsurgeon_env  addindel.py         "--help"
run_probe bamsurgeon_env  picard              "SortSam --version"
run_probe bio_work        gatk                "--list"
run_probe optitype_stable OptiTypePipeline.py "--help"
run_probe optitype_stable razers3             "--version"
run_probe mhc_env         mhcflurry-predict   "--help"
run_probe lohhla_env      Rscript             "--version"
run_probe vep_env         vep                 "--help"
run_probe pvac_env        pvacseq             "--version"

echo ""
echo "======================================================================"
echo " result"
echo "======================================================================"
echo ""
echo "  environments   ${OK} built, ${FAILED} failed"
echo "  tools          ${TOOLS_OK} run, ${TOOLS_BAD} did not"
echo ""

if [ "${FAILED}" -gt 0 ]; then
    echo "  unresolved: ${BROKEN[*]}"
    echo ""
    echo "  When a pinned build disappears from a channel the fix is to"
    echo "  relax the pins and accept current versions of the same"
    echo "  packages:"
    echo ""
    echo "    sed -E 's/=[0-9][^=]*\$//' setup/environments/NAME.yml > relaxed.yml"
    echo "    conda env create -f relaxed.yml -n NAME"
    echo ""
    echo "  The pipeline will run. The numbers may shift a little, since"
    echo "  that is no longer the exact stack these results came from."
fi

echo ""
echo "  built environments are at ${PREFIX}"
echo "  remove them with:  rm -rf ${PREFIX}"
echo ""
echo " finished $(date '+%F %T')"

rm -rf "${WORK}"
[ "${FAILED}" -gt 0 ] && exit 1
exit 0
