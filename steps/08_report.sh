#!/bin/bash
#
# Step 8 — the report.
#
# Runs on a login node in seconds; everything it needs is already on disk.
#
# Usage:
#   ./08_report.sh
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${HERE}/../config.sh"

use_env "${ENV_GATK}" || exit 1
python "${DEMO_ROOT}/lib/report.py" "${OUT_DIR}" "${SAMPLE_DIR}" ${SAMPLES}
