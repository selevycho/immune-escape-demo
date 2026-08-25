#!/bin/bash
#SBATCH --job-name=demo_step
#SBATCH --partition=cpu-single
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.log
#
# BAMSurgeon runs one bwa process per thread and each loads the whole hg38
# index, about 6 GB. Two threads inside 16 GB is the combination that held
# across forty samples; four threads on a login node get killed partway
# and bwa returns an empty SAM, which addsnv reports as a site that would
# not realign.
set -uo pipefail
export DATA_DIR="${DATA_DIR:-$(ws_find immune_demo)}"
export THREADS=2
exec "$@"
