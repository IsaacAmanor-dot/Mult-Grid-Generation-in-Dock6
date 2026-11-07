#!/bin/bash
#SBATCH --job-name=2ITY_multigrid
#SBATCH --partition=rn-long-40core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/multigrid_.out
#SBATCH --error=logs/multigrid_.err

set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"
python3 -u 001.prepare_multigrid.py generate
