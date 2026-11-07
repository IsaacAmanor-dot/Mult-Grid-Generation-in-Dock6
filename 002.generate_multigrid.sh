#!/bin/bash
#SBATCH --job-name=2ITY_multigrid
#SBATCH --partition=rn-long-40core
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --chdir=/gpfs/projects/rizzo/iamanor/DOCK6_Development/Similarity_Island/Multi-Grid_testing
#SBATCH --output=logs/multigrid_.out
#SBATCH --error=logs/multigrid_.err

set -euo pipefail

python3 -u 001.prepare_multigrid.py generate
