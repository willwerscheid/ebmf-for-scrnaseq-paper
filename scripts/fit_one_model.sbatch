#!/bin/bash

#SBATCH --partition=broadwl
#SBATCH --account=pi-mstephens
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --exclusive
#SBATCH --time=36:00:00

# This script allocates computing resources (CPUs, memory), loads R,
# and runs fit_one_model.R.

# Get the command-line arguments.
DATFILE=${1}
METHOD=${2}
K=${3}
SELGENES=${4}
OUTFILE=${5}
echo datfile=${DATFILE}
echo method=${METHOD}
echo k=${K}
echo select.genes=${SELGENES}
echo outfile=${OUTFILE}

# For reproducibility, I add version numbers to the "module load" calls.
module load R/3.6.1

# Run the R script.
echo "Running fit_one_model.R"
export OPENBLAS_NUM_THREADS=1
./fit_one_model.R --datfile ../data/${DATFILE} --method ${METHOD} -k ${K} --selgenes ${SELGENES} \
  -o ../output/${OUTFILE}.rds > ../output/${OUTFILE}.Rout
