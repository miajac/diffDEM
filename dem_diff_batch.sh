#!/bin/bash
#SBATCH --job-name=demDiff_batch
#SBATCH --nodes=3 # should be one node per DEM pair
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G # memory per node
#SBATCH --time=04:00:00 # kill time
#SBATCH --output=demDiff_%j.log # log file (%j = job ID)
#SBATCH --error=demDiff_%j.err

# Load required modules 
source ~/miniforge3/etc/profile.d/conda.sh
module load OpenMPI/4.1.4

# Activate conda environment
conda activate demDiff

# Run the batch differencing script
# Config file contains all DEM paths, pairs mode, and options
mpirun python dem_diff_batch.py config_batch_toy.yml