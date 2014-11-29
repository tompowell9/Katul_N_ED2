#!/bin/bash
#SBATCH -J Kelsea_test
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -o Kelsea_test.%j
#
module load openmpi
cd $SLURM_SUBMIT_DIR
srun -n 1 ./ed2 -f ED2IN_Kelsea_test
