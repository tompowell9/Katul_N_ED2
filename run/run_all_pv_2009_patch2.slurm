#!/bin/bash
#SBATCH -J PV_2009_patch2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -o PV_2009_patch2.%j
#
module load openmpi
cd $SLURM_SUBMIT_DIR
srun -n 1 ./ed2 -f ED2IN_all_PV_2009_patch2
