#!/bin/bash -l
#SBATCH -J trans-couple
#SBATCH -p debug
#SBATCH -o data-trans/output
#SBATCH -e data-trans/error
#SBATCH -N 16
#SBATCH -t 00:05:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
sbatch -d afterany:$SLURM_JOB_ID debug-pei.sh
srun -n 512 ./trans_er_couple -d data-trans    
