#!/bin/bash -l
#SBATCH -J ave_pei
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 00:05:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
sbatch -d afterany:$SLURM_JOB_ID debug-elm.sh 
srun -n 1 ipython pei_save.py 




 
