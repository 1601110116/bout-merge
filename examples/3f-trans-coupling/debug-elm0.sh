#!/bin/bash -l
#SBATCH -J elm-couple
#SBATCH -p debug
#SBATCH -o data-elm/output
#SBATCH -e data-elm/error
#SBATCH -N 16
#SBATCH -t 00:05:00

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

#run the application:
sbatch -d afterany:$SLURM_JOB_ID debug-trans0.sh
srun -n 512 ./elm_pb_couple -d data-elm   
