#!/bin/bash
#SBATCH -p debug
#SBATCH -N 1
#SBATCH -t 1:00
#SBATCH -J 0_data
#SBATCH -L SCRATCH

#OpenMP settings:
export OMP_NUM_THREADS=2
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n 4 -c 2  ./qparwcoll_bndry -d data
