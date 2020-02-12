#!/bin/bash -l
#SBATCH -J OMFIT6f
#SBATCH -p debug
#SBATCH -o BOUT.out
#SBATCH -e BOUT.err
#SBATCH -D .
#SBATCH -C haswell
#SBATCH -N 32
#SBATCH -t 00:30:00
#SBATCH -c 1

module swap PrgEnv-intel PrgEnv-gnu
module swap gcc gcc/4.9.3

chmod +x ./bout.exe
srun -n 1024 ./bout.exe restart append