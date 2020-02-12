#!/bin/bash -l
#SBATCH -J OMFITTrans
#SBATCH -p debug
#SBATCH -o BOUT.out
#SBATCH -e BOUT.err
#SBATCH -D .
#SBATCH -C haswell
#SBATCH -N 8
#SBATCH -t 00:30:00
#SBATCH -c 1

module swap PrgEnv-intel PrgEnv-gnu
module swap gcc gcc/4.9.3

chmod +x ./bout.exe
srun -n 256 ./bout.exe restart append