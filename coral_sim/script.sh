#!/bin/bash
#SBATCH --job-name=SC0
#SBATCH --ntasks-per-node=40
#SBATCH --nodelist=mrcd[105-106]
#SBATCH -N 2
#SBATCH --time=02-00:00:00
#SBATCH --partition=pi.skhatri3
#SBATCH --output output

mpirun ./main3d input3d
