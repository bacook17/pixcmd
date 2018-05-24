#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 0-01:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 2000 #memory pool for cores
#SBATCH -o logs/datatest_2.out
#SBATCH -e logs/datatest_2.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_integrate.py setup_files/datatest_1.py
