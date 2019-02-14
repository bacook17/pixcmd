#!/bin/bash
#SBATCH -p fas_gpu,holyseasgpu,gpu
#SBATCH -J postprocess_mocks
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 0-12:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 8000 #memory pool for cores
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python postprocess_mocks.py