#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 8 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:8 #Number of GPUs requested
#SBATCH --constraint=cuda-7.5 #require CUDA
#SBATCH -t 0-06:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 1000 #memory pool for cores
#SBATCH -o logs/batch_timing4.out
#SBATCH -e logs/batch_timing4.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_fitter.py setup_files/timing4.py
