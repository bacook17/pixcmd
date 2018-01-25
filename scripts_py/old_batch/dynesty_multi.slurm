#!/bin/bash
#SBATCH -p gpu_requeue
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH --constraint=cuda-7.5 #require CUDA
#SBATCH -t 6-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 4000 #memory pool for cores
#SBATCH -o logs/dynesty_%a.out
#SBATCH -e logs/dynesty_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_integrate.py setup_files/dynesty_${SLURM_ARRAY_TASK_ID}.py
