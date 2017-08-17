#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:2 #Number of GPUs requested
#SBATCH --constraint=cuda-7.5 #require CUDA
#SBATCH -t 0-00:01 #runtime in D-HH:MM
#SBATCH --mem 256 #memory pool for cores
#SBATCH -o logs/batch_gputest.out
#SBATCH -e logs/batch_gputest.err
#SBATCH --mail-type=END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python gpu_test.py
sjob "${SLURM_JOB_ID}"