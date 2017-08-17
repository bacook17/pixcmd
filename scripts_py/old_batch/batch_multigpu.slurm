#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:2 #Number of GPUs requested
#SBATCH --constraint=cuda-7.5 #require CUDA
#SBATCH -t 0-00:03 #runtime in D-HH:MM
#SBATCH --mem 1000 #memory pool for cores
#SBATCH -o logs/batch_multigpu.out
#SBATCH -e logs/batch_multigpu.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python test_fitter.py --N_scale=128 --N_walkers=20 --N_burn=0 --N_sample=2 --fixed_seed --SSP --sample_ball --require_gpu --require_cudac --append=_multigpu --N_threads=2
