#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH --constraint=cuda-7.5 #require CUDA
#SBATCH -t 2-00:00 #runtime in D-HH:MM
#SBATCH --mem 4000 #memory pool for cores
#SBATCH -o results/batch_SSP.out
#SBATCH -e results/batch_SSP.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python test_fitter.py --N_scale=1024 --N_walkers=20 --N_burn=500 --N_sample=5000 --SSP --require_gpu --require_cudac --append=_gpuSSP