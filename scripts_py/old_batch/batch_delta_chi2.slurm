#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 0-12:00 #runtime in D-HH:MM
#SBATCH --mem 4000 #memory pool for cores
#SBATCH -o results/batch_delta_chi2.out
#SBATCH -e results/batch_delta_chi2.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python loglike_variance.py
