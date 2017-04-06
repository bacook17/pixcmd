#!/bin/bash
#SBATCH -p conroy
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH -t 0-06:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 4000 #memory pool for cores
#SBATCH -o results/batch_short_noGPU.out
#SBATCH -e results/batch_short_noGPU.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python test_fitter.py --N_threads=1 --N_scale=256 --N_walkers=20 --N_burn=0 --N_sample=100 --no_gpu --SSP --append=_short_noGPU
