#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 2-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 1000 #memory pool for cores
#SBATCH -o logs/nestle_ssp.out
#SBATCH -e logs/nestle_ssp.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_integrate.py setup_files/nestle_ssp.py
sacct -j "${SLURM_JOB_ID}" --format=JOBID,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State
