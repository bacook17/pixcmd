#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 2 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:2 #Number of GPUs requested
#SBATCH -t 7-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 1000 #memory pool for cores
#SBATCH -o logs/batch_ball_large.out
#SBATCH -e logs/batch_ball_large.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_fitter.py setup_files/ball_large.py
sacct -j "${SLURM_JOB_ID}" --format=JOBID,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State
