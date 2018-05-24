#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 4 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:4 #Number of GPUs requested
#SBATCH -t 0-06:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 1000 #memory pool for cores
#SBATCH -o logs/batch_timing3.out
#SBATCH -e logs/batch_timing3.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_fitter.py setup_files/timing3.py
sacct -j "${SLURM_JOB_ID}" --format=JOBID,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
