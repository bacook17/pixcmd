#!/bin/bash
#SBATCH -p gpu
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH --constraint=cuda-7.5 #require CUDA
#SBATCH -t 0-03:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 2000 #memory pool for cores
#SBATCH -o logs/timetest_GPU.out
#SBATCH -e logs/timetest_GPU.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python timetest.py
RESULT=$?
sacct -j "${SLURM_JOB_ID}" --format=JOBID,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State
exit $RESULT
