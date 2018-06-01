#!/bin/bash
#SBATCH -p holyseasgpu
#SBATCH -J m49_test_nomu
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 6-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 4000 #memory pool for cores
#SBATCH -o logs/%x.out
#SBATCH -e logs/%x.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

pcmd_integrate.py --config "setup_files/${SLURM_JOB_NAME}.py" --data ../data/m49_testannulus_z-gz.dat --results "results/${SLURM_JOB_NAME}_tau.csv"
RESULT=${PIPESTATUS[0]}
sacct -j $SLURM_JOB_ID --format=JOBID%20,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
exit $RESULT
