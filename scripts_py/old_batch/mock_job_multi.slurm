#!/bin/bash
#SBATCH -p holyseasgpu
#SBATCH -J mock_job
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 6-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 4000 #memory pool for cores
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

pcmd_integrate_twostage.py --config "setup_files/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.py" --results-init "results/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}_tau.csv" --results-final "results/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}_nonparam.csv"
RESULT=${PIPESTATUS[0]}
sacct -j $SLURM_JOB_ID --format=JOBID%20,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
exit $RESULT
