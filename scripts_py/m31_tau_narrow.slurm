#!/bin/bash
#SBATCH -p holyseasgpu
#SBATCH -J m31_tau_narrow
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 7-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 4000 #memory pool for cores
#SBATCH -o logs/%x_r%a.out
#SBATCH -e logs/%x_r%a.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

pcmd_integrate.py --config "setup_files/m31_taufull_narrow.py" --data "../data/m31/m31_r${SLURM_ARRAY_TASK_ID}.pcmd" --results "results/${SLURM_JOB_NAME}_r${SLURM_ARRAY_TASK_ID}.csv"
RESULT=${PIPESTATUS[0]}
sacct -j $SLURM_JOB_ID --format=JOBID%20,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
exit $RESULT
