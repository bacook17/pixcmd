#!/bin/bash
#SBATCH -p gpu_requeue
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 1-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 2000 #memory pool for cores
#SBATCH -o logs/mock_tau.out
#SBATCH -e logs/mock_tau.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python ../pcmdpy/pcmdpy/pcmd_integrate.py --config setup_files/mock_tau.py --results results/mock_tau.csv
RESULT=${PIPESTATUS[0]}
sacct -j "${SLURM_JOB_ID}".batch --format=JOBID%20,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
exit $RESULT
