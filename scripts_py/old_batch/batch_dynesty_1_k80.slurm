#!/bin/bash
#SBATCH --reservation=conroygpu
#SBATCH -p aagk80
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 0-02:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 1000 #memory pool for cores
#SBATCH -o logs/dynesty_1_k80.out
#SBATCH -e logs/dynesty_1_k80.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

python pcmd_integrate.py setup_files/dynesty_1_k80.py
sacct -j "${SLURM_JOB_ID}" --format=JOBID,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State
