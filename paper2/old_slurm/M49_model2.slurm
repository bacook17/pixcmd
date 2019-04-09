#!/bin/bash
#SBATCH -p fas_gpu,holyseasgpu,gpu
#SBATCH -J M49_m2
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 7-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 4000 #memory pool for cores
#SBATCH -o logs/%x_r%a.out
#SBATCH -e logs/%x_r%a.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

MPROF_FILE="logs/${SLURM_JOB_NAME}_r${SLURM_ARRAY_TASK_ID}.mem"
command rm $MPROF_FILE

CONFIG_FILE="setup_files/M49/M49_model2.py"
DATA_FILE="../data/M49/pcmds/M49_z_gz_${SLURM_ARRAY_TASK_ID}.pcmd"
RESULTS_FILE="results/${SLURM_JOB_NAME}_r${SLURM_ARRAY_TASK_ID}.csv"

COMMAND="pcmd_integrate --config ${CONFIG_FILE} --data ${DATA_FILE} --results ${RESULTS_FILE}"

mprof run --nopython -T 10 -o $MPROF_FILE $COMMAND 

RESULT=${PIPESTATUS[0]}
sacct -j $SLURM_JOB_ID --format=JOBID%20,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
exit $RESULT
