#!/bin/bash
#SBATCH -p holyseasgpu
#SBATCH -J paper1_mock
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 #ensure that all cores are on one machine
#SBATCH --gres=gpu:1 #Number of GPUs requested
#SBATCH -t 14-00:00 #runtime in D-HH:MM
#SBATCH --mem-per-cpu 8000 #memory pool for cores
#SBATCH -o logs/%x_%a.out
#SBATCH -e logs/%x_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL #alert when done
#SBATCH --mail-user=bcook@cfa.harvard.edu #Email to send to

MPROF_FILE="logs/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.mem"
command rm $MPROF_FILE

CONFIG_FILE="setup_files/mocks_paper1/mock_${SLURM_ARRAY_TASK_ID}.py"
RESULTS_FILE="results/${SLURM_JOB_NAME}_${SLURM_ARRAY_TASK_ID}.csv"

COMMAND="pcmd_integrate --config ${CONFIG_FILE} --results ${RESULTS_FILE}"

mprof run --nopython -T 10 -o $MPROF_FILE $COMMAND 

RESULT=${PIPESTATUS[0]}
sacct -j $SLURM_JOB_ID --format=JOBID%20,JobName,NTasks,AllocCPUs,AllocGRES,Partition,Elapsed,MaxRSS,MaxVMSize,MaxDiskRead,MaxDiskWrite,State 
exit $RESULT
