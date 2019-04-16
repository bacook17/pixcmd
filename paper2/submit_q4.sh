#!/bin/bash

echo 'M87'
sbatch --array=[203,127,43] slurm/M87_m7.slurm

sleep 1

echo 'M87v2'
sbatch --array=[203,127,43] slurm/M87v2_m7.slurm

sleep 1

echo 'M49'
sbatch --array=[203,123,39] slurm/M49_m7.slurm

sleep 1

echo 'NGC 3377'
sbatch --array=[176,100,44] slurm/NGC3377_m7.slurm

sleep 1

echo 'NGC 4993'
sbatch --array=[202,142,82] slurm/NGC4993_m7.slurm
