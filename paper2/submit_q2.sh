#!/bin/bash

echo 'M87'
sbatch --array=[201,125,41] slurm/M87_m7.slurm

sleep 1

echo 'M87v2'
sbatch --array=[201,125,41] slurm/M87v2_m7.slurm

sleep 1

echo 'M49'
sbatch --array=[201,121,37] slurm/M49_m7.slurm

sleep 1

echo 'NGC 3377'
sbatch --array=[174,98,42] slurm/NGC3377_m7.slurm

sleep 1

echo 'NGC 4993'
sbatch --array=[204,144,84] slurm/NGC4993_m7.slurm
