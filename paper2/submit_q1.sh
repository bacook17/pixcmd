#!/bin/bash

echo 'M87'
sbatch --array=[204,128,44] slurm/M87_m7.slurm

sleep 1

echo 'M87v2'
sbatch --array=[204,128,44] slurm/M87v2_m7.slurm

sleep 1

echo 'M49'
sbatch --array=[204,124,40] slurm/M49_m7.slurm

sleep 1

echo 'NGC 3377'
sbatch --array=[173,97,41] slurm/NGC3377_m7.slurm

sleep 1

echo 'NGC 4993'
sbatch --array=[203,143,83] slurm/NGC4993_m7.slurm

sleep 1

echo 'M31'
sbatch --array=[1-5] slurm/M31_m7.slurm

sleep 1

echo 'M51'
sbatch --array=[1-5] slurm/M51_m7.slurm
