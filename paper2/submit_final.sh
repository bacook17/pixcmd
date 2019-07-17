#!/bin/bash

echo 'Models 20'
sbatch --array=[204,128,44] slurm/M87v2_m20.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m20.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m20.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m20.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m20.slurm
sleep 1

echo 'Models 21'
sbatch --array=[204,128,44] slurm/M87v2_m21.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m21.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m21.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m21.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m21.slurm
sleep 1

echo 'Models 22'
sbatch --array=[204,128,44] slurm/M87v2_m22.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m22.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m22.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m22.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m22.slurm
sleep 1

