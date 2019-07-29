#!/bin/bash

echo 'Models 28'
sbatch --array=[204,128,44] slurm/M87v2_m28.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m28.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m28.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m28.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m28.slurm
sleep 1

echo 'Models 29'
sbatch --array=[204,128,44] slurm/M87v2_m29.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m29.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m29.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m29.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m29.slurm
sleep 1

echo 'Models 30'
sbatch --array=[204,128,44] slurm/M87v2_m30.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m30.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m30.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m30.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m30.slurm
sleep 1

echo 'Models 31'
sbatch --array=[204,128,44] slurm/M87v2_m31.slurm
sleep 1
sbatch --array=[204,124,40] slurm/M49_m31.slurm
sleep 1
sbatch --array=[173,97,41] slurm/NGC3377_m31.slurm
sleep 1
sbatch --array=[16-17] slurm/NGC4993_m31.slurm
sleep 1
sbatch --array=[1-5] slurm/M31_m31.slurm
sleep 1

