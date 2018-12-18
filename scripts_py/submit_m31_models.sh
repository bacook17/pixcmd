#!/bin/bash

echo "submitting Summer Model 1"
sbatch --array=[1-9] m31_summer_model1.slurm 

sleep 2

echo "submitting Winter Model 1"
sbatch --array=[1-10] m31_winter_model1.slurm

sleep 2

echo "submitting Summer Model 2"
sbatch --array=[1-9] m31_summer_model2.slurm 

sleep 2

echo "submitting Winter Model 2"
sbatch --array=[1-10] m31_winter_model2.slurm

sleep 2

#echo "submitting Summer Model 3"
#sbatch --array=[1-9] m31_summer_model3.slurm 

#sleep 2

#echo "submitting Winter Model 3"
#sbatch --array=[1-10] m31_winter_model3.slurm

#sleep 2

echo "submitting Summer Model 4"
sbatch --array=[1,2,3,8,9] m31_summer_model4.slurm 

sleep 2

echo "submitting Winter Model 4"
sbatch --array=[1,2,3,7,8,9,10] m31_winter_model4.slurm
