echo "Submitting M87 Model 3, Region 44"
sbatch --array=[44] M87_model3.slurm
sleep 1

echo "Submitting M87 Model 4, Region 101-104"
sbatch --array=[101-104] M87_model4.slurm
sleep 1

echo "Submitting M87 Model 5, Region 204"
sbatch --array=[204] M87_model5.slurm
sleep 1

echo "Submitting M87 Model 6, Region 264"
sbatch --array=[264] M87_model6.slurm
sleep 1

echo "Submitting M87 Model 7, Region 104"
sbatch --array=[104] M87_model7.slurm
sleep 1

echo "Submitting M87 Model 8, Region 104"
sbatch --array=[104] M87_model8.slurm
sleep 1

echo "Submitting M87 Model 9, Region 104"
sbatch --array=[104] M87_model9.slurm
sleep 1

echo "Submitting M49 Model 3, Region 40"
sbatch --array=[40] M49_model3.slurm
sleep 1

echo "Submitting M49 Model 4, Region 97-100"
sbatch --array=[97-100] M49_model4.slurm
sleep 1

echo "Submitting M49 Model 5, Region 204"
sbatch --array=[204] M49_model5.slurm
sleep 1

echo "Submitting M49 Model 6, Region 256"
sbatch --array=[256] M49_model6.slurm
sleep 1

echo "Submitting M49 Model 7, Region 100"
sbatch --array=[100] M49_model7.slurm
sleep 1

echo "Submitting M49 Model 8, Region 100"
sbatch --array=[100] M49_model8.slurm
sleep 1

echo "Submitting M49 Model 9, Region 100"
sbatch --array=[100] M49_model9.slurm
sleep 1

echo "Submitting NGC3377 Model 3, Region 41"
sbatch --array=[41] NGC3377_model3.slurm
sleep 1

echo "Submitting NGC3377 Model 4, Region 97-100"
sbatch --array=[97-100] NGC3377_model4.slurm
sleep 1

echo "Submitting NGC3377 Model 5, Region 173"
sbatch --array=[173] NGC3377_model5.slurm
sleep 1

echo "Submitting NGC3377 Model 6, Region 241"
sbatch --array=[241] NGC3377_model6.slurm
sleep 1

echo "Submitting NGC3377 Model 7, Region 97"
sbatch --array=[97] NGC3377_model7.slurm
sleep 1

echo "Submitting NGC3377 Model 8, Region 97"
sbatch --array=[97] NGC3377_model8.slurm
sleep 1

echo "Submitting NGC3377 Model 9, Region 97"
sbatch --array=[97] NGC3377_model9.slurm
sleep 1

echo "Submitting NGC4993 Model 3, Region 35"
sbatch --array=[35] NGC4993_model3.slurm
sleep 1

echo "Submitting NGC4993 Model 4, Region 81-84"
sbatch --array=[81-84] NGC4993_model4.slurm
sleep 1

echo "Submitting NGC4993 Model 5, Region 103"
sbatch --array=[103] NGC4993_model5.slurm
sleep 1

#echo "Submitting NGC4993 Model 6, Region 256"
#sbatch --array=[256] NGC4993_model6.slurm
#sleep 1

echo "Submitting NGC4993 Model 7, Region 83"
sbatch --array=[83] NGC4993_model7.slurm
sleep 1

echo "Submitting NGC4993 Model 8, Region 83"
sbatch --array=[83] NGC4993_model8.slurm
sleep 1

echo "Submitting NGC4993 Model 9, Region 83"
sbatch --array=[83] NGC4993_model9.slurm
sleep 1

