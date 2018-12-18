#!/bin/bash

echo "submitting Paper1 Mocks"
sbatch --array=[1-28] paper1_mocks.slurm

sleep 2

echo "submitting Paper1 Mismatches"
sbatch --array=[1-11] paper1_mismatches.slurm
