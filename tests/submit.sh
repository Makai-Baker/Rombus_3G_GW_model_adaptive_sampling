#!/bin/bash
#SBATCH --job-name=test_rombus_parallel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --output=test.log

ml mamba && conda activate rombus

python test_PhenomP.py