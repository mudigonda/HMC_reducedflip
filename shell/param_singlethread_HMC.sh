#!/bin/bash -l
source /usr/Modules/init/bash

#SBATCH --job-name=HMC_oneparam
#SBATCH -p cortex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=03:30:00
#SBATCH -o /global/home/users/mayur/Logs/hmc_1pt2_15_pt03.o
#SBATCH -e /global/home/users/mayur/Errors/hmc_1pt2_15_pt03.e
cd /global/home/users/mayur/HMC_reducedflip/
module load matlab 
Epsilon=1.2
LeapSize=15
Beta=0.03
echo "Epsilon = ",$Epsilon
echo "Leap Size = ",$LeapSize
echo "Beta = ",$Beta
matlab -nodisplay -nojvm -r "make_figures_fneval_cluster $LeapSize $Epsilon $Beta"
