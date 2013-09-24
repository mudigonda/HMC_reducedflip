#!/bin/bash -l
source /usr/Moudles/init/bash

#SBATCH -p cortex
#SBATCH --time=12:30:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --ntasks-per-node=2

cd /global/home/users/mayur/HMC_reducedflip/
echo "Epsilon = ",$Epsilon
echo "Leap Size = ",$LeapSize
echo "Beta = ",$Beta
module load matlab 
matlab -nodisplay -nojvm -r "make_figures_fneval_cluster $LeapSize $Epsilon $Beta"
exit
