#!/bin/bash -l

#SBATCH -p cortex
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu 1GB
#SBATCH --ntasks-per-node=1
#SBATCH -o outfile
#SBATCH -e errorfile

cd $HOME/HMC_reducedflip

module load matlab
ac_path="/clusterfs/cortex/scratch/mayur/HMC_reducedflip/2D/"
matlab -nodesktop -r "ac_summary_stats $ac_path" 
