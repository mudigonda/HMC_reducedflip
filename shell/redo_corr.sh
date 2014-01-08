#!/bin/bash -l
outfile=/clusterfs/cortex/scratch/mayur/Logs/redo.out
errorfile=/clusterfs/cortex/scratch/mayur/Errors/redo.err

#SBATCH -p cortex
#SBATCH --time=150:00:00
#SBATCH --mem-per-cpu 2GB
#SBATCH --ntasks-per-node=3
#SBATCH --constraint=nogpu
#SBATCH -o $outfile
#SBATCH -e $errorfile
#SBATCH --constraint=nogpu

cd $HOME/HMC_reducedflip

module load matlab
redo_path="/clusterfs/cortex/scratch/mayur/HMC_reducedflip/2D/"
matlab -nodesktop -r "redo_autocorr $redo_path" 
