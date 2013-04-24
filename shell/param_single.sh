#!/bin/bash

#PBS -q cortex
#PBS -l nodes=1:ppn=4:gpu
#PBS -l walltime=22:35:00
#PBS -o /global/home/users/mayur/Logs
#PBS -e /global/home/users/mayur/Errors
cd /global/home/users/mayur/HMC_reducedflip/
module load matlab 
matlab -nodisplay -nojvm -r "make_figures_fneval_cluster_wrapper" 
