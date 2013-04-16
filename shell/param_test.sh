#!/bin/bash

#PBS -q cortex
#PBS -l nodes=1:ppn=10
#PBS -o /global/home/users/mayur/Logs
#PBS -e /global/home/users/mayur/Errors
cd /global/home/users/mayur/HMC_reducedflip/
module load matlab 
