#!/bin/sh
#Leap Size
param1=11
param2=1.3
param3=.04
JobFile="param_test.sh"
#LeapSize
for LeapSize in 14 16
do
#Epsilon
		for Epsilon in $(seq .8 .1 $param2);
		do
			#Beta
			for Beta in $(seq .01 .01 $param3);
			do
				echo $LeapSize,$Epsilon,$Beta
				Outfile=/global/home/users/mayur/Logs/HMCredflip_$LeapSize-$Epsilon-$Beta.o
				Errorfile=/global/home/users/mayur/Errors/HMCredflip_$LeapSize-$Epsilon-$Beta.e
				export LeapSize
				export Epsilon
				export Beta
				sbatch -o $Outfile -e $Errorfile -J HMCredflip $JobFile
			done
		done
done
