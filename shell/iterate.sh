#!/bin/sh
#Leap Size
param1=11
param2=1.6
param3=.5
JobFile="param_test.sh"
count=0
#LeapSize
for LeapSize in 1 10 50 100 
do
#Epsilon
		for Epsilon in $(seq .7 .1 $param2);
		do
			#Beta
			for Beta in $(seq .01 .1 $param3);
			do
				echo $LeapSize,$Epsilon,$Beta
				Outfile=/global/home/users/mayur/Logs/HMCredflip_$LeapSize-$Epsilon-$Beta.o
				Errorfile=/global/home/users/mayur/Errors/HMCredflip_$LeapSize-$Epsilon-$Beta.e
				export LeapSize
				export Epsilon
				export Beta
				sbatch -o $Outfile -e $Errorfile -J HMC$LeapSize$Epsilon$Beta $JobFile
				count=$(($count+1))
			done
		done
done
echo $count jobs you maniacal fuck!
