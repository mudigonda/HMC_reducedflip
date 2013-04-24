#!/bin/sh
#Leap Size
param1=11
param2=1.3
param3=.05 #Mayur put this to limit jobs <12
#LeapSize
for i in 1 5 10 15 20
do
#Epsilon
		for j in $(seq .8 .1 $param2);
		do
			#Beta
			for k in $(seq .01 .01 $param3);
			do
				echo $i,$j,$k
				qsub param_test.sh  -v "LeapSize=$i,Epsilon=$j,Beta=$k"
			done
		done
done
