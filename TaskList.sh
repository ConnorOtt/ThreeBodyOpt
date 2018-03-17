#!/bin/bash

# This bash script runs the threeBody C script for a number of
# conditions.  

# Compile three body
gcc -o threeBody threeBody.c -lm

# set for loops for conditions

obj_array=( 1, 2 )
Clearance_array=( 0, 10, 100, 1000, 5000, 10000, 50000, 100000 )
for obj in "${obj_array[@]}"
do
	for Cle in "${Clearance_array[@]}"
	do
		echo "Omptimizing Objective $obj with clearance $Cle m and accuracy 0.5"
		# call threebody with the inputs
		./threeBody $obj $Cle 0.5
	done
done
