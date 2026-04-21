#!/bin/bash

mol=$1

cd $mol-asfe/$mol
pwd
for i in */
do
	sbatch $i/simulation.sh $i
done
