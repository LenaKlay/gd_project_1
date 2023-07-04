#! /bin/bash

for ((K=3; K<=8; K++))
do
	for ((s=1; s<=9; s++))
	do
		qsub -q long.q -cwd -V -N simu_K_${K}_s_${s} -b y "python3 simu_K_${K}_s_${s}.py"
	done
done
