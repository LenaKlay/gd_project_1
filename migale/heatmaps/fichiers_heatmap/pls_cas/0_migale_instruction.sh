#! /bin/bash

echo "Precision ?  "
read precision

for cas in "a" "b_pos" "b_neg" "c" "d"
do
	for ((i=1; i<= precision; i++))
	do
        	qsub -q long.q -cwd -V -N heatmap_$i -b y "python3 cas_$cas/python/heatmap_$i.py"
	done
done
