#! /bin/bash

echo "nb_step ?"
read nb

for ((i=0; i<nb; i++))
do
	qsub -q long.q -cwd -V -N noeud_$i -b y "python3 file_$i.py"
done
