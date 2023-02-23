#! /bin/bash

echo "Precision ?  "
read precision

for ((i=45; i<= precision; i++))
do
        qsub -q long.q -cwd -V -N heatmap_$i -b y "python3 heatmap_$i.py"
done
