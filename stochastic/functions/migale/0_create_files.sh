#! /bin/bash

# Create one file.py per unit of precision. 
# Each file is dealing with one row of the heatmap).

#echo "Precision ?  "
#read precision

for ((K=3; K<=8; K++))
do
	for ((s=1; s<=9; s++))
	do
		cp 0_simu_migale.py simu_K_${K}_s_${s}.py
        	sed -i "s/K_val/${K}/g" simu_K_${K}_s_${s}.py
		sed -i "s/s_val/${s}/g" simu_K_${K}_s_${s}.py
	done
done

exit 0

# Pour exécuter

# chmod +x create_heatmap_file.sh
# ./create_heatmap_file.sh
