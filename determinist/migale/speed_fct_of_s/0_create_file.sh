#! /bin/bash

# Create one file.py per unit of precision. 
# Each file is dealing with one row of the heatmap).

echo "nb_step ?"
read nb

for ((i=0; i<nb; i++))
do
	cp 0_speed_fct_of_s.py file_$i.py
	sed -i "s/numero/${i}/g" file_$i.py
done

exit 0

# Pour exécuter

# chmod +x create_heatmap_file.sh
# ./create_heatmap_file.sh
