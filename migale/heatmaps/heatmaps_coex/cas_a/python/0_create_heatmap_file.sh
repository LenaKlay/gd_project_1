#! /bin/bash

# Create one file.py per unit of precision. 
# Each file is dealing with one row of the heatmap).

echo "Precision ?  "
read precision

for ((i=1; i<=precision; i++))
do
	cp 0_heatmap_numero.py heatmap_$i.py
        sed -i "s/quelle_precision/${precision}/g" heatmap_$i.py
	sed -i "s/numero/${i}/g" heatmap_$i.py
done

exit 0

# Pour exécuter

# chmod +x create_heatmap_file.sh
# ./create_heatmap_file.sh
