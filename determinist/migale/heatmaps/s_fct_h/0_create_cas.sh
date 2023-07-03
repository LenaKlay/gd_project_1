#! /bin/bash

# Create one file.py per unit of precision. 
# Each file is dealing with one row of the heatmap.
echo "Precision ?  "
read precision

# Create directories

for cas in "a" "b_pos" "b_neg" "c" "d"
do
	mkdir cas_$cas
	mkdir cas_$cas/python
	mkdir cas_$cas/speed_coex

	cp 0_heatmap_numero.py 0_heatmap_numero_$cas.py
	sed -i "s/quel_cas/${cas}/g" 0_heatmap_numero_$cas.py

	for ((i=1; i<=precision; i++))
	do	
		cp 0_heatmap_numero_$cas.py heatmap_$i.py
        	sed -i "s/quelle_precision/${precision}/g" heatmap_$i.py
		sed -i "s/numero/${i}/g" heatmap_$i.py
		mv heatmap_$i.py cas_$cas/python
	done
	mv 0_heatmap_numero_$cas.py cas_$cas
done

exit 0

# Pour exécuter

# chmod +x create_heatmap_file.sh
# ./create_heatmap_file.sh
