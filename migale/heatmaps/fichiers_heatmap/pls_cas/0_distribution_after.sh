#! /bin/bash

# Move results to right directories

for cas in "a" "b_pos" "b_neg" "c" "d"
do
	for ((i=1; i<=precision; i++))
	do	
		mv ${cas}_speed_$i.txt speed_$i.txt
		mv speed_$i.txt cas_${cas}/speed_coex
		mv ${cas}_coex_$i.txt coex_$i.txt
		mv coex_$i.txt cas_${cas}/speed_coex
	done
done

exit 0

# Pour exécuter

# chmod +x create_heatmap_file.sh
# ./create_heatmap_file.sh
