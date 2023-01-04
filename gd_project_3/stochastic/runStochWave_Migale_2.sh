#! /bin/bash

#
# This script submits jobs for the stochastic simulations.
# Is it meant to be run on the Migale cluster https://migale.inrae.fr
#

# Compile the simulation script
gcc stochwave_migale.c -lm -o stochwave

# Number of replicates for each set of parameters
export NREPS=1

# Create execution files
# Make the execution files executable
# and submit them
for s in 0.3000000 0.3189655 0.3379310 0.3568966 0.3758621 0.3948276 0.4137931 0.4327586 0.4517241 0.4706897 0.4896552 0.5086207 0.5275862 0.5465517 0.5655172 0.5844828 0.6034483 0.6224138 0.6413793 0.6603448 0.6793103 0.6982759 0.7172414 0.7362069 0.7551724 0.7741379 0.7931034 0.8120690 0.8310345 0.8500000
do
  for r in 0.4  0.8  1.2  1.6  2.0  2.4  2.8  3.2  3.6  4.0  4.4  4.8  5.2 5.6  6.0  6.4  6.8  7.2  7.6  8.0  8.4  8.8  9.2  9.6 10.0 10.4 10.8 11.2 11.6 12.0
  do
    for mig in 0.1
    do
      # Create execution file
      echo -e "#!/bin/bash\n./stochwave 1000 ${s} ${r} ${mig} 10000 ${NREPS} > data/stoch2_s-${s}_r-${r}_mig-${mig}.csv" > sim_s-${s}_r-${r}_mig-${mig}.sh
      # Make the file executable
      chmod +x sim_s-${s}_r-${r}_mig-${mig}.sh
      # Submit the file
      qsub -cwd -q short.q sim_s-${s}_r-${r}_mig-${mig}.sh
    done
  done
done
