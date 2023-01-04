This repository contains the codes to run and plot the numerical solutions (deterministic model) and the stochastic simulations presented in our article. 

# Contents

## Folder `deterministic/`

-  `heatmap.m` MATLAB / Octave code to generate the heatmaps. 
    The various numerical simulations presented in the manuscript are all produced by this Octave code or slight variations of it.
-  `evolution.m` plots the snapshots of Fig 1 and 4

## Folder `stochastic/`

-  `stochwave_migale.c` is the main C code of the simulation;  
-  `runStochWave_Migale_2.sh` is a code to run the code on a cluster;
-  `data/` contains the output of the simulations (the files were imported back from the cluster); 
-  `R/` contains the code to plot the outcome, and a subfolder `Pics/` where the figure(s) are saved.  


  
