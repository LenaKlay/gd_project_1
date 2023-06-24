#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies

import numpy as np
import matplotlib.pyplot as plt
import os


### Seed for reproductibility
#np.random.seed(28)
    
### Parameters

ini = 12                    # Carrying capacity for a spatial interval of size 1
T = 100                     # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)
dt = 0.001                  # Time interval
r = 0                       # Intrasic growth rate
s = 0.4                     # Disadvantage for drive
nb_i = 1000                 # Number of runs


extinction = np.ones((2,nb_i))*(-1)

# Loop on the alleles (0 for wt, 1 for drive)
for allele in range(2): 
    
    # Fecondity
    f = [0,1-s][allele]
    
    # Loop on the runs    
    for i in range(nb_i) :  
    
        ### Initialization        
        n = ini
        
        ### Evolution in time
        for t in np.arange(0,T+dt,dt): 
            t = np.round(t,3)           
            
            ### Birth and Death             
            # Stop when the population goes extinct
            if n==0 : 
                extinction[allele,i] = t
                break
            # Add births, substract deaths (mortality = 1)
            n = n + np.random.poisson(f*n*dt) - np.random.poisson(n*dt) 
            # Transform negative number of individuals into 0
            if n<0 : n=0 
            


### Save
       
# Create directory 
        
directory = f"1_without_space_dt_{dt}_s_{s}"
if not os.path.exists(directory): os.mkdir(directory)

# Save data

np.savetxt(f"{directory}/extinction.txt", extinction)

# Save figure

fig, ax = plt.subplots()
bins = [x - 0.5 for x in range(40)]
ax.hist([extinction[0,:], extinction[1,:]], bins = bins, histtype = 'bar', label = ['wt','drive'], density=True)
ax.vlines(np.mean(extinction[0,:]), 0, 0.5, color="black", linewidth=2, linestyle="-.")
ax.vlines(np.mean(extinction[1,:]), 0, 0.5, color="black", linewidth=2, linestyle="-.")
ax.set(xlabel='Time', xlim = [0,40], ylim = [0,0.35])
ax.set_title(f"Without space : extinction time.\n dt = {dt}, ini={ini}, runs = {len(extinction[0,:])}, mean wt = {np.round(np.mean(extinction[0,:]),3)} and drive = {np.round(np.mean(extinction[1,:]),3)}")
fig.savefig(f"{directory}/without_space_extinction_time.png", format='png')
plt.show() 

# Save para

file = open(f"{directory}/parameters.txt", "w") 
file.write(f"Parameters : \nini = {ini} \nT = {T} \ndt = {dt} \nr = {r} \ns = {s} \nf = {f} \nnb_i = {nb_i}")  
file.close() 











