#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 15 12:05:39 2023

@author: lena
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies

import numpy as np


# Initialization with datas
def snapshot(start, end, n_matrix, max_density, nb_drive, x_graph_values):
    # Time of the snapshot
    ti = (start+end)//2 
    # We conserve sites before (max_density) drive density.    
    last_index = np.where(n_matrix[ti,:]>max_density)[0][0]
    exp = n_matrix[ti, last_index-x_graph_values:last_index].astype(int) 
    # exp_nb is the first part of exp, before density reaches nb_drive   
    exp_nb = exp[:np.where(exp>nb_drive)[0][0]+1]
    return(exp, exp_nb)
    
# Initialization without data
def ini_exp(nb_sites, lambda_back, nb_drive, max_density, dx, x_graph_values):        
    abscisse = np.arange(-x_graph_values,1)*dx 
    exp = np.round(max_density*np.exp(lambda_back*abscisse),0).astype('int')
    exp_nb = exp[:np.where(exp>nb_drive)[0][0]+1]
    return(exp, exp_nb)
     

# Galton-watson in space
def gw(exp, exp_nb, allele, nb_i, T, dt, dx, r, s, c, h, m, dir_load, dir_save):
    
    ## Fitness according to the population simulated (idem for zyg or ger)
    f = [(r+1)*(1-s*h)*(1-c), (r+1)*(1-s)][allele]
    title = ["wt","drive"][allele]
    
    # Save parameters
    file = open(f"{dir_save}/gw_parameters.txt", "w") 
    file.write(f"Parameters : \nT = {T} \ndt = {dt} \ndx = {dx} \nr = {r} \ns = {s} \nm = {m} \nf = {f} \nnb_i = {nb_i}")  
    file.close()          
       
    # Extinction times
    extinction = np.ones(nb_i)*(-1)    
           
    # Loop on the runs    
    for i in range(nb_i) : 
            
       # Initialization   
       n = exp            
       
       # Evolution in time
       for t in np.arange(0,T+dt,dt): 
                 
            ### extinction vector : extinction times for each site of exp_nb
            
            # -1 means that there is still individuals on the site
            # a positive value is the time at which the last individual has been seen            
            if extinction[i] == -1 and n[len(exp_nb)-1]==0:
                extinction[i] = t 
            if extinction[i] != -1 and n[len(exp_nb)-1]!=0:
                extinction[i] = -1
                
            ### Stop the simulation if there is nobody left    
            
            if np.sum(n) == 0 : 
                #print("t =", t)
                break
            
            ### Birth and Death   
        
            # Add births, substract deaths (mortality = 1)
            n = n + np.random.poisson(f*n*dt) - np.random.poisson(n*dt)    
            # Transform negative number of individuals into 0
            n[np.where(n<0)[0]]=0
                    
            ### Migration  
    
            # Number of migrants in each site
            n_mig = np.random.binomial(n,m)
            # Half migrate to the right, half to the left
            n_mig_left = np.random.binomial(n_mig,0.5); n_mig_right = n_mig - n_mig_left
            # Substract the migrants leaving
            n -= n_mig
            # ... except for those going outside the windows (they stay home)
            n[0] += n_mig_left[0]
            n[-1] += n_mig_right[-1]
            # Add the migrants in the neighboor sites
            n[1:] += n_mig_right[:-1]
            n[:-1] += n_mig_left[1:]
                
    # Save data  
    np.savetxt(f"{dir_save}/gw_{title}_extinction.txt", extinction)
    return(extinction, exp)













