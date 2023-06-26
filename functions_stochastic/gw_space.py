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
import matplotlib.pyplot as plt
import os


# Space initialization

def homo_exp(max_density, m_range, nb_i, wt_or_drive, dx, T, dt, r, s, snapshot, dir_load, dir_save):
    
    ## Fitness according to the population simulated
    
    if wt_or_drive == "nD": f = 0
    elif wt_or_drive == "nW": f = 1-s
    
        
    ## Initial condition    
        
    # We take a snapshop of the wave for the CI
    if snapshot :  
        # Download the densities
        n_matrix = np.loadtxt(f"{dir_load}/{wt_or_drive}_matrix.txt")   # Matrix for the density : row = time, column = space
        # Download the time
        time = np.loadtxt(f"{dir_load}/time.txt")                       # Time vector
        # Time of the snapshot
        ti = (int(3*time[-1]/5)+ int(time[-1]-10) )//2 
        # We conserve sites before (max_density) drive density.    
        exp = n_matrix[ti, np.where(n_matrix[ti,:]>0)[0][0]-100:np.where(n_matrix[ti,:]>max_density)[0][0]].astype(int) 
        # exp100 is the first part of exp, before density reaches 100   
        exp100 = n_matrix[ti, np.where(n_matrix[ti,:]>0)[0][0]-100:np.where(n_matrix[ti,:]>100)[0][0]].astype(int)     
        print("Density close to 100 :", exp100[-1]); index100 = len(exp100)-1
        
    # We use the theoritical exponential profil for the CI
    #else :
        # ...... à compléter
                

    
    # Plot initial conditions
    fig, ax = plt.subplots()
    bins = [x - 0.5 for x in range(len(exp)+1)]
    plt.hist(np.arange(len(exp)), bins = bins, weights = exp,
             histtype = 'barstacked', color = 'cornflowerblue') 
    ax.set(xlabel=f'Spatial sites', ylabel='Number of individuals')
    ax.set_title(f"Density : {wt_or_drive}")
    fig.savefig(f"{dir_save}/{wt_or_drive}_ini.png", format='png')
    plt.show()

    # Loop on migration rates
    for m in m_range : 
        
        # Extinction times
        extinction = np.ones(nb_i)*(-1)    
        extinction_tps_entier = np.ones(nb_i)*(-1)  
        #graph_counter = 0  
        print("m =", m, wt_or_drive)
           
        # Loop on the runs    
        for i in range(nb_i) : 
            print(i)
            
            # Initialization   
            n = exp            
            
            # Evolution in time
            for t in np.arange(0,T+dt,dt): 
                
                ## Draft
                # Density previously zero and still zero : time stay inchanged
                # old_zeros = np.intersect1d(np.where(n==0)[0], np.where(extinction[i,:]>0)[0])
                # Density just turned to zero 
                #new_zeros  = np.intersect1d(np.where(n==0)[0], np.where(extinction[i,:]==-1)[0])  
                #extinction[i, new_zeros] = t 
                # Strictly positive density
                #extinction[i, np.where(n!=0)[0]] = -1
                
                
                ### extinction vector : extinction times for each site of exp100
                # -1 means that there is still individuals on the site
                # a positive value is the time at which the last individual has been seen
                
                if extinction[i] == -1 and n[index100]==0:
                    extinction[i] = t 
                if extinction[i] != -1 and n[index100]!=0:
                    extinction[i] = -1
                    
                if t == int(t) : 
                    if extinction_tps_entier[i] == -1 and n[index100]==0:
                        extinction_tps_entier[i] = t 
                    if extinction_tps_entier[i] != -1 and n[index100]!=0:
                        extinction_tps_entier[i] = -1
    
        
                ### Stop the simulation if there is nobody left
                
                if np.sum(n) == 0 : # np.sum(n[0:index100+100]) == 0 : 
                    #print("t =", t)
                    #print("time extinction =", extinction[i])
                    if extinction_tps_entier[i] == -1 :
                        print(i, "on a tps d'ext entier de -1")
                        extinction_tps_entier[i] = 1
                    #print("time extinction temps entier =", extinction_tps_entier[i])
                    break
            
                ### Birth and Death   
            
                # Check that all fecundity values are numbers.
                #if len(np.argwhere(np.isnan(f))) != 0 : 
                #    print("Houston, we have a problem")
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
                
                # Plot
                #if i==0 and t >= graph_counter*(T/nb_graph):
                #    vect = n[index100-30:index100+30]
                #    #vect = np.ones(len(n))*(-1); survive = np.where(n!=0)
                #    #vect[survive] = np.log(n[survive])
                #    fig, ax = plt.subplots()
                #    bins = [x - 0.5 for x in range(len(vect)+1)]
                #    plt.hist(np.arange(len(vect)), bins = bins, weights = vect,
                #             histtype = 'barstacked', color = 'cornflowerblue')  
                #    ax.set(xlabel=f'Spatial sites', ylabel='Number of individuals')
                #    ax.set_title(f"Density at time {t}")
                #    ax.set(ylim = [0,15])
                #    fig.savefig(f"test_ext_time/{t}.png", format='png')  
                #    ax.set(ylim = [0,15])
                #    plt.show()
                #    graph_counter += 1

        
        # Save data  
        
        #extinction = extinction[:,:len(exp100)]
        np.savetxt(f"{dir_save}/{wt_or_drive}_m_{m}_extinction.txt", extinction)
        np.savetxt(f"{dir_save}/{wt_or_drive}_m_{m}_extinction_tps_entier.txt", extinction_tps_entier)
        
        # Save figure            
        fig, ax = plt.subplots()
        ax.hist(extinction, bins=range(40), histtype = 'bar', label = 'wt', density=True)
        ax.vlines(np.mean(extinction), 0, 0.5, color="black", linewidth=2, linestyle="-.")
        ax.set(xlabel='Time', xlim = [0,40], ylim = [0,0.5])
        ax.set_title(f"extinction time of {wt_or_drive}.\n m={m}, mean={np.round(np.mean(extinction),3)} ")
        fig.savefig(f"{dir_save}/{wt_or_drive}_m_{m}_extinction_time.png", format='png')
        plt.show() 
        
        fig, ax = plt.subplots()
        ax.hist(extinction_tps_entier, bins=range(40), histtype = 'bar', label = 'wt', density=True)
        ax.vlines(np.mean(extinction_tps_entier), 0, 0.5, color="black", linewidth=2, linestyle="-.")
        ax.set(xlabel='Time', xlim = [0,40], ylim = [0,0.5])
        ax.set_title(f"Extinction time of {wt_or_drive}.\n m={m}, mean={np.round(np.mean(extinction_tps_entier),3)} ")
        fig.savefig(f"{dir_save}/{wt_or_drive}_m_{m}_extinction_tps_entier.png", format='png')
        plt.show() 
        # Save para        
        file = open(f"{dir_save}/parameters.txt", "w") 
        file.write(f"Parameters : \nT = {T} \ndt = {dt} \ndx = {dx} \nr = {r} \ns = {s} \nm = {m} \nf = {f} \nnb_i = {nb_i}")  
        file.close() 
            





### Parameters and datas ###
        
        
## From the datas
    
# Choose dx, s and m to load datas 
dx = 0.1                    # Space interval
s = 0.4                     # Disadvantage for drive 
m = 0.2                     # Migration rate


# Load the other parameters
dir_load = f"dx_{dx}_s_{s}_m_{m}"
file = open(f"{dir_load}/parameters.txt", "r")
para = file.read()
para_list = para.replace(' ', '').split("\n")[:-1]
print(para_list)
K = float(para_list[1].replace('K=', ''))               # Carrying capacity in one unit of space 
dt = float(para_list[5].replace('dt=', ''))             # Time interval 
r = float(para_list[7].replace('r=', ''))               # Intrasic growth rate   
nb_sites = float(para_list[2].replace('nb_sites=', ''))               # Intrasic growth rate   
file.close()  

## Simulation

# Parameters
T = 100                     # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)                               
nb_i = 100                  # Number of replicas
#nb_graph = 100             # Number of graphs
m_range = [0.2]             # Migration rates tested in this simumation (could be a few)
wt_or_drive = "nW"          # Wild-type (nW) or Drive (nD) population simulated 
snapshot = True             # CI created from a snapshot of the original simulation, of from the theortical exponential value
max_density = 10**7         # Maximal density in the exponential CI

# Create a directory to save data and figures                 
dir_save = f"2_GW_exp_dx_{dx}_s_{s}"
if not os.path.exists(dir_save): os.mkdir(dir_save)

# Run the simulation
homo_exp(max_density, m_range, nb_i, wt_or_drive, dx, T, dt, r, s, snapshot, dir_load, dir_save)







