#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os

# Change font to serif
plt.rcParams.update({'font.family':'serif'})

# Graph parameters 
#title_size = 20
label_size = 20
legend_size = 17
number_y_size = 17


def erad_s(conv_timing, r, c, h) : 
    print(f"\nEradication pour {np.round(r/(1+r),3)} < s < 1") 
    s_1 = c/(1-h*(1-c))   
    if conv_timing == "zyg" :
        s_2 = c/(2*c + h*(1-c))
        if s_1 < s_2 : 
            print("Coexistence", ":", np.round(s_1,3), "< s <", np.round(s_2,3))
        else :  
            print("Bistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        lin = c*(1-2*s)-(1-c)*s*h
        if lin > 0 : print("Linear speed :", 2*np.sqrt(lin))               
    if conv_timing == "ger" : 
        s_2 = c/(2*c*h + h*(1-c))
        if s_1 < s_2 :
            print("Coexistence", ":", np.round(s_1,3) , "< s <", np.round(s_2,3))
        else :  
            print("Bistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        lin = c*(1-2*s*h)-(1-c)*s*h 
        if lin > 0 : print("Linear speed :", 2*np.sqrt(lin))
    return(2*np.sqrt(lin))


# /!\ r=0 empêche une vrai explosion après recolonisation !

K = 10**8        # Carrying capacity for a spatial interval of size 1
dx = 1                      # Spatial interval 
T = 1000                     # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)
m = 0.2                     # Migration probability
dt = np.round(m*dx**2/2,10) # Time interval
nb_graph = 100                 # Number of graphs shown if T is reached 
nb_graph_kimo = 2000
nb_save = int(T/dt) #1000   # Number of values (nD, nW) saved in between 0...T (In zoom all values are saved)


conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
r = 0.1                     # Intrasic growth rate
c = 0.9                     # Conversion rate
h = 0.4                     # Dominance
s = 0.3                     # Disadvantage for drive

v_cont = erad_s(conv_timing, r, c, h) 
nb_sites = 1000 # int(((v_cont*T*2/dx)//1000)*1000+1000) 

kymo = False

# Where to save datas
dir_save = f"../../../stoch_not_save/datas/1D_{conv_timing}_K_{int(np.log10(K))}_dx_{dx}_s_{s}_r_{r}"
if not os.path.exists(dir_save ): os.mkdir(dir_save)


### Initialization

nD = np.zeros(nb_sites).astype(int); nD[:nb_sites//2] = K*dx
nW = np.zeros(nb_sites).astype(int); nW[nb_sites//2:] = K*dx
nD_matrix = np.zeros((nb_save+1,nb_sites)).astype(int)
nW_matrix = np.zeros((nb_save+1,nb_sites)).astype(int)
time = np.zeros(nb_save+1)
fD = np.zeros(nb_sites)
fW = np.zeros(nb_sites)
graph_counter = 0
graph_counter_kimo = 0
save_counter = 0

# kimo
ND = np.zeros([nb_graph_kimo,nb_sites])
NW = np.zeros([nb_graph_kimo,nb_sites])


### Evolution in time
for t in np.arange(0, T+dt, dt): 
    t = np.round(t,3)  
        
    ### GRAPH 1 Histogram (density function of space) or lines
    
    if t >= graph_counter*(T/nb_graph):
        fig, ax = plt.subplots(figsize = (10, 5))
        bins = [x - 0.5 for x in range(0, nb_sites+1)]
        #plt.hist([np.arange(nb_sites), np.arange(nb_sites)], bins = bins, weights = [nD, nW], 
        #          histtype = 'barstacked', label = ['Drive','Wild-type'], color = ['crimson','cornflowerblue']) 
        plt.plot(np.arange(nb_sites), nW, label ='Wild-type', color = 'cornflowerblue', linewidth = 5) 
        plt.plot(np.arange(nb_sites), nD, label = 'Drive', color = 'crimson', linewidth = 5) 
        #ax.set_xlabel('Space', fontsize = label_size)
        #ax.set_ylabel('Number of alleles', fontsize = label_size)
        ax.set(ylim = [0,1.1*K*dx], xlim = [0,nb_sites])          
        #ax.yaxis.set_tick_params(labelsize=number_y_size)
        #ax.set_title(f"Time {t}")
        plt.tick_params(labelbottom = False, bottom = False, labelleft = False, left = False,)
        #plt.legend(loc="upper left", fontsize = legend_size)
        plt.tight_layout()
        fig.savefig(f"{dir_save}/1D_K_{int(np.log10(K))}_s_{s}_t_{t}.png", format='png')
        #fig.savefig(f"{dir_save}/t_{t}.svg", format='svg')
        plt.close()
        graph_counter += 1  
        
    ### GRAPH 2 kimograph (time function of space)
    
    if kymo : 
        if t >= graph_counter_kimo*(T/nb_graph_kimo):
            ND[nb_graph_kimo - graph_counter_kimo-1,:] = nD
            NW[nb_graph_kimo - graph_counter_kimo-1,:] = nW
            graph_counter_kimo += 1        
         
               
    ### Save values (nD, nW) in a matrix, in between 0 ... T
        
    if t >= save_counter*(T/nb_save):
        nD_matrix[save_counter,:] = nD
        nW_matrix[save_counter,:] = nW
        time[save_counter] = t
        save_counter += 1
        #print(t)
        
        
    ### Stop the simulation if the wave goes outside the window 
        
    if np.where(nD==max(nD))[0][0] > len(nD)-10 : 
        print("t =",t)
        break
    
     
        ### Birth and Death   
        
        # Index for empty and non empty sites
    extinct_index = np.where(nD+nW==0)[0]
    survive_index = np.delete(np.arange(nb_sites), extinct_index)
    # For non empty site, the fecundity is given by the following.
    sv_pop = nD[survive_index] + nW[survive_index]; sv_nD = nD[survive_index]; sv_nW = nW[survive_index]
    
    if conv_timing == "zyg" :
        fD[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-s)*sv_nD + (1-s*h)*(1-c)*sv_nW +2*c*(1-s)*sv_nW ) /sv_pop
        fW[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-c)*(1-s*h)*sv_nD + sv_nW ) /sv_pop           
    if conv_timing == "ger" : 
        fD[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-s)*sv_nD + (1-s*h)*(1+c)*sv_nW ) /sv_pop
        fW[survive_index] = ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-c)*(1-s*h)*sv_nD + sv_nW ) /sv_pop
            
    # For empty site, the fecundity is 0.
    fD[extinct_index] = 0
    fW[extinct_index] = 0
    # Check that all fecundity values are numbers.
    if len(np.argwhere(np.isnan(fD))) != 0 or len(np.argwhere(np.isnan(fW))) != 0 : 
        print("Houston, we have a problem")
    # Add births, substract deaths (mortality = 1)
    nD = nD + np.random.poisson(fD*nD*dt) - np.random.poisson(nD*dt)            
    nW = nW + np.random.poisson(fW*nW*dt) - np.random.poisson(nW*dt)
    # Transform negative number of alleles into 0
    nD[np.where(nD<0)[0]]=0
    nW[np.where(nW<0)[0]]=0
    
     
    ### Migration  
    
    # Number of migrants in each site
    nD_mig = np.random.binomial(nD,m)
    nW_mig = np.random.binomial(nW,m)
    # Half migrate to the right, half to the left
    nD_mig_left = np.random.binomial(nD_mig,0.5); nD_mig_right = nD_mig - nD_mig_left
    nW_mig_left = np.random.binomial(nW_mig,0.5); nW_mig_right = nW_mig - nW_mig_left
    # Substract the migrants leaving
    nD -= nD_mig 
    nW -= nW_mig
    # ... except for those going outside the windows (they stay home)
    nD[0] += nD_mig_left[0]; nW[0] += nW_mig_left[0]
    nD[-1] += nD_mig_right[-1]; nW[-1] += nW_mig_right[-1]
    # Add the migrants in the neighboor sites
    nD[1:] += nD_mig_right[:-1]; nW[1:] += nW_mig_right[:-1] 
    nD[:-1] += nD_mig_left[1:]; nW[:-1] += nW_mig_left[1:]
    
    

# Save datas 
time[0]=-1; time = time[np.where(time!=0)[0]]; time[0]=0  #remove all the zeros at the end, when T is not reached.
nD_matrix = nD_matrix[:len(time),:]
nW_matrix = nW_matrix[:len(time),:]
np.savetxt(f"{dir_save}/time.txt", time)
np.savetxt(f"{dir_save}/nD_matrix.txt", nD_matrix)
np.savetxt(f"{dir_save}/nW_matrix.txt", nW_matrix)
file = open(f"{dir_save}/parameters.txt", "w") 
file.write(f"Parameters : \nK = {K} \nnb_sites = {nb_sites} \ndx = {dx} \nT = {T} \ndt = {dt} \nm = {m} \nconv_timing = {conv_timing} \nr = {r} \ns = {s}  \nc = {c} \nh = {h} \nnb_graph = {nb_graph} \nnb_save = {nb_save}")  
file.close() 
    


if kymo : 
    
## New color scales
    # cornblueflower
    N = 256
    rgba = np.array([100/255,149/255,237/255, 1])
    vals = np.ones((N, 4))
    for i in range(3):
        vals[:, i] = np.linspace(1,rgba[i], N)
    vals[:, 3] = np.linspace(0, 1, N)
    newcmpblue = ListedColormap(vals)
    #crimson
    rgba = np.array([200/255,0/255,10/255, 1])
    vals = np.ones((N, 4))
    for i in range(3):
        vals[:, i] = np.linspace(1,rgba[i], N)
    vals[:, 3] = np.linspace(0, 1, N)
    newcmpred = ListedColormap(vals)


## Kymograph visualisation
    fig, ax = plt.subplots()
    im1 = ax.imshow(NW, cmap = newcmpblue,vmin = 0, vmax = K*dx/10, alpha=0.9, aspect='auto')
    ax.figure.colorbar(im1, ax=ax) 
    im2 = ax.imshow(ND, cmap = newcmpred, vmin = 0, vmax = K*dx/20,alpha=1, aspect='auto')
    ax.figure.colorbar(im2, ax=ax) 
    plt.xlabel('Space', fontsize=label_size)
    plt.ylabel('Time', fontsize=label_size)
    plt.tick_params(labelleft = False, labelbottom = False, left = False, bottom = False) 
    plt.tight_layout()
    #plt.show()
    fig.savefig(f"{dir_save}/kymograph_s_{s}.png", format='png')
    #fig.savefig(f"{dir_save}/kymograph_s_{s}.svg", format='svg')
    plt.close()
    
    
    
    
