#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies

import numpy as np

### Parameters

K = 10**(3) # 1000000000       # Carrying capacity for a spatial interval of size 1
dx = 1                      # Spatial interval 
T = 1000                    # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)
m = 0.2                     # Migration probability
dt = np.round(m*dx**2/2,10) # Time interval

conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
r = 0.1                     # Intrasic growth rate
c = 0.9                     # Conversion rate
h = 0.4                     # Dominance
s = 0.1*(7)                   # Disadvantage for drive

replicats = 100               # Number of runs for the simulation


print("----- K =", K, ", s =", s)
        
# Speed and space for the wave not to go outside the window    
s_1 = c/(1-h*(1-c))   
s_2 = c/(2*c*h + h*(1-c))
lin = c*(1-2*s*h)-(1-c)*s*h     
v_cont = 2*np.sqrt(lin)
nb_sites = int(((v_cont*T*2/dx)//1000)*1000+1000)
       
for run in range(replicats):
            print(run)

            ### Initialization
        
            nD = np.zeros(nb_sites).astype(int); nD[:nb_sites//2] = K*dx
            nW = np.zeros(nb_sites).astype(int); nW[nb_sites//2:] = K*dx
            fD = np.zeros(nb_sites)
            fW = np.zeros(nb_sites)
            chasing = 0
            
            ### Evolution in time
            
            for t in np.arange(0, T, dt): 
                t = np.round(t,3)  
            
        
                ### Stop the simulation if the wave goes outside the window 
                
                if np.where(nD==max(nD))[0][0] > len(nD)-10 or t==T-dt : 
            
                    print("t =",t)
                    
                    if chasing < 0.05*(t-T/2)/dt : chasing = 0
                    else : chasing = 1
                    print ("chasing =", chasing)
                
                    # Save datas 
                    file = open(f"chasing/chasing_K_{int(np.log10(K))}_s_{s}.txt", "a") 
                    file.write(f"\n{chasing}")  
                    file.close() 
                    
                    #if run == 0 : 
                        #fig, ax = plt.subplots()
                        #bins = [x - 0.5 for x in range(0, nb_sites+1)]
                        #plt.hist([np.arange(nb_sites), np.arange(nb_sites)], bins = bins, weights = [nD, nW], 
                        #              histtype = 'barstacked', label = ['drive','wt'], color = ['crimson','cornflowerblue'])  
                        #ax.set(xlabel='Space', ylabel='Number of individuals', ylim = [0,1.1*K*dx])
                        #ax.set_title(f"Time {t}, K = 10^{int(np.log10(K))}, s = {s}")
                        #plt.legend(loc="upper left")
                        #plt.show()
        
                    break
            
                ### Count the chasing times
                
                if T > T/2 and np.where(nW>0)[0][0] < np.where(nD>0)[0][0]:
                    chasing += 1
    
     
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
                # Transform negative number of individuals into 0
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
                
