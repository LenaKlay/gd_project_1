#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la

### Graphs parameters 

plt.rcParams.update({'font.family':'serif'})
label_size = 12
legend_size = 10

    
### Parameters
    
dx = 1                      # Spatial interval 
T = 1000                    # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)
m = 0.2                     # Migration probability
dt = np.round(m*dx**2/2,10) # Time interval

conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
r = 0.1                     # Intrasic growth rate
c = 0.9                     # Conversion rate
h = 0.4                     # Dominance

K_range = np.arange(4,9)
s_range = np.arange(1,10,0.5)
     
nb_drive = 100
discrete = False
distance_matrix = np.zeros((len(K_range), len(s_range)))

if discrete : print("\nDiscrete\n")
else : print("\nContinuous\n")

for i in range(0, len(K_range)):
    for j in range(len(s_range)) :

        K = 10**(K_range[i])      # Carrying capacity for a spatial interval of size 1
        s = np.round(0.1*(s_range[j]), 3)# Disadvantage for drive
        
        print("----- K =", K, ", s =", s)
            
        # Speed and space for the wave not to go outside the window    
        s_1 = c/(1-h*(1-c))   
        s_2 = c/(2*c*h + h*(1-c))
        lin = c*(1-2*s*h)-(1-c)*s*h     
        v_cont = 2*np.sqrt(lin)
        nb_sites = int(((v_cont*T*2/dx)//1000)*1000+1000)
        
        distance_nb = np.zeros(int((T//2)/dt))
        nb_save=0
            
        nD = np.zeros(nb_sites).astype(int); nD[:nb_sites//2] = K*dx
        nW = np.zeros(nb_sites).astype(int); nW[nb_sites//2:] = K*dx
        fD = np.zeros(nb_sites)
        fW = np.zeros(nb_sites)
        
        if not discrete :
            # Matrix
            theta = 0.5
            C0 = -2*np.ones(nb_sites-2); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
            C1 = np.ones(nb_sites-2) 
            A = sp.spdiags([C1,C0,C1],[-1,0,1], nb_sites-2, nb_sites-2)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)             
            B = sp.identity(nb_sites-2)+((1-theta)*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme         
            B_ = sp.identity(nb_sites-2)-(theta*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  
          
            
        ### Evolution in time
        
        for t in np.arange(0, T, dt): 
                t = np.round(t,3)  
            
        
                ### Stop the simulation if the wave goes outside the window 
                
                if np.where(nD==max(nD))[0][0] > len(nD)-10 or t==T-dt : 
            
                    print("t =",t)
                    
                    distance_nb = distance_nb[np.where(distance_nb!=0)[0]]
                    print("distance =", np.mean(distance_nb))
                    distance_matrix[i,j] = np.mean(distance_nb)
                    print(distance_matrix) 
                
                    # Save datas 
                    file = open(f"dis_K_{int(np.log10(K))}_s_{s}.txt", "a") 
                    file.write(f"\n{np.mean(distance_nb)}")  
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
            
                
                if t > T/2 :
                    distance_nb[nb_save] = (np.where(nW>nb_drive)[0][0] - np.where(nD>nb_drive)[0][0])*dx
                    nb_save = nb_save+1
    
     
                if discrete: 
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
                
                else:
                    # Index for empty and non empty sites
                    extinct_index = np.where(nD+nW==0)[0]
                    survive_index = np.delete(np.arange(nb_sites), extinct_index)
                    # For non empty site, the fecundity is given by the following.
                    sv_pop = nD[survive_index] + nW[survive_index]; sv_nD = nD[survive_index]; sv_nW = nW[survive_index]
                
                    # For empty site, the fecundity is 0.
                    fD[extinct_index] = 0
                    fW[extinct_index] = 0
                    
                    #fD[survive_index] = sv_nD*((r*(K*dx-sv_pop)/(K*dx)+1)*((1-s)*(sv_nD/sv_pop)+(1-s*h)*(1+c)*(sv_nW/sv_pop))-1)
                    #fW[survive_index] = sv_nW*((r*(K*dx-sv_pop)/(K*dx)+1)*((sv_nW/sv_pop)+(1-s*h)*(1-c)*(sv_nD/sv_pop))-1)
                    
                    
                    fD[survive_index] = sv_nD*( ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-s)*sv_nD + (1-s*h)*(1+c)*sv_nW ) /sv_pop - 1 )
                    fW[survive_index] = sv_nW*( ( 1 + r*(1-sv_pop/(K*dx)) ) * ( (1-c)*(1-s*h)*sv_nD + sv_nW ) /sv_pop  - 1 )
                    
                    
                    nD[1:-1] = la.spsolve(B_, B.dot(nD[1:-1]) + dt*fD[1:-1])
                    nW[1:-1] = la.spsolve(B_, B.dot(nW[1:-1]) + dt*fW[1:-1])
                    
                    nW[0] = nW[1]; nD[0] = nD[1]    # Neumann condition alpha=0
                    nW[-1] = nW[-2]; nD[-1] = nD[-2]  # Neumann condition beta=0
                     
        ### Graph   

        graph = False
        if graph: 
             fig, ax = plt.subplots()
             bins = [x - 0.5 for x in range(0, nb_sites+1)]
             plt.hist([np.arange(nb_sites), np.arange(nb_sites)], bins = bins, weights = [nD, nW], 
                       histtype = 'barstacked', label = ['drive','wt'], color = ['crimson','cornflowerblue'])  
             ax.set(xlabel='Space', ylabel='Number of individuals', ylim = [0,1.1*K*dx])
             ax.set_title(f"Time {t}")
             plt.legend(loc="upper left")
             plt.show()
            
            
graph = True
if graph : 
            col = ["plum", "orchid", "m", "darkviolet", "indigo", "navy"]    
            fig, ax = plt.subplots()
            for i in (range(len(K_range)))[::-1]:
                ax.plot(s_range*0.1, distance_matrix[i,:], label = f"$K=10^{int(K_range[i])}$", color = col[i])
                ax.scatter(s_range*0.1, distance_matrix[i,:], alpha=0.5, color = col[i])
            ax.set(xlabel='s (fitness disadvantage for drive)', ylabel='Distance') 
            ax.set_ylim([0,100])
            ax.xaxis.label.set_size(label_size)
            ax.yaxis.label.set_size(label_size)
            plt.rc('legend', fontsize=legend_size)  
            plt.legend()   
            fig.savefig(f"distance_s.png", format='png')
            plt.show()
    

file = open(f"save_dis_{discrete}.txt", "a") 
file.write(f"{distance_matrix}")  
file.close() 
