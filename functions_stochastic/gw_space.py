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

# Space initialization

def gw(time, start, end, n_matrix, allele, max_density, left_margin, nb_drive, m, nb_i, dx, T, dt, r, s, c, h, lambda_back, snapshot, dir_load, dir_save):
    
    ## Fitness according to the population simulated (idem for zyg or ger)
    f = [(r+1)*(1-s*h)*(1-c), (r+1)*(1-s)][allele]
    title = ["wt","drive"][allele]
    col = ["cornflowerblue", "crimson"][allele]
    ## Initial condition    
    
    # Save parameters
    file = open(f"{dir_save}/gw_parameters.txt", "w") 
    file.write(f"Parameters : \nT = {T} \ndt = {dt} \ndx = {dx} \nr = {r} \ns = {s} \nm = {m} \nf = {f} \nnb_i = {nb_i}")  
    file.close() 
          
    # We take a snapshop of the wave for the CI
    if snapshot :  
        # Time of the snapshot
        ti = (int(3*time[-1]/5)+ int(time[-1]-10) )//2 
        # We conserve sites before (max_density) drive density.    
        exp = n_matrix[ti, np.where(n_matrix[ti,:]>0)[0][0]-left_margin:np.where(n_matrix[ti,:]>max_density)[0][0]].astype(int) 
        abscisse = np.arange(-(len(exp)-1),1)*dx
        # exp_nb is the first part of exp, before density reaches nb_drive   
        exp_nb = n_matrix[ti, np.where(n_matrix[ti,:]>0)[0][0]-left_margin:np.where(n_matrix[ti,:]>nb_drive)[0][0]].astype(int)     
        print(f"Density close to {nb_drive} :", exp_nb[-1])
        
    # We use the theoritical exponential profil for the CI
    else :
        abscisse = np.arange(-len(n_matrix[0,:]),1)*dx 
        exp = np.round(max_density*np.exp(lambda_back*abscisse),0).astype('int')
        abscisse = abscisse[np.where(exp>0)[0][0]-left_margin:]; exp = exp[np.where(exp>0)[0][0]-left_margin:]
        exp_nb = exp[:np.where(exp>nb_drive)[0][0]]
       
    
    #fig, ax = plt.subplots()
    #plt.plot(abscisse, exp, color = col)
    #ax.set(xlabel='Space', ylabel='Number of individuals')
    #ax.set_title(f"Density : {title}")
    #fig.savefig(f"{dir_save}/gw_{title}_ini.png", format='png')
    #plt.show()
     
    # Plot initial conditions
    fig, ax = plt.subplots()
    bins = [x - 0.5 for x in range(len(exp)+1)]
    plt.hist(np.arange(len(exp)), bins = bins, weights = exp,
             histtype = 'barstacked', color = col) 
    ax.set(xlabel=f'Spatial sites', ylabel='Number of individuals')
    ax.set_title(f"Density : {title}")
    fig.savefig(f"{dir_save}/gw_{title}_ini.png", format='png')
    ax.set_xticks(np.linspace((len(exp)-1)%10, len(exp)-1, len(np.arange(abscisse[0], abscisse[-1]+1, 10))))    
    ax.set_xticklabels(np.arange(abscisse[(len(exp)-1)%10], abscisse[-1]+1, 10))    
    plt.show()
       
    # Extinction times
    extinction = np.ones(nb_i)*(-1)    
    #extinction_tps_entier = np.ones(nb_i)*(-1)  
    #graph_counter = 0  
           
    # Loop on the runs    
    for i in range(nb_i) : 
            
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
                
                
                ### extinction vector : extinction times for each site of exp_nb
                # -1 means that there is still individuals on the site
                # a positive value is the time at which the last individual has been seen
                
            if extinction[i] == -1 and n[len(exp_nb)-1]==0:
                extinction[i] = t 
            if extinction[i] != -1 and n[len(exp_nb)-1]!=0:
                extinction[i] = -1
                
            #if t == int(t) : 
            #    if extinction_tps_entier[i] == -1 and n[len(exp_nb)-1]==0:
            #        extinction_tps_entier[i] = t 
            #    if extinction_tps_entier[i] != -1 and n[len(exp_nb)-1]!=0:
            #        extinction_tps_entier[i] = -1
    
        
                ### Stop the simulation if there is nobody left
                
            if np.sum(n) == 0 : # np.sum(n[0:len(exp_nb)-1+100]) == 0 : 
                #print("t =", t)
                #print("time extinction =", extinction[i])
                #if extinction_tps_entier[i] == -1 :
                #    print(i, "on a tps d'ext entier de -1")
                    #extinction_tps_entier[i] = 1
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
            #    vect = n[len(exp_nb)-1-30:len(exp_nb)-1+30]
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
        
    #extinction = extinction[:,:len(exp_nb)]
    np.savetxt(f"{dir_save}/gw_{title}_extinction.txt", extinction)
    #np.savetxt(f"{dir_save}/gw_{wt_or_drive}_extinction_tps_entier.txt", extinction_tps_entier)
    
    # Save figure            
    #fig, ax = plt.subplots()
    #ax.hist(extinction, bins=range(40), histtype = 'bar', label = 'wt', density=True)
    #ax.vlines(np.mean(extinction), 0, 0.5, color="black", linewidth=2, linestyle="-.")
    #ax.set(xlabel='Time', xlim = [0,40], ylim = [0,0.5])
    #ax.set_title(f"extinction time of {wt_or_drive}.\n m={m}, mean={np.round(np.mean(extinction),3)} ")
    #fig.savefig(f"{dir_save}/gw_{wt_or_drive}_m_{m}_extinction_time.png", format='png')
    #plt.show() 
    
    #fig, ax = plt.subplots()
    #ax.hist(extinction_tps_entier, bins=range(40), histtype = 'bar', label = 'wt', density=True)
    #ax.vlines(np.mean(extinction_tps_entier), 0, 0.5, color="black", linewidth=2, linestyle="-.")
    #ax.set(xlabel='Time', xlim = [0,40], ylim = [0,0.5])
    #ax.set_title(f"Extinction time of {wt_or_drive}.\n m={m}, mean={np.round(np.mean(extinction_tps_entier),3)} ")
    #fig.savefig(f"{dir_save}/gw_{wt_or_drive}_m_{m}_extinction_tps_entier.png", format='png')
    #plt.show() 
    # Save para        
    
    
    return(extinction)













