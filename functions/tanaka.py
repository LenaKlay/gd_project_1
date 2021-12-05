#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:04:17 2021

@author: lena
"""


# Libraries
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg  as la
import matplotlib.pyplot as plt

                  

def tanaka(s,model,model_para,num_para,graph_para): 
    T,L,M,N,theta = num_para
    CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating = model_para
    wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion = graph_para

    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    X = np.linspace(0,N,N+1)*dx   
    
    # Initialization       
    if CI == "center" : 
        P = np.zeros(N+1); P[0:N//2] = 0.95*max_capacity                # Proportion of drive individuals at t=0
    if CI == "left" :      
        P = np.zeros(N+1); P[0:N//10] = 0.95*max_capacity               # Proportion of drive individuals at t=0
    
    if graph_type != None :
        fig, ax = plt.subplots()  
        ax.plot(X, P, label=f'Tanaka {model}')
        plt.legend();plt.grid()
        plt.show()
    nb_graph = 1

    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    B = sp.identity(N+1)+((1-theta)*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    B_ = sp.identity(N+1)-(theta*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  

    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.
    vitesse_en_fct_du_tps = np.zeros((2,T//mod))   # first line : time, second line : speed of the wave at the corresponding time
    
    # Evolution
    for t in np.linspace(dt,T,M) : 
        
        t = round(t,2)
        if model == "cubic" : 
            f = s*P*(1-P)*(P-(2*s-1)/s)
        if model == "fraction" : 
            f = (s*P*(1-P)*(P-(2*s-1)/s))/(1-s+s*(1-P)**2)       
      
        P = la.spsolve(B_, B.dot(P) + dt*f)

        if t>=mod*nb_graph and graph_type != None :  
            fig, ax = plt.subplots()  
            ax.plot(X, P, label=f'Tanaka {model}')
            plt.legend();plt.grid()
            plt.show()
            nb_graph += 1
            
        if np.isin(True, (1-P)>0.5) and np.isin(True, (1-P)<0.99) :             # wild pop is still present somewhere in the environment
            position = np.append(position,np.where((1-P)>0.5)[0][0])            # List containing the first position where the 
                                                                                # proportion of wild alleles is lower than 0.5.               
        if t%mod == 0 : 
            #print("\nTemps :", t) ; print("Vitesse :", np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt, "\n")
            vitesse_en_fct_du_tps[:,int(t)//mod-1] = np.array([t,np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt])
            # first line : time, second line : speed of the wave at the corresponding time
                                                                                                                                                 
    if np.shape(position)[0] != 0 :
        position = position[position!=0]                          
        speed = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave     
    else :
        speed = None
        print(f"Can't determine the speed of the wave for s = {s}.")                                                 
            
    return(P, speed, vitesse_en_fct_du_tps)

