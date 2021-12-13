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

# Graph parameters 
legend_size = 15
label_size = 20
number_on_x_axe = False
number_x_size = 10
number_y_size = 20
line_size = 4
                  

def tanaka(s,model,model_para,num_para,graph_para): 
    T,L,M,N,mod,theta = num_para
    CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating = model_para
    graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_figure, speed_proportion = graph_para

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
        graph_tanaka(X,P,0,model,graph_para)
    nb_graph = 1

    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    B = sp.identity(N+1)+((1-theta)*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    B_ = sp.identity(N+1)-(theta*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  

    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.
    speed_fct_of_time = np.zeros((2,T//mod))   # first line : time, second line : speed of the wave at the corresponding time
    
    # Evolution
    for t in np.linspace(dt,T,M) : 
        
        t = round(t,2)
        if model == "cubic" : 
            f = s*P*(1-P)*(P-(2*s-1)/s)
        if model == "fraction" : 
            f = (s*P*(1-P)*(P-(2*s-1)/s))/(1-s+s*(1-P)**2)       
      
        P = la.spsolve(B_, B.dot(P) + dt*f)

        if t>=mod*nb_graph and graph_type != None :  
            graph_tanaka(X,P,t,model,graph_para)
            nb_graph += 1
            
        if np.isin(True, (1-P)>0.5) and np.isin(True, (1-P)<0.99) :             # wild pop is still present somewhere in the environment
            position = np.append(position,np.where((1-P)>0.5)[0][0])            # List containing the first position where the 
                                                                                # proportion of wild alleles is lower than 0.5.               
        if t%mod == 0 : 
            #print("\nTemps :", t) ; print("Speed :", np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt, "\n")
            speed_fct_of_time[:,int(t)//mod-1] = np.array([t,np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt])
            # first line : time, second line : speed of the wave at the corresponding time
                                                                                                                                                 
    if np.shape(position)[0] != 0 :
        position = position[position!=0]                          
        speed = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave     
    else :
        speed = None
        print(f"Can't determine the speed of the wave for s = {s}.")                                                 
            
    return(P, speed, speed_fct_of_time)
    
    
    
    
    
    
    
def graph_tanaka(X,P,t,model,graph_para):
    
        graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_figure, speed_proportion = graph_para

        fig, ax = plt.subplots()        
        ax.plot(X, P, label=f'Drive', color = "deeppink", linewidth=line_size)
        ax.set(xlabel='Space', ylabel='Proportion', ylim= (-0.03,1.03))
        if grid == True : 
            ax.grid()
        ax.xaxis.label.set_size(label_size)
        ax.yaxis.label.set_size(label_size)     
        if number_on_x_axe :
            ax.xaxis.set_tick_params(labelsize=number_x_size)
        else : 
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            ax.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax.yaxis.set_tick_params(labelsize=number_y_size)
        ax.yaxis.set_ticks(np.arange(0, 2, 1))
        plt.rc('legend', fontsize=legend_size)        
        ax.legend()
        plt.show() 