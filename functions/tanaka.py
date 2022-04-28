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

# External functions
from graph import save_fig_or_data

def tanaka(s,model,model_para,num_para,graph_para): 
    T,L,M,N,mod,theta = num_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion, show_graph_ini, show_graph_fin = graph_para

    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    X = np.linspace(0,N,N+1)*dx   
    #save_fig_or_data(f"tanaka/{model}/s_{s}", [], X, "abscisse_X", None, num_para)
    
    # Initialization       
    if CI == "center" : 
        P = np.zeros(N+1); P[0:N//2] = 1               # Proportion of drive individuals at t=0
    if CI == "left" :      
        P = np.zeros(N+1); P[0:N//10] = 1              # Proportion of drive individuals at t=0
    
    if graph_type != None and show_graph_ini :
        graph_tanaka(s,X,P,0,model,graph_para,num_para)
    nb_graph = 1

    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    B = sp.identity(N+1)+((1-theta)*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    B_ = sp.identity(N+1)-(theta*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  

    treshold = 0.5            # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    position = np.array([])   # list containing the first position where the proportion of wild alleles is higher than the treshold value.
    speed_fct_of_time = np.array([])      # speed computed... 
    time = np.array([])                   # ...for each value of time in this vector.
 
    # Evolution
    for t in np.linspace(dt,T,M) : 
        
        t = round(t,2)
        if model == "cubic" : 
            f = s*P*(1-P)*(P-(2*s-1)/s)
        if model == "fraction" : 
            f = (s*P*(1-P)*(P-(2*s-1)/s))/(1-s+s*(1-P)**2)  
        if model == "KPP" : 
            f = P*(1-P)*(1-2*s)
      
        P = la.spsolve(B_, B.dot(P) + dt*f)

        if t>=mod*nb_graph and graph_type != None :
            graph_tanaka(s,X,P,t,model,graph_para,num_para)
            nb_graph += 1
            
        if np.isin(True, (1-P)>treshold) and np.isin(True, (1-P)<0.99) and np.where((1-P)>treshold)[0][0] != 0 :             
            # first position where the Wild-type wave is over the treshold value
            position = np.append(position,np.where((1-P)>treshold)[0][0])  
        elif np.isin(True, (1-P)<treshold) and np.isin(True, (1-P)>0.01) and np.where((1-P)<treshold)[0][0] != 0 :             
            # first position where the Wild-type wave is below the treshold value
            position = np.append(position,np.where((1-P)<treshold)[0][0])
                
        # compute the speed
        if len(position) > 20 : 
            time = np.append(time, t)
            speed_fct_of_time = np.append(speed_fct_of_time, np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt)
        # if the treshold value of the wave is outside the window, stop the simulation  
        if not(np.isin(False, (1-P)>treshold) and np.isin(False, (1-P)<treshold) ) :
            print("t =",t)
            break
        
    # last graph
    if show_graph_fin :   
        graph_tanaka(s,X,P,t,model,graph_para,num_para)
        
    # plot the speed function of time    
    if len(speed_fct_of_time) != 0 : 
        if graph_type!= None :
            fig, ax = plt.subplots()
            ax.plot(time, speed_fct_of_time) 
            ax.set(xlabel='Time', ylabel='Speed', title = f'Speed function of time')   
            if save_fig :
                save_fig_or_data(f"tanaka/{model}/s_{s}", fig, speed_fct_of_time, "speed_fct_of_time", None, num_para)
            plt.show() 
    else :
        print('No wave')     
        
    return(P, time, speed_fct_of_time)
    
    
    
    
    
    
def graph_tanaka(s,X,P,t,model,graph_para,num_para):
    
        graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion, show_graph_ini, show_graph_fin = graph_para

        fig, ax = plt.subplots()        
        ax.plot(X, P, label=f'Drive', color = "deeppink", linewidth=line_size)
        ax.set(xlabel='Space', ylabel='Proportion', ylim= (-0.03,1.03), title=f"Tanaka {model}, t={t}, s={s}")
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
        
        # Saving figures and data
        if save_fig : 
            save_fig_or_data(f"tanaka/{model}/s_{s}", fig, [], f"t_{t}", None, num_para)
        
        
        
        
        
        
        
        