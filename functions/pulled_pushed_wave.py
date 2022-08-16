#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:03:42 2021

@author: lena
"""

# Libraries
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg  as la
import matplotlib.pyplot as plt

# External functions 
from graph import save_fig_or_data
from evolution import evolution 

# Change font to serif
plt.rcParams.update({'font.family':'serif'})

# Graph parameters 
title_size = 15
label_size = 17
legend_size = 12
line_size = 0.7
number_on_x_axe = False
number_x_size = 10
number_y_size = 20

# Parameters used for the figure : 
#r = 2              # growth rate
#s = 0.2            # when selection only acts on survival with a fitness cost s (b=1 and d=1) 
#c = 1              # homing rate
#homing = "zygote"  # "zygote" or "germline"
#CI = "left"        # "center" for having the border in the center, "left" for having the border on the left
#T = 50             # final time
#L = 200            # length of the spatial domain
#M = T*6            # number of time steps
#N = L              # number of spatial steps




def pulled_pushed_graph(X,W,D1,D2,t,bio_para,num_para):
       
        r,s,h,a,difW,difH,difD,c,homing = bio_para
        T,L,M,N,mod,theta = num_para
        
        fig, ax = plt.subplots()       
        plt.fill_between(X, D1, color="#afe88e", label="back",linestyle="None")
        plt.fill_between(X, D1, D1+D2, color="#509369", label="front",linestyle="None")
        plt.plot(X,D1+D2, color="black", linewidth=line_size)                 
                         
        # Graphic size, title and labels
        ax.set(xlabel='Space', ylabel="population densities", ylim = (0,D1[0]+0.1*D1[0]), xlim=(X[0],X[-1]))
        #ax.set_title(f"t = {t}", fontsize = title_size, loc='right')
            
        # Grid
        ax.grid()
            
        # Labels and legend sizes
        ax.xaxis.label.set_size(label_size)
        ax.yaxis.label.set_size(label_size)     
        # No axis values
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        plt.rc('legend', fontsize=legend_size)  
        ax.legend(bbox_to_anchor=(0.553,1.13), ncol=2)
        
        # Show the graph      
        plt.show()
        
        # Saving figures and datas
        directory = f"pulled_pushed/{homing}_r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"
        
        num = str(int(t)//mod)
        if len(num)==1: num = '0'+'0'+num
        if len(num)==2: num = '0'+num
        save_fig_or_data(directory, fig, [], f"{num}", bio_para, num_para)
        #save_fig_or_data(directory, fig, [], f"t_{t}", bio_para, num_para)
        #columns = [X,W,D]; np.savetxt(f"../outputs/{directory}/t_{t}.txt", np.column_stack(columns), fmt='%.3e', delimiter="  ") 
        
    


    
# Main evolution function      
def pulled_pushed_wave(bio_para, model_para, num_para, graph_para, what_to_do) :  
    
    graph_para[0] = None; graph_para[-1] = False; graph_para[-2] = False  #; graph_para[6] = False
    W, H, D, time, speed = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
    
    r,s,h,a,difW,difH,difD,c,homing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    T,L,M,N,mod,theta = num_para
    graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, WT_proportion_wave, show_graph_ini, show_graph_fin = graph_para
   
    # Two sub-population of drive
    D1 = np.zeros(len(D))
    D2 = np.zeros(len(D))
    # Very front of the wave index
    wave_front = np.where(D<0.01)[0][0]
    # Cut the wave in two parts, and replace it on the left of the domain
    
    a=N//40  # size of the front part
    b=3  # N//b where to replace the wave
    D1[0:N//b] = D[wave_front-a-N//b:wave_front-a]
    D2[N//b:N//b+a] = D[wave_front-a:wave_front] 
                              
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    X = np.linspace(0,N,N+1)*dx   
    pulled_pushed_graph(X,W,D1,D2,0,bio_para,num_para)
    nb_graph = 1

    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    Bw = sp.identity(N+1)+((1-theta)*difW*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    Bh = sp.identity(N+1)+((1-theta)*difD*dt/dx**2)*A            # for W, H and D.
    Bd = sp.identity(N+1)+((1-theta)*difD*dt/dx**2)*A     

    Bw_ = sp.identity(N+1)-(theta*difW*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  
    Bh_ = sp.identity(N+1)-(theta*difD*dt/dx**2)*A               # for W, H and D.
    Bd_ = sp.identity(N+1)-(theta*difD*dt/dx**2)*A        

    # Evolution
    for t in np.linspace(dt,T,M) :
        
        t = round(t,2)
        
        f1 = (r*(1-(W+D1+D2))+1) * (W**2)/(D1+D2+W) - W
        f2 = (1-s)*(r*(1-(W+D1+D2))+1) * (D1**2+2*D1*W+D1*D2)/(D1+D2+W) - D1
        f3 = (1-s)*(r*(1-(W+D1+D2))+1) * (D2**2+2*D2*W+D1*D2)/(D1+D2+W) - D2
        
        W = la.spsolve(Bw_, Bw.dot(W) + dt*f1)
        D1 = la.spsolve(Bh_, Bh.dot(D1) + dt*f2)
        D2 = la.spsolve(Bd_, Bd.dot(D2) + dt*f3)
        
        # Graph
        if t>=mod*nb_graph : 
            pulled_pushed_graph(X,W,D1,D2,t,bio_para,num_para)
            nb_graph += 1

        # if the treshold value of the wave is outside the window, stop the simulation  
        if not(np.isin(False, W>0.5) and np.isin(False, W<0.5) ) :
            print("t =",t)
            break
                    
        # NB : problems encountered when computing the speed : 
        # - decreasing section   ->     r=1.08 s=0.668 c=h=0.5 homing="germline"  T=1000  L=5000  M=6000  N=5000 
        # - small perturbation at the border, which can create fake decreasing sections 
        # - coexistance with WT equilibrium under the 0.2 treshold value   ->    r=0.36 s=0.332 c=h=0.25 homing="zygote"  T=5000  L=10000  M=T*6  N=L 

    return(W,D1,D2)



        
