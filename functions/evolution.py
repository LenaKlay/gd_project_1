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

# External functions 
from graph import graph

# equ is the allele concerned by the equation.
# oth1 and oth2 are the two other alleles, indifferently.

# Growth function
def growth(equ, oth1, oth2, bio_para, model_para):
    # Parameters
    r,s,h,a,difW,difH,difD,c,homing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    # Function
    n = (equ+oth1+oth2)
    if growth_dynamic == "exponential":
        return(r+1)
    if growth_dynamic == "allee_effect" :
        return(np.max(r*(1-n)*(n-a)+1,0))
    if growth_dynamic == "logistical" and linear_growth == False :
        return(r*(1-n)+1)
    if growth_dynamic == "logistical" and linear_growth == True :  
        return(r*np.minimum(1-(equ+oth1+oth2),equ) + equ)         


# Death function    
def death(equ, oth1, oth2, bio_para, model_para):
    # Parameters
    r,s,h,a,difW,difH,difD,c,homing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    # Function
    n = (equ+oth1+oth2)
    if death_dynamic == "exponential" :
        return(1)
    if death_dynamic == "logistical" :
        return(1+r*n)
    if death_dynamic == "allee_effect" :
        return(r+1+r*(1-n)*(n-a))
    
    
# Mating function
def mating(W, H, D, bio_para, model_para):  
    # Parameters
    r,s,h,a,difW,difH,difD,c,homing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    
    # Function
    if linear_mating == False : 
        
        WW = (W**2)/(W+H+D)
        HH = (H**2)/(W+H+D)
        DD = (D**2)/(W+H+D)    
        WH = (W*H)/(W+H+D)   # *2
        WD = (W*D)/(W+H+D)   # *2 
        HD = (D*H)/(W+H+D)   # *2
    
        if homing == "zygote" :       
            mat1 = WW + WH + 0.25*HH 
            mat2 = (1-c)*(WH + 2*WD + 0.5*HH + HD) 
            mat3 = c*WH + 2*c*WD + (0.5*c+0.25)*HH + (1+c)*HD + DD  
            
        if homing == "germline":
            mat1 = WW + (1-c)*WH + 0.25*(1-c)**2*HH 
            mat2 = (1+c)*WH + 2*WD + 0.5*(1-c**2)*HH + (1-c)*HD
            mat3 = 0.25*(c+1)**2*HH + (1+c)*HD + DD
            
    if linear_mating == True and homing == "zygote" : 
        mat1 = (W>D)
        mat2 = 0
        mat3 = ((W>D)+1)  

    return([mat1,mat2,mat3])          
        
    
    
# Main evolution function      
def evolution(bio_para, model_para, num_para, graph_para, what_to_do) :  

    r,s,h,a,difW,difH,difD,c,homing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    T,L,M,N,mod,theta = num_para
    graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_figure, speed_proportion, show_graph_ini = graph_para
    
    # Saving and figures parameters
    directory = f"evolution/{homing}/r_{np.round(r,3)}/s_{np.round(s,3)}/h_{np.round(h,2)}_c_{np.round(c,2)}"
    file = "evo"

    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    X = np.linspace(0,N,N+1)*dx   
    
    # Initialization             
    if CI == "center" : 
        W = np.ones(N+1); W[0:N//2] = 0    # Wild individuals at t=0  
        H = np.zeros(N+1)                                               # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//2] = 1                # Drive individuals at t=0
    if CI == "left" : 
        W = np.ones(N+1); W[0:N//10] = 0    # Wild individuals at t=0  
        H = np.zeros(N+1)                                                # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//10] = 1                # Drive individuals at t=0
    
    if graph_type != None and show_graph_ini :
        graph(X,W,H,D,0,graph_para,bio_para,num_para,directory,file,"t = 0")
    nb_graph = 1
        
    increasing_WT_wave = True # indicate if the wave is monotone increasing 
    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.
    speed_fct_of_time = np.zeros((2,T//mod))   # first line : time, second line : speed of the wave at the corresponding time
    
    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    Bw = sp.identity(N+1)+((1-theta)*difW*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    Bh = sp.identity(N+1)+((1-theta)*difH*dt/dx**2)*A            # for W, H and D.
    Bd = sp.identity(N+1)+((1-theta)*difD*dt/dx**2)*A     

    Bw_ = sp.identity(N+1)-(theta*difW*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  
    Bh_ = sp.identity(N+1)-(theta*difH*dt/dx**2)*A               # for W, H and D.
    Bd_ = sp.identity(N+1)-(theta*difD*dt/dx**2)*A        

    # Evolution
    for t in np.linspace(dt,T,M) : 
        
        t = round(t,2)
        
        f1 = growth(W,H,D,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[0] - death(W,H,D,bio_para,model_para)*W
        f2 = (1-s*h)*growth(H,W,D,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[1] - death(W,H,D,bio_para,model_para)*H
        f3 = (1-s)*growth(D,H,W,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[2] - death(W,H,D,bio_para,model_para)*D
      
        W = la.spsolve(Bw_, Bw.dot(W) + dt*f1)
        H = la.spsolve(Bh_, Bh.dot(H) + dt*f2)
        D = la.spsolve(Bd_, Bd.dot(D) + dt*f3)
        
        # Graph
        if t>=mod*nb_graph and graph_type != None : 
            graph(X,W,H,D,t,graph_para,bio_para,num_para,directory,file,f"t = {np.round(t,2)}")
            nb_graph += 1
                            
        # Compute the speed on WT proportion or density
        if speed_proportion == True : WT = W/(W+H+D)
        else : WT = W
        # Check if the WT wave is partially decreasing
        epsilon = -0.000000001
        index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        # if the wave is increasing OR if the first increasing section come above 0.2 : the position will be correct. 
        if len(index_neg) == 0 or WT[index_neg[0]] > 0.2 :
            # if we realize only now that the wave is not increasing, it means that the previously recorded positions are false. We erase it and start again.
            if len(index_neg) != 0 and increasing_WT_wave :
                position = np.array([])   
                speed_fct_of_time = np.zeros((2,T//mod))   
                increasing_WT_wave = False
            # we recorde the position only if the WT wave is still in the environment
            if np.isin(True, WT>0.2) and np.isin(True, WT<0.99) :       
                # List containing the first position of WT where the proportion of wild alleles is higher than 0.2.  
                position = np.append(position,np.where(WT>0.2)[0][0])   
             
        if t%mod == 0 and what_to_do == "speed function of time" : 
            #print("\nTemps :", t) ; print("Speed :", np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt, "\n")
            speed_fct_of_time[:,int(t)//mod-1] = np.array([t,np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt])
            # first line : time, second line : speed of the wave at the corresponding time


    # Compute the wave speed with the position vector
    if np.shape(position)[0] != 0 :
        # Do not count positions when the wave was partially outside the window.
        position = position[position!=0]   
        # Speed compute on the last fifth of the position vector.
        speed = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave   
    else :
        speed = None
        print(f"Can't determine the speed of the wave : For r = {r} and s = {s}, the environment is empty or not mixed enought (at least one place where W<0.5 and one where W>0.1), or T is too short.")                                                 

    if what_to_do == "speed function of time" : 
        return(speed,W,H,D,speed_fct_of_time)
    if what_to_do == "system+evolution" : 
        return(speed,position,W,H,D)
    else :
        return(speed,W,H,D)




        
