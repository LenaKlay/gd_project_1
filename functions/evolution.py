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
    graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_figure, WT_proportion_wave, show_graph_ini = graph_para
      
    # Check for a coexitence stable state 
    if homing == "zygote" and (1-h)*(1-c) > 0.5 : 
        p_star = (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h)))
        #mean_fitness = (1-s)*p**2+2*(c*(1-s)+(1-c)*(1-s*h))*p*(1-p)+(1-p)**2                
    if homing == "germline" and h < 0.5 : 
        p_star = ((1-s*h)*(1+c)-1)/(s*(1-2*h))
        #mean_fitness = (1-s)*p**2+2*(1-s*h)*p*(1-p)+(1-p)**2 
    # If an admissible coexistence stable state exists, we change the treshold value.
    if  p_star > 0 and p_star < 1 : 
        #n_star =  1-(1-mean_fitness)/(mean_fitness*r)
        WT_allele_wave = True
        treshold = (1-p_star) + 0.001
        T = 5000; L = 5000; M = T*6; N = L 
             
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    X = np.linspace(0,N,N+1)*dx   
    
    # Initialization             
    if CI == "center" : 
        W = np.ones(N+1); W[0:N//2] = 0    # Wild individuals at t=0  
        H = np.zeros(N+1)                  # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//2] = 1   # Drive individuals at t=0
    if CI == "left" : 
        W = np.ones(N+1); W[0:N//10] = 0    # Wild individuals at t=0  
        H = np.zeros(N+1)                   # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//10] = 1   # Drive individuals at t=0
    
    if graph_type != None and show_graph_ini :
        graph(X,W,H,D,0,graph_para,bio_para,num_para)
    nb_graph = 1
        
    increasing_WT_wave = True # indicates if the WT wave is monotone increasing (can be decreasing for a short transition phase at the beginning)
    p_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    treshold = 0.2 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    WT_allele_wave = False
    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.
    speed = np.array([])    # speed computed at each step of time : 
    #speed_fct_of_time = np.zeros((2,T//mod))   # to plot the speed function of time, with mod steps (first line : time, second line : speed of the wave at the corresponding time)
    
    
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
            graph(X,W,H,D,t,graph_para,bio_para,num_para)
            nb_graph += 1
                                          
        # Compute the speed on WT proportion or density
        if WT_allele_wave : 
            if homing == "zygote" : WT = (W+0.5*H)/(W+H+D) 
            if homing == "germline" : WT = (W+0.5*(1-c)*H)/(W+H+D) 
        if WT_proportion_wave : 
            WT = W/(W+H+D)
        else : WT = W
        
        # If the WT wave is partially decreasing, index_neg is not empty.
        epsilon = -0.0001 ; index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        
        # if the wave is increasing OR if the first increasing section come above the treshold : the position will be correct.       
        if len(index_neg) == 0 or WT[index_neg[0]] > treshold :
            # if we realize only now that the wave is not increasing, it means that the previously recorded positions are false. We erase it and start again. index_neg[0]>3 is there to assure 
            # that the decreasing is not a numerical error due to borders. 
            if len(index_neg) != 0 and increasing_WT_wave and index_neg[0]>3 :
                position = np.array([])   
                speed = np.array([])   
                #speed_fct_of_time = np.zeros((2,T//mod))   
                increasing_WT_wave = False
            # we recorde the position only if the WT wave is still in the environment. We do not recorde the 0 position since the treshold value of the wave might be outside the window.
            if np.isin(True, WT>treshold) and np.isin(True, WT<0.99) and np.where(WT>treshold)[0][0] != 0:  
                # List containing the first position of WT where the proportion of wild alleles is higher than the treshold.  
                position = np.append(position, np.where(WT>treshold)[0][0])
                # Wait a minimum time and a minimum number of positions to compute the speed
                if np.shape(position)[0] > 5 : 
                    # The speed is computed on the last fifth of the position vector.
                    speed = np.append(speed, np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt)
                    # Stop the simulation if the speed stabilies
                    #if np.shape(speed)[0] > 500 and np.mean(np.abs((speed[1:]-speed[:-1])[-501:-1]))<0.00001 : 
                    #    print("ça y est !")
                    #    print(position)
                    #    print(speed)
                    #    print("Time :",t)
                    #    break
                    
        # NB : problems encountered when computing the speed : 
        # - decreasing section   ->     r=1.08 s=0.668 c=h=0.5 homing="germline"  T=1000  L=5000  M=6000  N=5000 
        # - small perturbation at the border, which can create fake decreasing sections 
        # - coexistance with WT equilibrium under the 0.2 treshold value   ->    r=0.36 s=0.332 c=h=0.25 homing="zygote"  T=5000  L=10000  M=T*6  N=L 


        # Speed fonction of time    
        #if t%mod == 0 and what_to_do == "speed function of time" : 
            # first line : time, second line : speed of the wave at the corresponding time
        #    speed_fct_of_time[:,int(t)//mod-1] = np.array([t,np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt])
            
        # if the treshold value of the wave is outside the window, stop the simulation  
        if not(np.any(WT>treshold) and np.any(WT<treshold)) : 
            break
        
    # Compute the wave speed with the position vector
    if np.shape(position)[0] == 0 :
        # Speed compute on the last fifth of the position vector.
        #speed[-1] = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave   
    #else :
        speed[-1] = None
        print(f"Can't determine the speed of the wave for r = {r} and s = {s} : the position vector is empty.")                                                 

    if what_to_do == "speed function of time" : 
        return(speed[-1],W,H,D,speed)
    if what_to_do == "system+evolution" : 
        return(speed[-1],position,W,H,D)
    else :
        return(speed[-1],W,H,D)




        
