#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""


# A checker avant de lancer : 
# h et c
# homing


num = 39

############################## Libraries ########################################

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
             
############################# Functions #########################################

# equ is the allele concerned by the equation.
# oth1 and oth2 are the two other alleles, indifferently.
    
def growth(equ,oth1,oth2,r,a):
    n = equ+oth1+oth2; N = len(n)
    if growth_dynamic == "constant":
        return(r+1)
    if growth_dynamic == "allee_effect" :
        return(np.maximum(r*(1-n)*(n-a)+1, np.zeros(N)))
    if growth_dynamic == "logistical" and linear_growth == False :
        return(r*(1-n)+1)
    if growth_dynamic == "logistical" and linear_growth == True :  
        return(r*np.minimum(1-(equ+oth1+oth2),equ) + equ)      
        

    
def death(equ,oth1,oth2,r,a):
    n = (equ+oth1+oth2)
    if death_dynamic == "constant" :
        return(1)
    if death_dynamic == "logistical" :
        return(1+r*n)
    if death_dynamic == "allee_effect" :
        return(r+1+r*(1-n)*(n-a))

           
      
def mating(W,H,D):
    
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
        
    
     

            
# Main evolution function (1D)     
def evolution(r,s,h,difWW,difDW,difDD,c,T,L,M,N,theta,conversion_timing) :  

    # Parameters initialization
    position = np.array([])   # list containing the first position where the proportion of wild alleles is higher than the threshold value.
    speed_fct_of_time = np.array([])      # speed computed... 
    time = np.array([])       # ...for each value of time in this vector.
    threshold = 0.2 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    pw_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    coex_density = -1  # made to be replace by the total population density in the coexistence state (if any)
    
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
      
    # Initialization             
    W = np.ones(N+1); W[0:N//2] = 0    # Wild individuals at t=0  
    H = np.zeros(N+1)                  # Heterozygous individuals at t=0
    D = np.zeros(N+1); D[0:N//2] = 1   # Drive individuals at t=0

        
    # Matrix
    C0 = -2*np.ones(N-1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N-1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N-1, N-1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    Bw = sp.identity(N-1)+((1-theta)*difWW*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    Bh = sp.identity(N-1)+((1-theta)*difDW*dt/dx**2)*A            # for W, H and D.
    Bd = sp.identity(N-1)+((1-theta)*difDD*dt/dx**2)*A     

    Bw_ = sp.identity(N-1)-(theta*difWW*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  
    Bh_ = sp.identity(N-1)-(theta*difDW*dt/dx**2)*A               # for W, H and D.
    Bd_ = sp.identity(N-1)-(theta*difDD*dt/dx**2)*A        

    # Evolution
    for t in np.linspace(dt,T,M) :
        
        t = round(t,2)
        
        f1 = growth(W,H,D,r,a)*mating(W,H,D)[0] - death(W,H,D,r,a)*W
        f2 = (1-s*h)*growth(H,W,D,r,a)*mating(W,H,D)[1] - death(W,H,D,r,a)*H
        f3 = (1-s)*growth(D,H,W,r,a)*mating(W,H,D)[2] - death(W,H,D,r,a)*D      
      
        W[1:-1] = la.spsolve(Bw_, Bw.dot(W[1:-1]) + dt*f1[1:-1])
        H[1:-1] = la.spsolve(Bh_, Bh.dot(H[1:-1]) + dt*f2[1:-1])
        D[1:-1] = la.spsolve(Bd_, Bd.dot(D[1:-1]) + dt*f3[1:-1])
             
        W[0] = W[1]; H[0] = H[1]; D[0] = D[1]   # Neumann condition alpha=0
        W[-1] = W[-2]; H[-1] = H[-2]; D[-1] = D[-2]  # Neumann condition beta=0
                                
        # Compute the speed on WT allele proportion or density
        if conversion_timing == "zygote" : 
            WT = (W+0.5*H)/(W+H+D) 
        if conversion_timing == "germline" : 
            WT = (W+0.5*(1-c)*H)/(W+H+D) 
        
        # Check if the WT wave is strictly increasing
        epsilon = -0.01 ; index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        if len(index_neg) != 0 : print("WT proportion decreasing : problem !") 
               
        # Check if there is a coexistence state (if so adapt threshold and erase all speed values)
        center_index = np.arange(N//2-N//10,N//2+N//10)
        if np.all((WT[center_index]>0.05) & (WT[center_index]<0.95)) and pw_star < 0 :    # if there exists a palliate in ]0,1[ around N//2 (and pw_star as not been changed yet)
            pw_star = np.mean(WT[center_index])                                            # if this palliate is above the threshold value
            threshold = (pw_star+1)/2
            coex_density = (W+H+D)[N//2]
            position = np.array([])   
            speed_fct_of_time = np.array([])   
            time = np.array([]) 
        # we recorde the position only if the WT wave is still in the environment. We do not recorde the 0 position since the threshold value of the wave might be outside the window.
        if np.isin(True, WT>threshold) and np.isin(True, WT<0.99) and np.where(WT>threshold)[0][0] != 0:  
            # List containing the first position of WT where the proportion of wild alleles is higher than the threshold.  
            position = np.append(position, np.where(WT>threshold)[0][0])
            # Wait a minimum time and a minimum number of positions to compute the speed
            if len(position) > 20  : 
                # The speed is computed on the last fifth of the position vector.
                time = np.append(time,t)
                speed_fct_of_time = np.append(speed_fct_of_time, np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt)
        # if the threshold value of the wave is outside the window, stop the simulation  
        if not(np.isin(False, WT>threshold) and np.isin(False, WT<threshold) ) :
            print("t =",t)
            break

    # return the speed, and 1 if there was a coexistence state (0 otherwise)
    return(speed_fct_of_time[-1], coex_density)





def heatmap(precision,hmin,hmax,smin,smax):

    # Range for s and h
    delta_h = (hmax-hmin)/precision    
    h_range = np.arange(hmin+delta_h/2,hmax+delta_h/2,delta_h)   
    delta_s = (smax-smin)/precision    
    s_range = np.arange(smin+delta_s/2,smax+delta_s/2,delta_s)   
        
    heatmap_vect = np.zeros(precision) 
    coex_vect = np.zeros(precision) 
 
    s_index = num-1
    s = s_range[s_index]
     
    for h_index in range(0, precision) :          
            h = h_range[h_index]
            print("\nh=", np.round(h,3)) 
            print("c=", np.round(c,3)) 
            print("r=", np.round(r,3)) 
            print("s=", np.round(s,3)) 
            heatmap_vect[h_index], coex_vect[h_index] = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta, homing)
            print("speed=", heatmap_vect[h_index]) 
            print("coex=", coex_vect[h_index]) 

    return(heatmap_vect, coex_vect) 
  




############################### Parameters ######################################


# Heatmap

precision = 50
  
# Biological
    
c = 0.85     # homing rate
r = 10
homing = "germline"   # "zygote" or "germline"

a = 0.2      # coefficient for allee effect (growth or death)
difW = 1   # diffusion coefficient for WW individuals
difH = 1   # diffusion coefficientrate for WD individuals
difD = 1   # diffusion coefficient rate for DD individuals

# How Population Grow and Decrease

growth_dynamic = "constant"      # exponential or logistical or allee_effect
death_dynamic = "logistical"      # exponential or logistical or allee_effect


# Linearization
  
linear_growth = False   
linear_mating = False          

theta = 0.5     # discretization in space : theta = 0.5 for Crank Nicholson
                # theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  

                 


############################### Results #########################################

T = 500; L = 2000; M = T*40; N = L*10
hmin = 0; hmax = 1 
smin = 0; smax = 1       
heatmap_vect, coex_vect = heatmap(precision,hmin,hmax,smin,smax)
np.savetxt(f'heatmap_{num}.txt', heatmap_vect)  
np.savetxt(f'coex_{num}.txt', coex_vect)  

        
    
    




