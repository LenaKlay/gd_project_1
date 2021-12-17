#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""


# A checker avant de lancer : 
# heatmap_type
# h et c
# homing

num = 10

############################## Libraries ########################################

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
             
############################# Functions #########################################

# equ is the allele concerned by the equation.
# oth1 and oth2 are the two other alleles, indifferently.
    
def growth(equ,oth1,oth2,r,a):
    n = (equ+oth1+oth2)
    if growth_dynamic == "exponential":
        return(r+1)
    if growth_dynamic == "allee_effect" :
        return(np.max(r*(1-n)*(n-a)+1,0))
    if growth_dynamic == "logistical" and linear_growth == False :
        return(r*(1-n)+1)
    if growth_dynamic == "logistical" and linear_growth == True :  
        return(r*np.minimum(1-(equ+oth1+oth2),equ) + equ)         
    
def death(equ,oth1,oth2,r,a):
    n = (equ+oth1+oth2)
    if death_dynamic == "exponential" :
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
        
    
                        
        
def evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta): 
    
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Initialization 
            
    if CI == "center" : 
        W = np.ones(N+1); W[0:N//2] = 0          # Wild individuals at t=0  
        H = np.zeros(N+1)                        # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//2] = 1         # Drive individuals at t=0
    if CI == "left" : 
        W = np.ones(N+1); W[0:N//10] = 0         # Wild individuals at t=0  
        H = np.zeros(N+1)                        # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//10] = 1        # Drive individuals at t=0
        
    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.

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
        
        f1 = growth(W,H,D,r,a)*mating(W,H,D)[0] - death(W,H,D,r,a)*W
        f2 = (1-s*h)*growth(H,W,D,r,a)*mating(W,H,D)[1] - death(W,H,D,r,a)*H
        f3 = (1-s)*growth(D,H,W,r,a)*mating(W,H,D)[2] - death(W,H,D,r,a)*D
      
        W = la.spsolve(Bw_, Bw.dot(W) + dt*f1)
        H = la.spsolve(Bh_, Bh.dot(H) + dt*f2)
        D = la.spsolve(Bd_, Bd.dot(D) + dt*f3)
        
        if np.isin(True, W/(W+H+D)>0.5) and np.isin(True, W/(W+H+D)<0.99) :  # wild pop is still present somewhere in the environment
            position = np.append(position,np.where(W/(W+H+D)>0.5)[0][0])        # List containing the first position where the 
                                                                                    # proportion of wild alleles is lower than 0.5.   
                                                                                    
    if np.shape(position)[0] != 0 :
        position = position[position!=0]                          
        speed = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave                                                   

    return(speed,W,H,D)



def tanaka(s,T,L,M,N,theta,model): 
    
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Initialization       
    if CI == "center" : 
        P = np.zeros(N+1); P[0:N//2] = 1            # Proportion of drive individuals at t=0
    if CI == "left" :      
        P = np.zeros(N+1); P[0:N//10] = 1           # Proportion of drive individuals at t=0

    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    B = sp.identity(N+1)+((1-theta)*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    B_ = sp.identity(N+1)-(theta*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  

    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.

    # Evolution
    for t in np.linspace(dt,T,M) : 
        
        t = round(t,2)
        if model == "cubic" : 
            f = s*P*(1-P)*(P-(2*s-1)/s)
        if model == "fraction" : 
            f = (s*P*(1-P)*(P-(2*s-1)/s))/(1-s+s*(1-P)**2)       
      
        P = la.spsolve(B_, B.dot(P) + dt*f)

        if np.isin(True, (1-P)>0.5) and np.isin(True, (1-P)<0.99) :             # wild pop is still present somewhere in the environment
            position = np.append(position,np.where((1-P)>0.5)[0][0])            # List containing the first position where the 
                                                                                # proportion of wild alleles is lower than 0.5.               
                                                                                                                                                
    if np.shape(position)[0] != 0 :
        position = position[position!=0]                          
        speed = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave        
    return(P, speed)





def heatmap(precision,smin,smax,rmin,rmax, what):
    s_range = np.linspace(smin,smax,precision)             # range of values for s
    r_range = np.linspace(rmin,rmax,precision)[num-1]             # range of values for r
    heatmap_values = np.zeros(precision)          # matrix containing all the wave speeds for the different values of s and r 
    zero_done = False
    zero_line = np.array([])
    
    for s_index in range(0, precision) :
            r = r_range
            s = s_range[s_index]
 
            if what == "classic" :
                heatmap_values[s_index] = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta)[0]
                if s_index != 0 and heatmap_values[s_index-1]*heatmap_values[s_index]<=0 and heatmap_values[s_index-1] != 0 and zero_done == False :
                    zero_line = np.array([s_index-0.5])
                    zero_done = True 
                    
            if what == "speed_cubic" : 
                heatmap_values[s_index] = np.abs(evolution(r,s,1,1,1,1,1,T,L,M,N,0.5)[0]-tanaka(s,T,L,M,N,0.5,"cubic")[1])
            if what == "speed_fraction" : 
                heatmap_values[s_index] = np.abs(evolution(r,s,1,1,1,1,1,T,L,M,N,0.5)[0]-tanaka(s,T,L,M,N,0.5,"fraction")[1])
            if what == "r_one_minus_n_cubic" : 
                speed_girardin, W, H, D = evolution(r,s,1,1,1,1,1,T,L,M,N,0.5)
                n = D+H+W
                heatmap_values[s_index] = np.max(abs(r*(1-n)))
            if what == "r_one_minus_n_fraction" :
                speed_girardin, W, H, D = evolution(r,s,1,1,1,1,1,T,L,M,N,0.5)
                p = D/(D+H+W); n = D+H+W
                heatmap_values[s_index] = np.max(abs(r*(1-n) - s*p*(2-p)/(1-s+s*(1-p)**2))) 
                   
    if what == "classic" and len(zero_line) != 0 : 
        return(heatmap_values, zero_line) 
    else : return(heatmap_values) 
    



############################### Parameters ######################################

# Biological
    
c = 1/4     # homing rate
h = 1/4     # and sh for heterozygous individuals
a = 0       # coefficient for allee effect (growth or death)

homing = "zygote"   # "zygote" or "germline"

difW = 1   # diffusion coefficient for WW individuals
difH = 1   # diffusion coefficientrate for WD individuals
difD = 1   # diffusion coefficient rate for DD individuals


# Initialization

CI = "center"       # "center" for having the border in the center, "left" for having the border on the left

# How Population Grow and Decrease

growth_dynamic = "logistical"     # exponential or logistical or allee_effect
death_dynamic = "exponential"      # exponential or logistical or allee_effect


# Linearization
  
linear_growth = False   
linear_mating = False          

theta = 0.5     # discretization in space : theta = 0.5 for Crank Nicholson
                # theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  

# Heatmap

heatmap_type = "classic"   #  "classic"   "speed_cubic"   "speed_fraction"   "r_one_minus_n_cubic"   "r_one_minus_n_fraction"                              
precision = 30
                   


############################### Results #########################################

what_to_do = "heatmap" 
       
if what_to_do == "heatmap" :
    CI = "center"
    theta = 0.5                

    if heatmap_type == "classic" :             
        T = 500; L = 2000; M = T*6; N = L 
        smin = 0.3; smax = 0.9; rmin = 0 ; rmax = 12 
        heatmap_values, zero_line = heatmap(precision,smin,smax,rmin,rmax, "classic")
        np.savetxt(f'{heatmap_type}_zero_line_{num}.txt', zero_line)     
    else :  
        print("yes")
        T = 1000; L = 4000; M = T*40; N = L 
        smin = 0.3; smax = 0.9; rmin = 50 ; rmax = 60 
        linear_growth = False ; linear_mating = False              
        homing = "zygote" 
        heatmap_values = heatmap(precision,smin,smax,rmin,rmax,heatmap_type)
    np.savetxt(f'{heatmap_type}_heatmap_{num}.txt', heatmap_values)  

        
    
    




