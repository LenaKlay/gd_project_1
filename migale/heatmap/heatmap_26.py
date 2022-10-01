#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""


# A checker avant de lancer : 
# h et c
# homing
# rlog

num = 26

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
        
    # Parameters initialization
    increasing_WT_wave = True # indicates if the WT wave is monotone increasing (can be decreasing for a short transition phase at the beginning)
    p_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    treshold = 0.2 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    WT_allele_wave = False
    position = np.array([])   # list containing the first position where the proportion of wild alleles is higher than the treshold value.
    speed_fct_of_time = np.array([])      # speed computed... 
    time = np.array([])       # ...for each value of time in this vector.
 
    # Check for a coexitence stable state 
    if (homing == "zygote" and (1-h)*(1-c) > 0.5) or (homing == "germline" and h < 0.5) : 
            #s_1 = c/(1-h*(1-c))
            if homing == "zygote" : 
                #s_2 = c/(2*c + h*(1-c))
                p_star = (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h)))
                #mean_fitness = (1-s)*p**2+2*(c*(1-s)+(1-c)*(1-s*h))*p*(1-p)+(1-p)**2 
            if homing == "germline" : 
                #s_2 = c/(2*c*h + h*(1-c))
                p_star = ((1-s*h)*(1+c)-1)/(s*(1-2*h))
                #mean_fitness = (1-s)*p**2+2*(1-s*h)*p*(1-p)+(1-p)**2    
            # If an admissible coexistence stable state exists, we change the treshold value.
            if p_star > 0 and p_star < 1 : 
                #n_star =  1-(1-mean_fitness)/(mean_fitness*r)
                WT_allele_wave = True
                treshold = ((1-p_star)+1)/2              
    
             
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial

    
    # Initialization   
    W = np.ones(N+1); W[0:N//2] = 0    # Wild individuals at t=0  
    H = np.zeros(N+1)                  # Heterozygous individuals at t=0
    D = np.zeros(N+1); D[0:N//2] = 1   # Drive individuals at t=0

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
                         
        # Compute the speed on WT proportion or density
        if WT_allele_wave : 
            if homing == "zygote" : WT = (W+0.5*H)/(W+H+D) 
            if homing == "germline" : WT = (W+0.5*(1-c)*H)/(W+H+D) 
        else : WT = W
        
        # If the WT wave is partially decreasing, index_neg is not empty.
        epsilon = -0.0001 ; index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        
        # if the wave is increasing OR if the first increasing section come above the treshold : the position will be correct.       
        if len(index_neg) == 0 or WT[index_neg[0]] > treshold :
            # if we realize only now that the wave is not increasing, it means that the previously recorded positions are false. We erase it and start again. index_neg[0]>3 is there to assure 
            # that the decreasing is not a numerical error due to borders. 
            if len(index_neg) != 0 and increasing_WT_wave and index_neg[0]>3 :
                position = np.array([])   
                speed_fct_of_time = np.array([])   
                time = np.array([])   
                increasing_WT_wave = False
            # we recorde the position only if the WT wave is still in the environment. We do not recorde the 0 position since the treshold value of the wave might be outside the window.
            if np.isin(True, WT>treshold) and np.isin(True, WT<0.99) and np.where(WT>treshold)[0][0] != 0:  
                # List containing the first position of WT where the proportion of wild alleles is higher than the treshold.  
                position = np.append(position, np.where(WT>treshold)[0][0])
                # Wait a minimum time and a minimum number of positions to compute the speed
                if len(position) > 20  : 
                    # The speed is computed on the last fifth of the position vector.
                    time = np.append(time,t)
                    speed_fct_of_time = np.append(speed_fct_of_time, np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt)
            # if the treshold value of the wave is outside the window, stop the simulation  
            if not(np.isin(False, WT>treshold) and np.isin(False, WT<treshold) ) :
                print("t =",t)
                break

    return(speed_fct_of_time[-1])



def heatmap(precision,smin,smax,rmin,rmax):

    # Range for r and s
    delta_s = (smax-smin)/precision    # delta_s is the size of a simulation pixel (size mesured with the s scale)  
    s_range = np.arange(smin+delta_s/2,smax+delta_s/2,delta_s)       # s values for simulations (NB : smin and smax are not simulated, we simulate values centered on the simulation pixels)                  
    if rlog : 
        r_range = np.logspace(-2, 1, num=precision)
    else : 
        delta_r = (rmax-rmin)/precision    # delta_r is the size of a simulation pixel (size mesured with the r scale)  
        r_range = np.arange(rmin+delta_r/2,rmax+delta_r/2,delta_r)       # r values for simulations (NB : rmin and rmax are not simulated, we simulate values centered on the simulation pixels)

    heatmap_vect = np.zeros(precision)    
 
    r_index = num-1
    r = r_range[r_index]
     
    for s_index in range(0, precision) :          
            s = s_range[s_index]
            print("\ns=", np.round(s,3))           
            heatmap_vect[s_index] = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta)

    return(heatmap_vect) 
  



############################### Parameters ######################################


# Heatmap

precision = 50
  
# Biological
    
c = 0.25     # homing rate
h = 0.75     # and sh for heterozygous individuals
rlog = True           # scale for r (normal or log)
homing = "germline"   # "zygote" or "germline"

a = 0       # coefficient for allee effect (growth or death)
difW = 1   # diffusion coefficient for WW individuals
difH = 1   # diffusion coefficientrate for WD individuals
difD = 1   # diffusion coefficient rate for DD individuals

# How Population Grow and Decrease

growth_dynamic = "logistical"     # exponential or logistical or allee_effect
death_dynamic = "exponential"      # exponential or logistical or allee_effect


# Linearization
  
linear_growth = False   
linear_mating = False          

theta = 0.5     # discretization in space : theta = 0.5 for Crank Nicholson
                # theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  

                 


############################### Results #########################################

T = 500; L = 2000; M = T*6; N = L 
smin = 0.1; smax = 0.9; 
if rlog: rmin = 0.01 ; rmax = 10    
else: rmin = 0 ; rmax = 12            
heatmap_vect = heatmap(precision,smin,smax,rmin,rmax)
np.savetxt(f'heatmap_{num}.txt', heatmap_vect)  

        
    
    




