#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""

num = 20

############################## Libraries ########################################

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
             
############################# Functions #########################################

# equ is the allele concerned by the equation.
# oth1 and oth2 are the two other alleles, indifferently.
    
def growth(equ,oth1,oth2,r,a):    
    n = equ+oth1+oth2
    if growth_dynamic == "constant":
        return(r+1)
    if growth_dynamic == "allee_effect" :
        #return(r*(1-n)*(n-a)+1)
        return(np.maximum(r*(1-n)*(n-a)+1, np.zeros(len(n))))
    if growth_dynamic == "logistical" :
        return(r*(1-n)+1)  
    
def death(equ,oth1,oth2,r,a):
    n = (equ+oth1+oth2)
    if death_dynamic == "constant" :
        return(1)
    if death_dynamic == "logistical" :
        return(1+r*n)
    if death_dynamic == "allee_effect" :
        return(r+1+r*(n-1)*(n-a))
        
      
def mating(W,H,D):
     
    WW = (W**2)/(W+H+D)
    HH = (H**2)/(W+H+D)
    DD = (D**2)/(W+H+D)    
    WH = (W*H)/(W+H+D)   # *2
    WD = (W*D)/(W+H+D)   # *2 
    HD = (D*H)/(W+H+D)   # *2

    if conversion_timing == "zygote" :       
        mat1 = WW + WH + 0.25*HH 
        mat2 = (1-c)*(WH + 2*WD + 0.5*HH + HD) 
        mat3 = c*WH + 2*c*WD + (0.5*c+0.25)*HH + (1+c)*HD + DD  
            
    if conversion_timing == "germline":
        mat1 = WW + (1-c)*WH + 0.25*(1-c)**2*HH 
        mat2 = (1+c)*WH + 2*WD + 0.5*(1-c**2)*HH + (1-c)*HD
        mat3 = 0.25*(c+1)**2*HH + (1+c)*HD + DD
            
    return([mat1,mat2,mat3])          
        
    
     

            
# Main evolution function (1D)     
def evolution(r,s,h,difWW,difDW,difDD,c,T,L,M,N,theta,conversion_timing) :  
    
    # Parameters initialization
    position = np.array([])   # list containing the first position where the proportion of wild alleles is higher than the threshold value.
    speed_fct_of_time = np.array([])      # speed computed... 
    time = np.array([])       # ...for each value of time in this vector.
    threshold = 0.5 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    pw_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    coex_density = -1       
      
    # Initial conditions for the system          
    W = np.ones(N+1); W[0:N//2] = 0    # Wild individuals at t=0  
    H = np.zeros(N+1)                  # Heterozygous individuals at t=0
    D = np.zeros(N+1); D[0:N//2] = 1   # Drive individuals at t=0

    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
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
                                
        
        # Compute the speed on WT allele density when coexistence density is close to zero
        # to avoid numerical errors (because then, proportions are very approximative), or when
        # there is no coexistence (coex_density = -1)         
        if coex_density < 0.1: 
            if conversion_timing == "zygote" : 
                WT = W+0.5*H
            if conversion_timing == "germline" : 
                WT = W+0.5*(1-c)*H
            p = WT/(W+H+D)
        # Else (coexistence state with final density > 0.1), compute the speed on WT allele density
        else : 
            if conversion_timing == "zygote" : 
                WT = (W+0.5*H)/(W+H+D) 
            if conversion_timing == "germline" : 
                WT = (W+0.5*(1-c)*H)/(W+H+D) 
            
        # Check if the WT wave is strictly increasing
        epsilon = -0.01 ; index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        if len(index_neg) > 4 : print("WT prop decreasing :", index_neg) 
               
        # Check if there is a coexistence state (if so adapt threshold and erase all speed values)
        # palliate centered
        center_index = np.arange(N//2-N//10,N//2+N//10); cond1 = np.all((p[center_index]>0.01) & (p[center_index]<0.99))
        # palliate left (if the wave move faster to the left)
        left_index = np.arange(N//2-N//5,N//2); cond2 = np.all((p[left_index]>0.01) & (p[left_index]<0.99))
        # palliate right (if the wave move faster to the left)
        right_index = np.arange(N//2,N//2+N//5); cond3 = np.all((p[right_index]>0.01) & (p[right_index]<0.99))
        if (cond1 or cond2 or cond3) and pw_star < 0 :    # if there exists a palliate in ]0,1[ around N//2 (and pw_star as not been changed yet)
            pw_star = np.mean(p[center_index])           
            coex_density = (W+H+D)[N//2]
            threshold = (pw_star+1)/2
            position = np.array([])   
            speed_fct_of_time = np.array([])   
            time = np.array([]) 
            print("coex ! p_w* =", pw_star)
            print("coex_density  =", coex_density) 
            print("threshold =", threshold)  
            
        # Check if there is a big jump in the speed (if so, erase all the precedent values)
        # This sometimes happens in coexistance, when the coex density is too close to 0.
        # Or in a preliminary phase where a decreasing section appears. (ex : germline r=10 s=0.92 h=0.51 c=0.85 cas a)
        if len(speed_fct_of_time) > 20 : 
            if abs(np.mean(speed_fct_of_time) - speed_fct_of_time[-1]) > 1: 
                print("Big jump in the speed for t =",t,"")
                position = np.array([])   
                speed_fct_of_time = np.array([])   
                time = np.array([]) 
                coex_density = -1
                pw_star = -1 
                threshold = 0.5
            
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
        
    # False coexistence detection
    if speed_fct_of_time[-1] < 0 : 
        print("False coex detection (speed<0)")
        coex_density = -1
        pw_star = -1 
        
    # return the speed, and 1 if there was a coexistence state (0 otherwise)
    return(speed_fct_of_time[-1], coex_density)




def heatmap_axes(y, rlog, precision):
    x_min = 0; x_max=1
    # x_delta is the size of a simulation pixel (size mesured with the s scale)  
    x_delta = (x_max-x_min)/precision    
    # x values for simulations (NB : x_min and x_max are not simulated, we simulate values centered on the simulation pixels)
    x_axis = np.arange(x_min+x_delta/2,x_max+x_delta/2,x_delta)      
    if y == "r" :                    
        if rlog : 
            y_min=0.01; y_max=10; y_delta = None
            y_axis = np.logspace(-2, 1, num=precision)
        else :
            y_min=0; y_max=12
            y_delta = (y_max-y_min)/precision  
            y_axis = np.arange(y_min+y_delta/2,y_max+y_delta/2,y_delta)               
    else :
        y_min = 0; y_max=1
        y_delta = (y_max-y_min)/precision    
        y_axis = np.arange(y_min+y_delta/2,y_max+y_delta/2,y_delta)       
    return(x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta)


            
def heatmap(precision, x, y, rlog, r, s, h, c):
    x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta = heatmap_axes(y, rlog, precision)
    speed_vect = np.zeros(precision)   
    coex_vect = np.zeros(precision) 
    y_index = num-1
    if y == "r" : r = y_axis[y_index]
    elif y == "s" : s = y_axis[y_index]
    elif y == "h" : h = y_axis[y_index]
    elif y == "c" : c = y_axis[y_index]
    
    for x_index in range(0, precision) :
        if x == "s" : s = x_axis[x_index]
        elif x == "h" : h = x_axis[x_index]
        elif x == "c" : c = x_axis[x_index]
        elif x == "r" : r = x_axis[x_index]
           
        print("\nh=", np.round(h,3)) 
        print("c=", np.round(c,3)) 
        print("r=", np.round(r,3)) 
        print("s=", np.round(s,3)) 
        speed_vect[x_index], coex_vect[x_index] = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,conversion_timing)
        print("speed=", np.round(speed_vect[x_index],3)) 
        print("coex=", coex_vect[x_index]) 

    return(speed_vect, coex_vect) 
  




############################### Parameters ######################################

## Biological
cas = "b_neg"                            # None "a" "b_pos" "b_neg" "c" "d"
x = "h"
y = "r"
rlog = True
r = 10                               # intrinsic growth rate
s = 0.9                              # fitness disadvantage for drive
h = 0.9                              # dominance coefficient
c = 0.85                             # conversion rate
conversion_timing = "germline"       # "zygote" or "germline"

# Particular growth/death terms
a = 0.2                              # coefficient for allee effect (growth or death)
growth_dynamic = "logistical"        # constant or logistical or allee_effect
death_dynamic = "constant"           # constant or logistical or allee_effect
# Different cases studied (if not concerned, use "cas = None" below)
if cas == "a" : growth_dynamic = "logistical"; death_dynamic = "constant"
if cas == "b_pos": growth_dynamic = "allee_effect"; death_dynamic = "constant"; a = 0.2
if cas == "b_neg": growth_dynamic = "allee_effect"; death_dynamic = "constant"; a = -0.2
if cas == "c" : growth_dynamic = "constant"; death_dynamic = "logistical"
if cas == "d" : growth_dynamic = "constant"; death_dynamic = "allee_effect"; a = 0.2

# Diffusion
difW = 1; difH = 1; difD = 1       # diffusion coefficient for resp. WW, WD or DD individuals

## Numerical
CI = "center"                      # Initial conditions : "center" for having the border in the center, "left" for having the border on the left
T = 1000                            # final time
L = 4000                           # length of the spatial domain
M = T*6                            # number of time steps
N = L                              # number of spatial steps
theta = 0.5                        # discretization in space : theta = 0.5 for Crank Nicholson, theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  
precision = 50       # heatmap precision
 
  



############################### Results #########################################

speed_vect, coex_vect = heatmap(precision, x, y, rlog, r, s, h, c)
np.savetxt(f'speed_{num}.txt', speed_vect) 
np.savetxt(f'coex_{num}.txt', coex_vect)  

  
        



