#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""

num = 149

############################### Load Stuff #######################################

# Load libraries
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg  as la

############################### Parameters ######################################


# Parameters
T = 3200; L = 4800; M = T*120; N = L*20; theta = 0.5 

# Set the x-axis (values of s)
s_min = 0.34 ; s_max = 0.4
s_values = np.linspace(s_min,s_max,200)
s = np.round(s_values[num],3)


# Steps
dt = T/M    # time
dx = L/N    # spatial
# Spatial domain (1D)
X = np.linspace(0,N,N+1)*dx   
    
# Initialization       
P = np.zeros(N+1); P[0:N//5] = 1              # Proportion of drive individuals at t=0
    
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
        if num%2:
            f = (s*P*(1-P)*(P-(2*s-1)/s))/(1-s+s*(1-P)**2)  
        else:
            f = P*(1-P)*(1-2*s)
      
        P = la.spsolve(B_, B.dot(P) + dt*f)
            
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
        
print(speed_fct_of_time[-1])
np.savetxt(f"{num}.txt", np.ones(1)*speed_fct_of_time[-1])  






                




