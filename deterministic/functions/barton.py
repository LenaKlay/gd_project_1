#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 10:06:31 2023

@author: lena
"""

import numpy as np
import matplotlib.pyplot as plt

T = 100                            # final time
L = 2000                           # length of the spatial domain
M = T*6                           # number of time steps
N = L                            # number of spatial steps

m = 0.1                   # migration rate
eps = L/N                 # spatial step size  
s = 0.2                   # fitness disadvantage
Mig = m/(eps**2)
rho = 1                   # density (cst in space ?)

p = np.ones(N)*0.8        # initialization
p_new = np.ones(N)*(-1)     
H = np.ones(M)*(-1)       # potential function
H[0] = np.sum( 2*eps*m*(rho**2) * (m/(2*eps**2)*(p[1:]-p[:-1])**2 + s*p[1:]**2*(1-p[1:])**2) )


for i in range(1,M) :
    t = np.round(np.linspace(0,T,M)[i],3)
    
    p_new[1:-1] = (1-2*Mig)*p[1:-1] + Mig*p[2:] + Mig*p[:-2] - s*p[1:-1]*(1-p[1:-1])*(1-2*p[1:-1])
    p_new[0] = (1-Mig)*p[0] + Mig*p[1] - s*p[0]*(1-p[0])*(1-2*p[0])
    p_new[-1] = (1-Mig)*p[-1] + Mig*p[-2] - s*p[-1]*(1-p[-1])*(1-2*p[-1])       
    p = p_new
        
    H[i] = np.sum( 2*eps*m*(rho**2) * (m/(2*eps**2)*(p[1:]-p[:-1])**2 + s*p[1:]**2*(1-p[1:])**2) )


plt.plot(np.linspace(0,T,M), H)
plt.show()
        
        
        
        
        

