#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 16:44:51 2022

@author: lena
"""

import numpy as np
import matplotlib.pyplot as plt

h = np.linspace(0,1,100)
fig, ax = plt.subplots()
for c in np.round(np.linspace(0,1,6),2) :
    ax.plot(h, c/(1-h*(1-c)), label=f's1 c={c}')
    #ax.plot(h, c/(2*c + h*(1-c)), label='zyg')
for c in np.round(np.linspace(0,1,6),2) :
    ax.plot(h, c/(2*c*h + h*(1-c)), label=f's2 c={c}')
ax.set(xlabel='s', ylabel="h", xlim = [0,1], ylim = [0,1]) 
ax.legend()
plt.show()



s = np.linspace(0.1,1,100)
fig, ax = plt.subplots()
col = ['blue','green','orange', 'deeppink']
i=0
for c in np.round(np.linspace(0.2,0.8,4),2) :
    ax.plot(s, (s-c)/(s*(1-c)), label=f'c={c}', color = col[i])
    ax.plot(s, c/(s*(1+c)), color = col[i])
    i=i+1
ax.plot(s, np.ones(100)*0.5, color = "black", linewidth = 1, linestyle='-.')
ax.plot(s, np.ones(100)*0.8, color = "black", linewidth = 1, linestyle='-.')
ax.plot(s, np.ones(100)*0.2, color = "black", linewidth = 1, linestyle='-.')
ax.set(xlabel='s', ylabel="h", xlim = [0,1.1], ylim = [-0.01,1]) 
ax.legend()
plt.show()



s = np.linspace(0.1,1,100)
fig, ax = plt.subplots()
col = ['blue','green','orange', 'deeppink']
i=0
for h in np.round(np.linspace(0.25,1,4),2) :
    ax.plot(s, s*(1-h)/(1-s*h), label=f'h={h}', color = col[i])
    ax.plot(s, s*h/(1-s*h), color = col[i])
    i=i+1
ax.plot(s, np.ones(100)*0.5, color = "black", linewidth = 1, linestyle='-.')
ax.plot(s, np.ones(100)*0.8, color = "black", linewidth = 1, linestyle='-.')
ax.plot(s, np.ones(100)*0.2, color = "black", linewidth = 1, linestyle='-.')
ax.set(xlabel='s', ylabel="c", xlim = [0,1.1], ylim = [-0.01,1]) 
ax.legend()
plt.show()

c = np.linspace(0,1,100)
fig, ax = plt.subplots()
col = ['blue','royalblue','lightskyblue','green','yellowgreen', 'gold', 'orange', 'deeppink','red']
i=0
for h in np.round(np.linspace(0.1,1,9),2) :
    ax.plot(c, c/(1-h*(1-c)), label=f'h={h}', color = col[i])
    i=i+1
ax.vlines(0.5, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.vlines(0.2, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.vlines(0.8, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.set(xlabel='c', ylabel="s", xlim = [0,1.1], ylim = [-0.01,1], title ="s1") 
ax.legend()
plt.show()

c = np.linspace(0,1,100)
fig, ax = plt.subplots()
col = ['blue','royalblue','lightskyblue','green','yellowgreen', 'gold', 'orange', 'deeppink','red']
i=0
for h in np.round(np.linspace(0.1,1,9),2) :
    ax.plot(c, c/(h*(1+c)), color = col[i])
    i=i+1
ax.vlines(0.5, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.vlines(0.2, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.vlines(0.8, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.set(xlabel='c', ylabel="s", xlim = [0,1.1], ylim = [-0.01,1], title ="s2") 
ax.legend()
plt.show()


s = 0.604
p_star = (2*s-1)/s
A = 1-s+s*(1-p_star)**2
fig, ax = plt.subplots()
col = ['blue','green','orange', 'deeppink']
i=0
for r in np.round(np.linspace(0,2,4),2) :
    n_star = 1 - (1-A)/(A*r)
    ax.plot(r, c/(h*(1+c)), color = col[i])
    i=i+1
ax.vlines(0.5, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.vlines(0.2, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.vlines(0.8, 0, 1, color = "black", linewidth = 1, linestyle='-.')
ax.set(xlabel='c', ylabel="s", xlim = [0,1.1], ylim = [-0.01,1]) 
ax.legend()
plt.show()


s = 0.604
p_star = (2*s-1)/s
A = 1-s+s*(1-p_star)**2
for r in [0.12,  0.36,  0.6 ,  0.84,  1.08,  1.32,  1.56,  1.8 ,  2.04,
        2.28,  2.52,  2.76,  3.  ,  3.24,  3.48,  3.72,  3.96,  4.2 ,
        4.44,  4.68,  4.92,  5.16,  5.4 ,  5.64,  5.88,  6.12,  6.36,
        6.6 ,  6.84,  7.08,  7.32,  7.56,  7.8 ,  8.04,  8.28,  8.52,
        8.76,  9.  ,  9.24,  9.48,  9.72,  9.96, 10.2 , 10.44, 10.68,
       10.92, 11.16, 11.4 , 11.64, 11.88] :
    n_star = 1 - (1-A)/(A*r)
    p_star = (s_range*(1-(1-c)*(1-h)) - c*(1-s_range))/(s_range*(1-2*(1-c)*(1-h)))  
                mean_fitness = (1-s_range)*p_star**2+2*(c*(1-s_range)+(1-c)*(1-s_range*h))*p_star*(1-p_star)+(1-p_star)**2 
    print(r)
    print('n',n_star)
    print('p',p_star)
