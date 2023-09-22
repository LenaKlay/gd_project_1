#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 18:53:37 2022

@author: lena
"""

### Librairies

import numpy as np
import matplotlib.pyplot as plt

### Graphs parameters 

plt.rcParams.update({'font.family':'serif'})
label_size = 12
legend_size = 10


### Parameters

#K = 10**8 # 1000000000       # Carrying capacity for a spatial interval of size 1
dx = 1                      # Spatial interval 
T = 1000                    # Time at which simulation stops (it can stop earlier if one of the type disappear from the environment)
m = 0.2                     # Migration probability
dt = np.round(m*dx**2/2,10) # Time interval

conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
r = 0.1                     # Intrasic growth rate
c = 0.9                     # Conversion rate
h = 0.4                     # Dominance
#s = 0.9                     # Disadvantage for drive

#col = ["darkturquoise", "deepskyblue", "dodgerblue", "royalblue", "blue", "navy"]        
col = ["plum", "orchid", "m", "darkviolet", "indigo", "navy"]          
K_range = [10**3, 10**4, 10**5, 10**6, 10**7, 10**8]
s_range = np.round(np.arange(0.1,1,0.05),3)
chasing_matrix = np.ones((len(K_range), len(s_range)))*(-1)

for i in range(len(K_range)):
    K = K_range[i]
    for j in range(len(s_range)):         
        s = s_range[j]
        chasing_matrix[i,j] = np.mean(np.loadtxt(f"res/res_K_{int(np.log10(K))}_s_{s}.txt"))


fig, ax = plt.subplots()
for i in (range(len(K_range)))[::-1]:
    ax.plot(s_range, chasing_matrix[i,:], label = f"$K=10^{int(np.log10(K_range[i]))}$", color = col[i])
    ax.scatter(s_range, chasing_matrix[i,:], alpha=0.5, color = col[i])
ax.set(xlabel='s (fitness disadvantage for drive)', ylabel='Proportion of chasing') 
#ax.set_title(f"Proportion of chasing among 100 runs of the stochastic simulation")
ax.xaxis.label.set_size(label_size)
ax.yaxis.label.set_size(label_size)
plt.rc('legend', fontsize=legend_size)  
plt.legend()   
fig.savefig(f"chasing_s.png", format='png')
plt.show()

fig, ax = plt.subplots()
for j in range(len(s_range)):
    ax.semilogx(K_range, chasing_matrix[:,j], label = f"s={s_range[j]}")
ax.set(xlabel='K (carrying capacity)', ylabel='Proportion of chasing') 
#ax.set_title(f"Proportion of chasing among 100 runs of the stochastic simulation")
ax.xaxis.label.set_size(label_size)
ax.yaxis.label.set_size(label_size)   
plt.rc('legend', fontsize=legend_size)  
plt.legend()
fig.savefig(f"chasing_k.png", format='png')
plt.show()


fig, ax = plt.subplots()
im = ax.imshow(chasing_matrix)  
ax.figure.colorbar(im, ax=ax)
plt.show()