#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 12:23:35 2022

@author: lena
"""


import numpy as np
import matplotlib.pyplot as plt

# Colors used
col_pink = ['indigo', 'purple', 'darkmagenta', 'm', 'mediumvioletred', 'crimson', 'deeppink', 'hotpink', 'lightpink', 'pink' ]    
col_blue = ['navy', 'blue','royalblue', 'cornflowerblue', 'lightskyblue']    


# Change font to serif
plt.rcParams.update({'font.family':'serif'})

title_size = 15
label_size = 17
legend_size = 12
line_size = 3

s_min = 0.34 ; s_max = 0.4
nb_points = 100
fct_of_s = np.zeros((4,nb_points)) 
fct_of_s[0,:] = np.linspace(s_min, s_max, nb_points)
fct_of_s[1,:] = 2*np.sqrt(1-2*fct_of_s[0,:])
    
for i in range(0,100):
    file_frac = open(f"{2*i}.txt", "r")
    fct_of_s[2,:][i] = file_frac.read()
    file_frac.close()
    file_kpp = open(f"{2*i+1}.txt", "r")
    fct_of_s[3,:][i] = file_kpp.read()
    file_kpp.close()

np.savetxt(f"0_speed_frac.txt", fct_of_s[2,:])  
np.savetxt(f"0_speed_kpp.txt", fct_of_s[3,:])  

fig, ax = plt.subplots()
ax.plot(fct_of_s[0,:], fct_of_s[1,:], label = 'kpp theo', color = col_pink[0]) 
ax.plot(fct_of_s[0,:], fct_of_s[2,:], label = 'kpp num', color = col_pink[4]) 
ax.plot(fct_of_s[0,:], fct_of_s[3,:], label = 'r inf', color = 'orange') 
ax.set(xlabel='s', ylabel='Speed') 
ax.xaxis.label.set_size(label_size); ax.yaxis.label.set_size(label_size)  
plt.legend() 
plt.grid()  
fig.savefig(f"speed.svg", format='svg')
plt.show() 
