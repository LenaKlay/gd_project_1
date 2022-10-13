#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 18 14:22:58 2021

@author: lena
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 

# A checker avant de lancer :
# heatmap_type
# c et h
# homing
# rlog
    

################################ Paramètres ###################################
precision = 50
heatmap_type = "classic"       # "classic"  "speed_cubic" "speed_fraction" "r_one_minus_n_cubic"  "r_one_minus_n_fraction"                                           

homing = "zygote"
c = 1/4
h = 1/10
rlog = True

smin = 0.1; smax=0.9
s_range = np.linspace(smin,smax,precision)    

if rlog: 
    rmin=0.01; rmax=10
    r_range = np.logspace(-2, 1, num=precision)
else: 
    if heatmap_type == "classic" : 
        rmin=0; rmax=12; r_range = np.linspace(rmin,rmax,precision) 
    else : 
        rmin=50;rmax=60
        
plt.rcParams.update({'font.family':'serif'})

###############################################################################


def print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line) : 
        
    fig, ax = plt.subplots()     
    # Color choice
    if heatmap_type == "classic" :
        colors1 = plt.cm.Blues(np.flip(np.linspace(0.3, 1, 128))) 
        colors2 = plt.cm.hot(np.flip(np.linspace(0, 0.75, 128)))
        colors = np.vstack((colors1, colors2))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        im = ax.imshow(heatmap_values,cmap=mymap, vmin=-4, vmax=4)
    else :
        if heatmap_type == "r_one_minus_n_cubic" or heatmap_type == "r_one_minus_n_fraction" : 
            colors = plt.cm.plasma(np.flip(np.linspace(0, 1, 128)))
            mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
            im = ax.imshow(heatmap_values,cmap=mymap,interpolation='bicubic', vmin=0, vmax = 9)
        if heatmap_type == "speed_cubic" or heatmap_type == "speed_fraction" : 
            colors = plt.cm.viridis(np.flip(np.linspace(0, 1, 128)))
            mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
            im = ax.imshow(heatmap_values,cmap=mymap,interpolation='bicubic', vmin=0, vmax = 0.6)
    
    
    
    # Size of a simulation pixel, where the distance 'center of the first pixel' to 'center of the last pixel' is 1. 
    delta_square = 1/(precision-1)
    
    # Ticks positions : we want the ticks to start from the bordure, not the center
    ax.set_xticks(np.linspace(0-delta_square/2,1+delta_square/2,len(np.arange(smin,smax+0.1,0.1)))*(precision-1))                  
    if rlog : ax.set_yticks(np.linspace(0,1,4)*(precision-1))    
    else : ax.set_yticks(np.linspace(0-delta_square/2,1+delta_square/2,len(np.arange(int(rmin),rmax+1,1)))*(precision-1))                
        
    # Ticks labels
    ax.set_xticklabels(np.around(np.arange(smin,smax+0.1,0.1),2))                                         
    if rlog : ax.set_yticklabels(np.around(np.logspace(-2, 1, num=4),2))  
    else : ax.set_yticklabels(np.around(np.arange(int(rmin),rmax+1,1),2))   
    
  
    # Colorbar
    ax.figure.colorbar(im, ax=ax)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    # Problem if s = 1
    if np.isin(1, s_range) :                              # to avoid dividing by 0
        s_range = np.resize(s_range,len(s_range)-1); r_range = np.resize(r_range,len(r_range)-1)
        
    abscisse = (s_range-smin)*(precision-1)/(smax-smin)   
    
    # Eradication line
    if heatmap_type == "classic" :
        #if np.shape(zero_line)[1] != 0 : 
        #    ax.plot(np.array(zero_line[0,:]).ravel(),np.array(zero_line[1,:]).ravel(),color="red",label="zero speed", linewidth = 1.5, linestyle='-.')
        eradication_drive = (s_range/(1-s_range)-rmin)*((precision-1)/(rmax-rmin))
        ax.plot(abscisse,eradication_drive, color='orangered',label="eradication drive", linewidth = 2)  
        #ax.legend()
  
    #ax.set_title(f"(c={c} and h={h}, {homing})",fontsize = 13)
    #fig.suptitle(f"Travelling wave speed", fontsize=14)
    #ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
    #ax.set_ylabel("r (growth rate)", fontsize=12)
    fig.tight_layout()
    plt.gca().invert_yaxis()
    plt.xlim([0,precision-1]);plt.ylim([0,precision-1]) 
    fig.savefig(f'heatmap_c_{c}_h_{h}.png')    
    fig.savefig(f"heatmap_c_{c}_h_{h}.pdf", format='pdf')
    fig.savefig(f"heatmap_c_{c}_h_{h}.svg", format='svg')
    #ax.xaxis.label.set_size(14)
    #ax.yaxis.label.set_size(14)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)
    plt.show()
    
    
heatmap_values = np.zeros((precision,precision))
zero_line = np.matrix([[],[]]) 
for num in range(1,precision+1):
    heatmap_values[num-1,:] = np.loadtxt(f'heatmap_{num}.txt')
    if heatmap_type == "classic" : 
        zero_line = np.append(zero_line,np.matrix([[int(np.loadtxt(f'zero_line_{num}.txt'))],[num-1]]), axis=1)
        
np.savetxt(f'heatmap.txt', heatmap_values)
if heatmap_type == "classic" : np.savetxt(f'zero_line.txt', heatmap_values)

# heatmap_values = np.loadtxt(f'zero_line.txt')

print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line)





