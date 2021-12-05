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
# if heatmap_type == "classic" -> style (else stype = None)
# c et h
# homing
    

################################ Paramètres ###################################
precision = 30
heatmap_type = "classic"       # "classic"  "speed_cubic" "speed_fraction" "r_one_minus_n_cubic"  "r_one_minus_n_fraction"                                           
style = "simple"               # if heatmap_type != "classic : simple", "eradication" or "collapse"
                               # else : None

homing = "zygote"
c = 3/4
h = 3/4
smin = 0.3; smax=0.9

if heatmap_type == "classic" : 
    rmin=0;rmax=12
else : 
    rmin=50;rmax=60
    
s_range = np.linspace(smin,smax,precision)          
r_range = np.linspace(rmin,rmax,precision)          
plt.rcParams.update({'font.family':'serif'})

###############################################################################


def print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, style) : 
        
    fig, ax = plt.subplots()     
    # Color choice
    if style == "simple" or style == "eradication" or style == "collapse" :
        colors1 = plt.cm.viridis(np.linspace(0, 0.8, 128))
        colors2 = plt.cm.plasma(np.flip(np.linspace(0., 1, 128)))
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

    ax.set_xticks(np.linspace(0,1,len(np.arange(smin,smax,0.1)))*(precision-1))        # (np.linspace(0,1,int((smax-smin)*10+1))*(precision-1))
    ax.set_yticks(np.linspace(0,1,len(np.arange(int(rmin),rmax+1,1)))*(precision-1))                      # ax.set_yticks(np.linspace(0,1,len(np.arange(int(rmin),rmax+1,2))))   
    ax.set_xticklabels(np.around(np.arange(smin,smax,0.1),2))                         # (np.around(np.linspace(smin,smax,int((smax-smin)*10+1)),3))                    
    ax.set_yticklabels(np.around(np.arange(int(rmin),rmax+1,1),2))             # ax.set_yticklabels(np.arange(int(rmin),rmax+1,2))   
    # Colorbar
    ax.figure.colorbar(im, ax=ax)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    # Problem if s = 1
    if np.isin(1, s_range) :                              # to avoid dividing by 0
        s_range = np.resize(s_range,len(s_range)-1); r_range = np.resize(r_range,len(r_range)-1)
        
    abscisse = (s_range-smin)*(precision-1)/(smax-smin)   
    
    # Eradication line
    if style == "simple" or style == "eradication" or style == "collapse" :
        if np.shape(zero_line)[1] != 0 : 
            ax.plot(np.array(zero_line[0,:]).ravel(),np.array(zero_line[1,:]).ravel(),color="red",label="zero speed", linewidth = 1.5, linestyle='-.')
        eradication_drive = (s_range/(1-s_range)-rmin)*((precision-1)/(rmax-rmin))
        ax.plot(abscisse,eradication_drive, color='orangered',label="eradication drive", linewidth = 2)  
        ax.legend()
    # Eradication zone
    if style == "eradication":
        abscisse_for_zero_line = np.unique(np.array(zero_line[0,:]).ravel(),return_index=True)
        unique_zero_line = np.array(zero_line[1, abscisse_for_zero_line[1]]).ravel()
        abscisse_for_zero_line_into_s = abscisse_for_zero_line[0]*(smax-smin)/(precision-1)+smin
        eradication_drive_for_zero_line = (abscisse_for_zero_line_into_s/(1-abscisse_for_zero_line_into_s)-rmin)*((precision-1)/(rmax-rmin))
        plt.fill_between(abscisse, eradication_drive, where=np.round(abscisse,2)<=zero_line[0,0], color='orangered')
        plt.fill_between(abscisse_for_zero_line[0], unique_zero_line, eradication_drive_for_zero_line, where= eradication_drive_for_zero_line>=unique_zero_line, color='orangered')
    
    # Collapse line
    if style == "collapse":
        if homing == "zygote":
            collapsing_drive = ((1/((c+1)*(1-s_range)+(1-c)*(1-s_range*h))-1)-rmin)*((precision-1)/(rmax-rmin))
        if homing == "germline":
            collapsing_drive = ((1/((1-s_range)+(c+1)*(1-s_range*h))-1)-rmin)*((precision-1)/(rmax-rmin))
        ax.plot(abscisse,collapsing_drive, color='deepskyblue',label="collapsing drive", linewidth = 2)
    # Collapsing zone 
        plt.fill_between(abscisse, collapsing_drive, color='deepskyblue')
        
    ax.set_title(f"(c={c} and h={h}, {homing})",fontsize = 13)
    fig.suptitle(f"Travelling wave speed", fontsize=14)
    ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
    ax.set_ylabel("r (growth rate)", fontsize=12)
    fig.tight_layout()
    plt.gca().invert_yaxis()
    plt.xlim([0,precision-1]);plt.ylim([0,precision-1])
    #plt.title("\n", f"Heatmap : {heatmap_type}")
     
    fig.savefig(f'{heatmap_type}_heatmap.png')
    #ax.xaxis.label.set_size(14)
    #ax.yaxis.label.set_size(14)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)
    plt.show()
    
    
heatmap_values = np.zeros((precision,precision))
zero_line = np.matrix([[],[]]) 
for num in range(1,precision+1) :
    heatmap_values[num-1,:] = np.loadtxt(f'{heatmap_type}_heatmap_{num}.txt')
    if heatmap_type == "classic" : 
        zero_line = np.append(zero_line,np.matrix([[int(np.loadtxt(f'{heatmap_type}_zero_line_{num}.txt'))],[num-1]]), axis=1)
        
np.savetxt(f'{heatmap_type}_heatmap.txt', heatmap_values)
if heatmap_type == "classic" : np.savetxt(f'{heatmap_type}_zero_line.txt', heatmap_values)

# heatmap_values = np.loadtxt(f'{heatmap_type}_zero_line.txt')

print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, style)









