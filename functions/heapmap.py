#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:32:47 2021

@author: lena
"""


############################## Libraries ########################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 

########################## Heatmap functions ####################################


def heatmap(heatmap_parameters, evolution, evolution_parameters, tanaka):
    # Parameters details
    precision, smin, smax, rmin, rmax, what = heatmap_parameters 
    r, s, h, difW, difH, difD, c, T, L, M, N, theta, mod, homing = evolution_parameters
    
    s_range = np.linspace(smin,smax,precision)             # range of values for s
    r_range = np.linspace(rmin,rmax,precision)             # range of values for r
    heatmap_values = np.zeros((precision,precision))       # matrix containing all the wave speeds for the different values of s and r 
    zero_line = np.matrix([[],[]])                         # separate the graph in two areas : positiv speed (right side) and negativ speed (left side)
    
    for r_index in range(0, precision) : 
        print("\n------- NEW r=", r_range[r_index])
        zero_done = False
        for s_index in range(0, precision) :
            r = r_range[r_index]
            s = s_range[s_index]
            
            print(f"\n s={np.round(s,2)},r={np.round(r,4)},c={np.round(c,2)},h={np.round(h,2)}")
            
            if what == "classic" :
                heatmap_values[r_index,s_index] = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)[0]
                if s_index != 0 and heatmap_values[r_index,s_index-1]*heatmap_values[r_index,s_index]<=0 and heatmap_values[r_index,s_index-1] != 0 and zero_done == False :
                    zero_line = np.append(zero_line,np.matrix([[s_index-0.5],[r_index]]), axis=1)
                    zero_done = True 
            if what == "speed_cubic" : 
                heatmap_values[r_index,s_index] = np.abs(evolution(r,s,1,1,1,1,1,T,L,M,N,0.5,mod)[0]-tanaka(s,T,L,M,N,0.5,mod,"cubic")[1])
            if what == "speed_fraction" : 
                heatmap_values[r_index,s_index] = np.abs(evolution(r,s,1,1,1,1,1,T,L,M,N,0.5,mod)[0]-tanaka(s,T,L,M,N,0.5,mod,"fraction")[1])
            if what == "r_one_minus_n_cubic" : 
                speed_girardin, W, H, D = evolution(r,s,1,1,1,1,1,T,L,M,N,0.5,mod)
                n = D+H+W
                heatmap_values[r_index,s_index] = np.max(abs(r*(1-n)))
            if what == "r_one_minus_n_fraction" : 
                speed_girardin, W, H, D = evolution(r,s,1,1,1,1,1,T,L,M,N,0.5,mod)
                p = D/(D+H+W); n = D+H+W
                heatmap_values[r_index,s_index] = np.max(abs(r*(1-n) - s*p*(2-p)/(1-s+s*(1-p)**2)))
            np.savetxt(f'heatmap_{homing}_c_{c}_h_{h}.txt', heatmap_values)  
            print(f"heatmap value ={heatmap_values[r_index,s_index]} \n")
                   
    np.savetxt(f'heatmap_{homing}_c_{c}_h_{h}_inside.txt', heatmap_values)               # heatmap_values = np.loadtxt('heatmap.txt') to get it back 
    if what == "classic" : np.savetxt(f'zero_line_{homing}_c_{c}_h_{h}_inside.txt', zero_line) 
    return(s_range, r_range, heatmap_values, zero_line) 
    
    
def print_heatmap(heatmap_parameters, print_heatmap_parameters) : 
    # Parameters details
    precision, smin, smax, rmin, rmax, what = heatmap_parameters 
    homing, c, h, s_range, r_range, heatmap_values, zero_line, style = print_heatmap_parameters 
        
    fig, ax = plt.subplots()     
    # Color choice
    if style == "simple" or style == "eradication" or style == "collapse" :
        colors1 = plt.cm.viridis(np.linspace(0, 0.8, 128))
        colors2 = plt.cm.plasma(np.flip(np.linspace(0., 1, 128)))
        colors = np.vstack((colors1, colors2))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        im = ax.imshow(heatmap_values,cmap=mymap, vmin=-4, vmax=4)
    else :
        im = ax.imshow(heatmap_values,cmap='Blues',interpolation='bicubic')

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
        
    #ax.set_title(f"Speed of the wave (c={c} and h={h}, {homing})",fontsize = 13)
    #fig.suptitle(f"Heatmap : {heatmap_type}", fontsize=14)
    ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
    ax.set_ylabel("r (growth rate)", fontsize=12)
    fig.tight_layout()
    plt.gca().invert_yaxis()
    plt.xlim([0,precision-1]);plt.ylim([0,precision-1])
    #plt.title("\n", f"Heatmap : {heatmap_type}")
     
    fig.savefig(f"heatmap_{homing}_c_{c}_h_{h}.png")
    #ax.xaxis.label.set_size(14)
    #ax.yaxis.label.set_size(14)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)
    plt.show()
    

        
# Pour retracer les heatmaps :
"""c=0.5;h=0.5;homing='germ line'; style = 'simple'
heatmap_type = 'classic'
heatmap_values = np.loadtxt(f'heatmap_{homing}_c_{c}_h_{c}.txt')  #_{homing}   'germ line'
heatmap_values = heatmap_values*(-1)
zero_line = np.loadtxt(f'zero_line_{homing}_c_{c}_h_{c}.txt')
precision = 20 
smin = 0.3; smax = 0.9; rmin = 0.1 ; rmax = 12    
s_range = np.linspace(smin,smax,precision)             
r_range = np.linspace(rmin,rmax,precision)
print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, style)"""

    