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
c = 1
h = 1
rlog = True



# Range for r and s
smin = 0.1; smax=0.9
delta_s = (smax-smin)/precision    # delta_s is the size of a simulation pixel (size mesured with the s scale)  
s_range = np.arange(smin+delta_s/2,smax+delta_s/2,delta_s)       # s values for simulations (NB : smin and smax are not simulated, we simulate values centered on the simulation pixels)                  
if rlog : 
    rmin=0.01; rmax=10
    r_range = np.logspace(-2, 1, num=precision)
else : 
    rmin=0; rmax=12
    delta_r = (rmax-rmin)/precision    # delta_r is the size of a simulation pixel (size mesured with the r scale)  
    r_range = np.arange(rmin+delta_r/2,rmax+delta_r/2,delta_r)       # r values for simulations (NB : rmin and rmax are not simulated, we simulate values centered on the simulation pixels)

plt.rcParams.update({'font.family':'serif'})

###############################################################################

    
def print_heatmap(homing, c, h, s_range, r_range, heatmap_values) : 
        
    # Figure
    fig, ax = plt.subplots()     
    # Color choice
    colors1 = plt.cm.Blues(np.flip(np.linspace(0.3, 1, 128))) #plt.cm.ocean(np.linspace(0.28, 0.95, 128))
    # Positive scale
    colors2 = plt.cm.hot(np.flip(np.linspace(0, 0.75, 128)))
    # Merge the two
    colors = np.vstack((colors1, colors2))
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
    # Plot heatmap values
    im = ax.imshow(heatmap_values,cmap=mymap, vmin=-4, vmax=4, aspect='auto')

    # Size of a simulation pixel, where the distance 'center of the first pixel' to 'center of the last pixel' is 1. (useful when rlog=False) 
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
    # Problem if s = 1
    if np.isin(1, s_range) :                              # Avoid dividing by 0
        s_range = np.resize(s_range,len(s_range)-1); r_range = np.resize(r_range,len(r_range)-1)
        
    # Abscisse to plot lines    
    abscisse = np.arange(precision)       # (s_range-smin)*(precision-1)/(smax-smin)   
    
    # Lines always there
    if heatmap_type == "classic" :  
        
        # Tool to draw precise line
        if rlog : r_range_precise = np.logspace(-2, 1, num=precision*100)
        else : r_range_precise = np.arange(rmin+delta_r/2, rmax+delta_r/2, delta_r/100)    

        # Eradication drive line
        if rlog :
            eradication_drive = np.zeros(precision)
            for i in range(precision):
                s = s_range[i]
                eradication_drive[i] = np.where(s/(1-s) < r_range_precise)[0][0]/100         
        else : 
            eradication_drive = (s_range/(1-s_range)-rmin)*((precision-1)/(rmax-rmin)) 
        ax.plot(abscisse,eradication_drive-0.5, color='#73c946ff', label="eradication drive", linewidth = 4) 
           
                
        # Values s1 and s2
        delta_s = s_range[1]-s_range[0]   # = (smax-smin)/precision
        s_1 = c/(1-h*(1-c))
        if homing == "zygote" : 
            s_2 = c/(2*c + h*(1-c))
        if homing == "germline" : 
            s_2 = c/(2*c*h + h*(1-c))
         
        # Add another x axis on the top of the heatmap to indicate s_1 and s_2
        axtop = ax.twiny()
        ticks = []
        ticklabels = []
        if smin < s_1 and s_1 < smax : 
            ticks.append((s_1-smin)/(smax-smin)); ticklabels.append("s1")
            # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
            ax.vlines((s_1-s_range[0])*(1/delta_s),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))        
        if smin < s_2 and s_2 < smax : 
            ticks.append((s_2-smin)/(smax-smin)); ticklabels.append("s2")
            # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
            ax.vlines((s_2-s_range[0])*(1/delta_s),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))
        axtop.set_xticks(ticks) 
        axtop.set_xticklabels(ticklabels)
          
        # Population eradication line in case of coexistence state
        if (homing == "zygote" and (1-h)*(1-c) > 0.5) or (homing == "germline" and h < 0.5) :
                        
             # neglect the values outside s_1 s_2 (which do not interest us)
             s1_s2_len = 100
             s1_s2_range = np.linspace(s_1,s_2,s1_s2_len) 
             
             if homing == "zygote" :    
                p_star = (s1_s2_range*(1-(1-c)*(1-h)) - c*(1-s1_s2_range))/(s1_s2_range*(1-2*(1-c)*(1-h)))  
                mean_fitness = (1-s1_s2_range)*p_star**2+2*(c*(1-s1_s2_range)+(1-c)*(1-s1_s2_range*h))*p_star*(1-p_star)+(1-p_star)**2 
             if homing == "germline" : 
                p_star = ((1-s1_s2_range*h)*(1+c)-1)/(s1_s2_range*(1-2*h))
                mean_fitness = (1-s1_s2_range)*p_star**2+2*(1-s1_s2_range*h)*p_star*(1-p_star)+(1-p_star)**2   
            
            
     ################################ Verifier pour rlog = False; si ok copier le tout dans heatmap.py ###################################
       
             eradication_pop = np.zeros(s1_s2_len)
             for i in range(s1_s2_len):
                 eradication_pop[i] = np.where((1-mean_fitness[i])/mean_fitness[i] < r_range_precise)[0][0]/100 
             eradication_pop = eradication_pop[0:np.where(eradication_pop==0)[0][0] + 1]
             abscisse_pop = ((s1_s2_range-s_range[0])*((precision-1)/(s_range[-1]-s_range[0])))[0:np.where(eradication_pop==0)[0][0] + 1]
             ax.plot(abscisse_pop, eradication_pop-0.5, color='#40720cff', linewidth = 4) 
         
             
    plt.gca().invert_yaxis()                 # inverse r axis (increase)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)
    
    fig.savefig(f"{path}/{homing}_c_{c}_h_{h}.png", format='png')
    fig.savefig(f"{path}/{homing}_c_{c}_h_{h}.pdf", format='pdf') #; fig.savefig(f"../outputs/{directories}/{title}.png") 
    fig.savefig(f"{path}/{homing}_c_{c}_h_{h}.svg", format='svg')
  
    plt.show()
    
path = f'../figures/heatmaps/{homing}_c_{c}_h_{h}'
heatmap_values = np.loadtxt(f'{path}/heatmap_c_{c}_h_{h}.txt') 
print_heatmap(homing, c, h, s_range, r_range, heatmap_values)
