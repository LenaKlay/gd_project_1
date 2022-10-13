#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:32:47 2021

@author: lena
"""


# Libraries
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
plt.rcParams.update({'font.family':'serif'})

# External functions 
from evolution import evolution
from graph import save_fig_or_data




#### Give the heatmap axes (r in log scale or not) ####
    
def heatmap_axes(rlog, precision):
    smin = 0.1; smax=0.9
    # delta_s is the size of a simulation pixel (size mesured with the s scale)  
    delta_s = (smax-smin)/precision    
    # s values for simulations (NB : smin and smax are not simulated, we simulate values centered on the simulation pixels)
    s_range = np.arange(smin+delta_s/2,smax+delta_s/2,delta_s)      
    # r in log scale                   
    if rlog : 
        rmin=0.01; rmax=10; delta_r = None
        r_range = np.logspace(-2, 1, num=precision)
    # r in normal scale
    else : 
        rmin=0; rmax=12
        # delta_r is the size of a simulation pixel (size mesured with the r scale)  
        delta_r = (rmax-rmin)/precision    
        # r values for simulations (NB : rmin and rmax are not simulated, we simulate values centered on the simulation pixels)
        r_range = np.arange(rmin+delta_r/2,rmax+delta_r/2,delta_r)       
    return(s_range, smin, smax, delta_s, r_range, rmin, rmax, delta_r)



#### Create an heatmap (speed of the traveling wave for couples (s,r) ####

def heatmap(bio_para, model_para, num_para, graph_para, rlog, precision):
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    s_range, smin, smax, delta_s, r_range, rmin, rmax, delta_r = heatmap_axes(rlog, precision)
    
    # Create a directory and save parameters.txt with r_range and s_range
    bio_para[0] = r_range; bio_para[1] = s_range
    save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", [], [], None, bio_para, num_para)
     
    # Print principal parameters
    print("conversion_timing =", bio_para[8]);  print("\nr =", r_range); print("s =", s_range); print("h =", h); print("c =", c)
            
    # Initialization
    # matrix containing all the wave speeds for the different values of s and r 
    heatmap_values = np.zeros((precision,precision))      

    # Loop to calculate the speed, for each couple (s,r)
    for r_index in range(0, precision) : 
        print("\n------- NEW r=", r_range[r_index])
        for s_index in range(0, precision) :
            print("\ns=", np.round(s_range[s_index],3))             
            # Update values in bio_para
            bio_para[0] = r_range[r_index]
            bio_para[1] = s_range[s_index]                    
            # Speed value given by evolution.py
            heatmap_values[r_index,s_index] = evolution(bio_para, model_para, num_para, graph_para)[4][-1] 
            print("speed :", heatmap_values[r_index,s_index])
                     
        # for each r, save the corresponding line
        save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", [], heatmap_values[r_index,:], f"r_line_{r_index}", bio_para, num_para)   
        
    # Save final data
    save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", [], heatmap_values, f"{conversion_timing}_c_{c}_h_{h}", bio_para, num_para)
    return(heatmap_values) 
 
    


       
#### Load an already existing heatmap ####
    
def load_heatmap(conversion_timing, c, h, rlog, precision, migale) : 
    # migale = True means we download each line of the heatmap (from the migale file)
    if migale :
        heatmap_values = np.zeros((precision,precision))
        for r_index in range(0,precision):
            heatmap_values[r_index,:] = np.loadtxt(f'../migale/heatmaps/heatmap_{r_index+1}.txt')
        if not os.path.exists(f'../outputs/heatmaps/{conversion_timing}_c_{c}_h_{h}'): 
            os.mkdir(f'../outputs/heatmaps/{conversion_timing}_c_{c}_h_{h}')
        np.savetxt(f'../outputs/heatmaps/{conversion_timing}_c_{c}_h_{h}/{conversion_timing}_c_{c}_h_{h}.txt', heatmap_values) 
    # migale = False means we download the full data from a .txt, either in the file : where = 'figures' or in the file : where = 'outputs'
    else :       
        heatmap_values = np.loadtxt(f'../figures/heatmaps/{conversion_timing}_c_{c}_h_{h}/{conversion_timing}_c_{c}_h_{h}.txt')  
    return(heatmap_values)
  
  
    
#### How do we store the values in heatmap_values ? ####
# The heatmap_values[r_index,s_index] correspond to the values s : s_range[s_index] and r : r_range[r_index]
#indice_r = np.where((5.3 < r_range) & (r_range < 5.5))[0] ; print("\nindice r :", indice_r)
#indice_s = np.where((0.58 < s_range) & (s_range < 0.6))[0] ; print("\nindice s :", indice_s)
#print(heatmap_values[indice_r,indice_s])


#### Print an heatmap from heatmap_values ####
    
def print_heatmap(heatmap_values, bio_para, num_para, rlog, precision) :    
    r, s, h, a, difWW, difDW, difDD, c, conversion_timing = bio_para
    s_range, smin, smax, delta_s, r_range, rmin, rmax, delta_r = heatmap_axes(rlog, precision)
    
    # Figure
    fig, ax = plt.subplots()     
    
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
     
    # Colorbar creation
    colors1 = plt.cm.Blues(np.flip(np.linspace(0.3, 1, 128))); colors2 = plt.cm.hot(np.flip(np.linspace(0, 0.75, 128)))  # Create the color scale   
    colors = np.vstack((colors1, colors2)) # Merge the two
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)  
    # Plot heatmap values
    im = ax.imshow(heatmap_values,cmap=mymap, vmin=-4, vmax=4, aspect='auto')  
    # Add the colorbar
    ax.figure.colorbar(im, ax=ax)
    
    # Add another x axis on the top of the heatmap to indicate s_1 and s_2
    s_1 = c/(1-h*(1-c))
    if conversion_timing == "zygote" : 
        s_2 = c/(2*c + h*(1-c))
    if conversion_timing == "germline" : 
        s_2 = c/(2*c*h + h*(1-c))  
    axtop = ax.twiny(); ticks = []; ticklabels = []
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
        
    # Plot lines    
    abscisse = np.arange(precision)    
    # Tool to draw precise line
    if rlog : r_range_precise = np.logspace(-2, 1, num=precision*100)
    else : r_range_precise = np.arange(rmin+delta_r/2, rmax+delta_r/2, delta_r/100)
        
    # Pure drive persistance line
    if rlog :
        eradication_drive = np.zeros(precision)
        for i in range(precision):
            s = s_range[i]
            eradication_drive[i] = np.where(s/(1-s) < r_range_precise)[0][0]/100         
    else : 
        eradication_drive = (s_range/(1-s_range)-rmin)*((precision-1)/(rmax-rmin)) 
    ax.plot(abscisse,eradication_drive-0.5, color='#73c946ff', label="eradication drive", linewidth = 4)           
                 
    # Composite persistance line
    if (conversion_timing == "zygote" and (1-h)*(1-c) > 0.5) or (conversion_timing == "germline" and h < 0.5) :                        
        # neglect the values outside s_1 s_2 (which do not interest us)
        s1_s2_len = 200
        s1_s2_range = np.linspace(s_1,s_2,s1_s2_len) 
             
        if conversion_timing == "zygote" :    
            p_star = (s1_s2_range*(1-(1-c)*(1-h)) - c*(1-s1_s2_range))/(s1_s2_range*(1-2*(1-c)*(1-h)))  
            mean_fitness = (1-s1_s2_range)*p_star**2+2*(c*(1-s1_s2_range)+(1-c)*(1-s1_s2_range*h))*p_star*(1-p_star)+(1-p_star)**2 
        if conversion_timing == "germline" :           
            p_star = ((1-s1_s2_range*h)*(1+c)-1)/(s1_s2_range*(1-2*h))
            mean_fitness = (1-s1_s2_range)*p_star**2+2*(1-s1_s2_range*h)*p_star*(1-p_star)+(1-p_star)**2   
      
        eradication_pop = np.zeros(s1_s2_len)
        for i in range(s1_s2_len):
            eradication_pop[i] = np.where((1-mean_fitness[i])/mean_fitness[i] < r_range_precise)[0][0]/100 
        eradication_pop = eradication_pop[0:np.where(eradication_pop==0)[0][0]+1]
        abscisse_pop = ((s1_s2_range-s_range[0])*((precision-1)/(s_range[-1]-s_range[0])))[0:np.where(eradication_pop==0)[0][0] + 1]
        ax.plot(abscisse_pop, eradication_pop-0.5, color='#40720cff', linewidth = 4) 
               
    # Axis and label sizes
    plt.gca().invert_yaxis()                 # inverse r axis (increasing values)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)

    # Always save figure
    save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", fig, [], f"{conversion_timing}_c_{c}_h_{h}", bio_para, num_para)

    plt.show()

