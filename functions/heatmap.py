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
    x_min = 0; x_max=1
    # x_delta is the size of a simulation pixel (size mesured with the s scale)  
    x_delta = (x_max-x_min)/precision    
    # x values for simulations (NB : x_min and x_max are not simulated, we simulate values centered on the simulation pixels)
    x_axis = np.arange(x_min+x_delta/2,x_max+x_delta/2,x_delta)      
    # r in log scale                   
    if rlog == True : 
        y_min=0.01; y_max=10; y_delta = None
        y_axis = np.logspace(-2, 1, num=precision)
    # r in normal scale
    if rlog == False : 
        y_min=0; y_max=12
        # y_delta is the size of a simulation pixel (size mesured with the r scale)  
        y_delta = (y_max-y_min)/precision    
        # y values for simulations (NB : y_min and y_max are not simulated, we simulate values centered on the simulation pixels)
        y_axis = np.arange(y_min+y_delta/2,y_max+y_delta/2,y_delta)               
    # y_axis do not represent r (more probably s)
    if rlog == None : 
        y_min = 0; y_max=1
        # y_delta is the size of a simulation pixel (size mesured with the s scale)  
        y_delta = (y_max-y_min)/precision    
        # y values for simulations (NB : y_min and y_max are not simulated, we simulate values centered on the simulation pixels)
        y_axis = np.arange(y_min+y_delta/2,y_max+y_delta/2,y_delta)       
    return(x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta)



#### Create an heatmap (speed of the traveling wave for couples (s,r) ####

def heatmap(bio_para, model_para, num_para, graph_para, rlog, precision):
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta = heatmap_axes(rlog, precision)
    
    # Create a directory and save parameters.txt with y_axis and x_axis
    bio_para[0] = y_axis; bio_para[1] = x_axis
    save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", [], [], None, bio_para, num_para, model_para)
     
    # Print principal parameters
    print("conversion_timing =", bio_para[8]);  print("\nr =", y_axis); print("s =", x_axis); print("h =", h); print("c =", c)
            
    # Initialization
    # matrix containing all the wave speeds for the different values of s and r 
    heatmap_values = np.zeros((precision,precision))      

    # Loop to calculate the speed, for each couple (s,r)
    for r_index in range(0, precision) : 
        print("\n------- NEW r=", y_axis[r_index])
        for s_index in range(0, precision) :
            print("\ns=", np.round(x_axis[s_index],3))             
            # Update values in bio_para
            bio_para[0] = y_axis[r_index]
            bio_para[1] = x_axis[s_index]                    
            # Speed value given by evolution.py
            heatmap_values[r_index,s_index] = evolution(bio_para, model_para, num_para, graph_para)[4][-1] 
            print("speed :", heatmap_values[r_index,s_index])
                     
        # for each r, save the corresponding line
        save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", [], heatmap_values[r_index,:], f"r_line_{r_index}", bio_para, num_para, model_para)   
        
    # Save final data
    save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_h_{h}", [], heatmap_values, f"{conversion_timing}_c_{c}_h_{h}", bio_para, num_para, model_para)
    return(heatmap_values) 
 
    


       
#### Load an already existing heatmap ####
    
def load_heatmap(conversion_timing, c, r, h, s, rlog, precision, migale, cas) : 
    heatmap_values = np.zeros((precision,precision))
    coex_values = np.ones((precision,precision))*(-2)
    # migale = True means we download each line of the heatmap (from the migale file)
    if migale :
        for r_index in range(0,precision):
            if os.path.exists(f'../migale/heatmaps/heatmap_{r_index+1}.txt') :
                heatmap_values[r_index,:] = np.loadtxt(f'../migale/heatmaps/heatmap_{r_index+1}.txt')
            else : print(f'Manque heatmap_{r_index+1}.txt')
            if os.path.exists(f'../migale/heatmaps/coex_{r_index+1}.txt') :
                coex_values[r_index,:] = np.loadtxt(f'../migale/heatmaps/coex_{r_index+1}.txt')
            else : print(f'Manque coex_{r_index+1}.txt')
        if not os.path.exists(f'../outputs/heatmaps/{conversion_timing}_c_{c}_r_{r}'): 
            os.mkdir(f'../outputs/heatmaps/{conversion_timing}_c_{c}_r_{r}')
        np.savetxt(f'../outputs/heatmaps/{conversion_timing}_c_{c}_r_{r}/{conversion_timing}_c_{c}_r_{r}.txt', heatmap_values) 
        np.savetxt(f'../outputs/heatmaps/{conversion_timing}_c_{c}_r_{r}/{conversion_timing}_c_{c}_r_{r}_coex.txt', coex_values) 
    # migale = False means we download the full data from a .txt, either in the file : where = 'figures' or in the file : where = 'outputs'
    else :  
        if cas == None :
            heatmap_values = np.loadtxt(f'../figures/heatmaps/{conversion_timing}_c_{c}_h_{h}/{conversion_timing}_c_{c}_h_{h}.txt')  
            if os.path.exists(f'../figures/heatmaps/{conversion_timing}_c_{c}_h_{h}/{conversion_timing}_c_{c}_h_{h}_coex.txt') :
                coex_values = np.loadtxt(f'../figures/heatmaps/{conversion_timing}_c_{c}_h_{h}/{conversion_timing}_c_{c}_h_{h}_coex.txt')     
        else : 
            heatmap_values = np.loadtxt(f'../outputs/heatmaps/cas_{cas}/{conversion_timing}_c_{c}_r_{r}.txt')  
            if os.path.exists(f'../outputs/heatmaps/cas_{cas}/{conversion_timing}_c_{c}_r_{r}_coex.txt') :
                coex_values = np.loadtxt(f'../outputs/heatmaps/cas_{cas}/{conversion_timing}_c_{c}_r_{r}_coex.txt')     
    return(heatmap_values, coex_values)
  
    
#### How do we store the values in heatmap_values ? ####
# The heatmap_values[r_index,s_index] correspond to the values s : x_axis[s_index] and r : y_axis[r_index]
#indice_r = np.where((5.3 < y_axis) & (y_axis < 5.5))[0] ; print("\nindice r :", indice_r)
#indice_s = np.where((0.58 < x_axis) & (x_axis < 0.6))[0] ; print("\nindice s :", indice_s)
#print(heatmap_values[indice_r,indice_s])


#### Print an heatmap from heatmap_values ####
    
def print_heatmap(heatmap_values, bio_para, num_para, model_para, rlog, precision) :    
    r, s, h, a, difWW, difDW, difDD, c, conversion_timing = bio_para
    x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta = heatmap_axes(rlog, precision)
    
    # Figure
    fig, ax = plt.subplots()     
    
    # Size of a simulation pixel, where the distance 'center of the first pixel' to 'center of the last pixel' is 1. (useful when rlog=False) 
    pixel = 1/(precision-1)   
    # Ticks positions : we want the ticks to start from the bordure, not the center
    ax.set_xticks(np.linspace(0-pixel/2,1+pixel/2,len(np.arange(x_min,x_max+0.1,0.1)))*(precision-1))                  
    if rlog == True : ax.set_yticks(np.linspace(0,1,4)*(precision-1))    
    elif rlog == False : ax.set_yticks(np.linspace(0-pixel/2,1+pixel/2,len(np.arange(int(y_min),y_max+1,1)))*(precision-1)) 
    else : ax.set_yticks(np.linspace(0-pixel/2,1+pixel/2,len(np.arange(x_min,x_max+0.1,0.1)))*(precision-1))                                        
    # Ticks labels
    ax.set_xticklabels(np.around(np.arange(x_min,x_max+0.1,0.1),2))                                         
    if rlog == True : ax.set_yticklabels(np.around(np.logspace(-2, 1, num=4),2))  
    elif rlog == False : ax.set_yticklabels(np.around(np.arange(int(y_min),y_max+1,1),2))   
    else : ax.set_yticklabels(np.around(np.arange(x_min,x_max+0.1,0.1),2))   
     
    # Colorbar creation
    colors1 = plt.cm.Blues(np.flip(np.linspace(0.3, 1, 128))); colors2 = plt.cm.hot(np.flip(np.linspace(0, 0.75, 128)))  # Create the color scale   
    colors = np.vstack((colors1, colors2)) # Merge the two
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)  
    # Plot heatmap values
    im = ax.imshow(heatmap_values,cmap=mymap, vmin=num_para[-1][0], vmax=num_para[-1][1], aspect='auto')  
    # Add the colorbar
    ax.figure.colorbar(im, ax=ax)
    
    # Add another x axis on the top of the heatmap to indicate s_1 and s_2
    #s_1 = c/(1-h*(1-c))
    #if conversion_timing == "zygote" : 
    #    s_2 = c/(2*c + h*(1-c))
    #if conversion_timing == "germline" :
    #    if h != 0 : s_2 = c/(2*c*h + h*(1-c))  
    #    else : s_2 = 100                      # inf in reality
    #axtop = ax.twiny(); ticks = []; ticklabels = []
    #if x_min < s_1 and s_1 < x_max : 
        #ticks.append((s_1-x_min)/(x_max-x_min)); ticklabels.append("s1")
        # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
        #ax.vlines((s_1-x_axis[0])*(1/x_delta),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))        
    #if x_min < s_2 and s_2 < x_max : 
        #ticks.append((s_2-x_min)/(x_max-x_min)); ticklabels.append("s2")
        # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
        #ax.vlines((s_2-x_axis[0])*(1/x_delta),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))
    #axtop.set_xticks(ticks) 
    #axtop.set_xticklabels(ticklabels)
        
    # Plot lines    
    #abscisse = np.arange(precision)    
    # Tool to draw precise line
    #if rlog : y_axis_precise = np.logspace(-2, 1, num=precision*100)
    #else : y_axis_precise = np.arange(y_min+y_delta/2, y_max+y_delta/2, y_delta/100)
        
    # Pure drive persistance line
    #if rlog :
    #    eradication_drive = np.zeros(precision)
    #    for i in range(precision):
    #        s = x_axis[i]
    #        eradication_drive[i] = np.where(s/(1-s) < y_axis_precise)[0][0]/100         
    #else : 
    #    eradication_drive = (x_axis/(1-x_axis)-y_min)*((precision-1)/(y_max-y_min)) 
    #ax.plot(abscisse,eradication_drive-0.5, color='#73c946ff', label="eradication drive", linewidth = 4)           
                 
    # Composite persistance line
    #if (conversion_timing == "zygote" and (1-h)*(1-c) > 0.5) or (conversion_timing == "germline" and h < 0.5) :                        
        # neglect the values outside s_1 s_2 (which do not interest us)
     #   s1_s2_len = 200
     #   s1_s2_range = np.linspace(s_1,s_2,s1_s2_len) 
             
        #if conversion_timing == "zygote" :    
        #    p_star = (s1_s2_range*(1-(1-c)*(1-h)) - c*(1-s1_s2_range))/(s1_s2_range*(1-2*(1-c)*(1-h)))  
        #    mean_fitness = (1-s1_s2_range)*p_star**2+2*(c*(1-s1_s2_range)+(1-c)*(1-s1_s2_range*h))*p_star*(1-p_star)+(1-p_star)**2 
        #if conversion_timing == "germline" :           
        #    p_star = ((1-s1_s2_range*h)*(1+c)-1)/(s1_s2_range*(1-2*h))
        #    mean_fitness = (1-s1_s2_range)*p_star**2+2*(1-s1_s2_range*h)*p_star*(1-p_star)+(1-p_star)**2   
      
        #eradication_pop = np.zeros(s1_s2_len)
        #for i in range(s1_s2_len):
        #    eradication_pop[i] = np.where((1-mean_fitness[i])/mean_fitness[i] < y_axis_precise)[0][0]/100 
        #eradication_pop = eradication_pop[0:np.where(eradication_pop==0)[0][0]+1]
        #abscisse_pop = ((s1_s2_range-x_axis[0])*((precision-1)/(x_axis[-1]-x_axis[0])))[0:np.where(eradication_pop==0)[0][0] + 1]
        #ax.plot(abscisse_pop, eradication_pop-0.5, color='#40720cff', linewidth = 4) 

      
    nb_precise = 1001
    abscisse_precise = np.linspace(0,precision, nb_precise)
    x_axis_precise = np.linspace(x_min, x_max, nb_precise)  
    y_axis_precise = np.linspace(y_min, y_max, nb_precise)  
    s1_line = np.ones(nb_precise)*1.1; s2_line = np.ones(nb_precise)*1.1
    s1_line[0] = c; s1_line[-1] = 1
    if conversion_timing == "zygote" : s2_line[0] = 1/2; s2_line[-1] = c/(1+c)
    if conversion_timing == "germline" : s2_line[-1] = c/(1+c)
    for i in range(1,nb_precise-1):
        h_loc = x_axis_precise[i]
        s1_line[i] = np.where(c/(1-h_loc*(1-c)) < y_axis_precise)[0][0]/nb_precise   
        if conversion_timing == "zygote" : 
                s2_line[i] = np.where(c/(2*c + h_loc*(1-c)) < y_axis_precise)[0][0]/nb_precise 
        if conversion_timing == "germline" :
            if c/(h_loc*2*c + h_loc*(1-c)) < 1:
                s2_line[i] = np.where(c/(h_loc*2*c + h_loc*(1-c)) < y_axis_precise)[0][0]/nb_precise 
                
    ax.plot(abscisse_precise-0.5, s1_line*(precision/(y_max-y_min))-0.5, label="s1", linewidth = 2)    
    ax.plot(abscisse_precise-0.5, s2_line*(precision/(y_max-y_min))-0.5, label="s2", linewidth = 2)    
                 
    # Axis and label sizes
    plt.gca().invert_yaxis()                 # inverse r axis (increasing values)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)

    # Always save figure
    save_fig_or_data(f"heatmaps/{conversion_timing}_c_{c}_r_{r}", fig, [], f"{conversion_timing}_c_{c}_r_{r}", bio_para, num_para, model_para)

    plt.show()

