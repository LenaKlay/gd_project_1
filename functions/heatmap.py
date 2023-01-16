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
    
def heatmap_axes(y, rlog, precision):
    x_min = 0; x_max=1
    # x_delta is the size of a simulation pixel (size mesured with the s scale)  
    x_delta = (x_max-x_min)/precision    
    # x values for simulations (NB : x_min and x_max are not simulated, we simulate values centered on the simulation pixels)
    x_axis = np.arange(x_min+x_delta/2,x_max+x_delta/2,x_delta)      
    if y == "r" :                    
        if rlog : 
            y_min=0.01; y_max=10; y_delta = None
            y_axis = np.logspace(-2, 1, num=precision)
        else :
            y_min=0; y_max=12
            y_delta = (y_max-y_min)/precision  
            y_axis = np.arange(y_min+y_delta/2,y_max+y_delta/2,y_delta)               
    else :
        y_min = 0; y_max=1
        y_delta = (y_max-y_min)/precision    
        y_axis = np.arange(y_min+y_delta/2,y_max+y_delta/2,y_delta)       
    return(x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta)



#### Create an heatmap (speed of the traveling wave for couples (s,r) ####

def heatmap(bio_para, num_para, graph_para, rlog, precision, x, y):
    r,s,h,difWW,difDW,difDD,c,conversion_timing,cas,a,growth_dynamic,death_dynamic = bio_para
    x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta = heatmap_axes(y, rlog, precision)
    
    # Create a directory and save parameters.txt with y_axis and x_axis
    rshc = [['r','s','h','c'], [0,1,2,6]] 
    for i in range(4) : 
        if x == rshc[0][i] : xnum = rshc[1][i]; bio_para[xnum] = x_axis 
        if y == rshc[0][i] : ynum = rshc[1][i]; bio_para[ynum] = y_axis
  
    # cstnum = which values are constant in r,s,h,c (given by their index in bio_para)
    rshc[0].remove(x); rshc[0].remove(y); rshc[1].remove(xnum); rshc[1].remove(ynum)
    save_fig_or_data(f"heatmaps/{conversion_timing}_{y}_fct_{x}/{rshc[0][0]}_{bio_para[rshc[1][0]]}_{rshc[0][1]}_{bio_para[rshc[1][1]]}", [], [], None, bio_para, num_para)
     
    # Print principal parameters
    print("conversion_timing =", conversion_timing); print("\n", x, "=", x_axis);  print(y, "=", y_axis); print("h =", h); print("c =", c)
            
    # Initialization
    # matrix containing all the wave speeds for the different values of s and r 
    heatmap_values = np.zeros((precision,precision))      

    # Loop to calculate the speed, for each couple (s,r)
    for y_index in range(0, precision) : 
        print("\n------- NEW", y, "=", y_axis[y_index])
        for x_index in range(0, precision) :
            print("\n", x, "=", np.round(x_axis[x_index],3))             
            # Update values in bio_para
            bio_para[xnum] = y_axis[y_index]
            bio_para[ynum] = x_axis[x_index]                    
            # Speed value given by evolution.py
            heatmap_values[y_index,x_index] = evolution(bio_para, num_para, graph_para)[4][-1] 
            print("speed :", heatmap_values[y_index,x_index])
                     
        # for each r, save the corresponding line
        save_fig_or_data(f"heatmaps/{conversion_timing}_{y}_fct_{x}", [], heatmap_values[y_index,:], f"y_line_{y_index}", bio_para, num_para)   
        
    # Save final data
    save_fig_or_data(f"heatmaps/{conversion_timing}_{y}_fct_{x}", [], heatmap_values, f"{conversion_timing}_{y}_fct_{x}", bio_para, num_para)
    return(heatmap_values) 
 
    


       
#### Load an already existing heatmap ####
    
def load_heatmap(conversion_timing, c, r, h, s, rlog, precision, migale, cas) : 
    heatmap_values = np.zeros((precision,precision))
    coex_values = np.ones((precision,precision))*(-2)
    # migale = True means we download each line of the heatmap (from the migale file)
    if migale :
        for r_index in range(0,precision):
            if os.path.exists(f'../migale/heatmaps/speed_{r_index+1}.txt') :
                heatmap_values[r_index,:] = np.loadtxt(f'../migale/heatmaps/speed_{r_index+1}.txt')
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
    
def print_heatmap(heatmap_values, bio_para, num_para, rlog, precision, x, y, file_name) :    
    r,s,h,difWW,difDW,difDD,c,conversion_timing,cas,a,growth_dynamic,death_dynamic = bio_para
    x_axis, x_min, x_max, x_delta, y_axis, y_min, y_max, y_delta = heatmap_axes(y, rlog, precision)
    
    # Figure
    fig, ax = plt.subplots()     
    
    # Size of a simulation pixel, where the distance 'center of the first pixel' to 'center of the last pixel' is 1. (useful when rlog=False) 
    pixel = 1/(precision-1)   
    # Ticks positions : we want the ticks to start from the bordure, not the center
    ax.set_xticks(np.linspace(0-pixel/2,1+pixel/2,len(np.arange(x_min,x_max+0.1,0.1)))*(precision-1))       
    if y == "r" :            
        if rlog : ax.set_yticks(np.linspace(0,1,4)*(precision-1))    
        else : ax.set_yticks(np.linspace(0-pixel/2,1+pixel/2,len(np.arange(int(y_min),y_max+1,1)))*(precision-1)) 
    else : ax.set_yticks(np.linspace(0-pixel/2,1+pixel/2,len(np.arange(x_min,x_max+0.1,0.1)))*(precision-1))                                        
    # Ticks labels
    ax.set_xticklabels(np.around(np.arange(x_min,x_max+0.1,0.1),2))      
    if y == "r" :                                      
        if rlog : ax.set_yticklabels(np.around(np.logspace(-2, 1, num=4),2))  
        else : ax.set_yticklabels(np.around(np.arange(int(y_min),y_max+1,1),2))   
    else : ax.set_yticklabels(np.around(np.arange(x_min,x_max+0.1,0.1),2))   
     
    # Colorbar creation
    colors1 = plt.cm.Blues(np.flip(np.linspace(0.3, 1, 128))); colors2 = plt.cm.hot(np.flip(np.linspace(0, 0.75, 128)))  # Create the color scale   
    colors = np.vstack((colors1, colors2)) # Merge the two
    mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)  
    # Plot heatmap values
    im = ax.imshow(heatmap_values,cmap=mymap, vmin=num_para[-1][0], vmax=num_para[-1][1], aspect='auto')  
    # Add the colorbar
    ax.figure.colorbar(im, ax=ax)
    
        
    # Plot lines    
    nb_precise = 1001
    abscisse_precise = np.linspace(0,precision, nb_precise)
    x_axis_precise = np.linspace(x_min, x_max, nb_precise)  
    if y == "r" and rlog : y_axis_precise = np.logspace(-2, 1, num=nb_precise)
    else : y_axis_precise = np.linspace(y_min, y_max, nb_precise)
       
    if y == "r" and x == "h" :
        # Vertical axis h_1 and h_2
        h_1 = (s-c)/(s*(1-c))
        if conversion_timing == "zygote" : 
            h_2 = c*(1-2*s)/(s*(1-c))
        if conversion_timing == "germline" :
            if h != 0 : h_2 = c/(s*(1+c))  
        axtop = ax.twiny(); ticks = []; ticklabels = []
        if x_min < h_1 and h_1 < x_max : 
            ticks.append((h_1-x_min)/(x_max-x_min)); ticklabels.append("h1")
            # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
            ax.vlines((h_1-x_axis[0])*(1/x_delta),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))        
        if x_min < h_2 and h_2 < x_max : 
            ticks.append((h_2-x_min)/(x_max-x_min)); ticklabels.append("h2")
            # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
            ax.vlines((h_2-x_axis[0])*(1/x_delta),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))
        axtop.set_xticks(ticks) 
        axtop.set_xticklabels(ticklabels)
             
    
    
    if y == "r" and x == "s" :     
        # Vertical axis s_1 and s_2
        s_1 = c/(1-h*(1-c))
        if conversion_timing == "zygote" : 
            s_2 = c/(2*c + h*(1-c))
        if conversion_timing == "germline" :
            if h != 0 : s_2 = c/(2*c*h + h*(1-c))  
            else : s_2 = 100                      # inf in reality
        axtop = ax.twiny(); ticks = []; ticklabels = []
        if x_min < s_1 and s_1 < x_max : 
            ticks.append((s_1-x_min)/(x_max-x_min)); ticklabels.append("s1")
            # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
            ax.vlines((s_1-x_axis[0])*(1/x_delta),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))        
        if x_min < s_2 and s_2 < x_max : 
            ticks.append((s_2-x_min)/(x_max-x_min)); ticklabels.append("s2")
            # - 0.5 to correct the biais in the heatmap (0 is centered in the middle of the first pixel)
            ax.vlines((s_2-x_axis[0])*(1/x_delta),-0.5,precision-0.5, color="black", linewidth = 1, linestyle=(0, (3, 5, 1, 5)))
        axtop.set_xticks(ticks) 
        axtop.set_xticklabels(ticklabels)
            

        # Pure drive persistance line        
        # index for which we have a pure drive persistance line in between y_min and y_max
        if cas in ['a','b_pos','b_neg','c'] :            
            index_eradication_drive = np.intersect1d(np.where(x_axis_precise[:-1]/(1-x_axis_precise[:-1]) >= ((a-1)/2)**2*y_min)[0], np.where(x_axis_precise[:-1]/(1-x_axis_precise[:-1]) <= ((a-1)/2)**2*y_max)[0])
        elif cas in ['d'] :
            index_eradication_drive = np.intersect1d(np.where(x_axis_precise[:-1]/(a-x_axis_precise[:-1]) >= y_min)[0], np.where(x_axis_precise[:-1]/(a-x_axis_precise[:-1]) <= y_max)[0])
        # values of the pure drive persistance line
        eradication_drive = np.ones(len(index_eradication_drive))*(-1)
        if cas in ['a','b_pos','b_neg','c'] :   
            for i in range(len(eradication_drive)):
                # s_loc = local value of s (inside the for loop)
                s_loc = x_axis_precise[index_eradication_drive[i]]
                eradication_drive[i] = np.where(s_loc/(1-s_loc) < ((a-1)/2)**2*y_axis_precise)[0][0]/nb_precise  
        elif cas in ['d'] :
            for i in range(len(eradication_drive)):
                # s_loc = local value of s (inside the for loop)
                s_loc = x_axis_precise[index_eradication_drive[i]]
                if a>s_loc : 
                    eradication_drive[i] = np.where(s_loc < (a-s_loc)*y_axis_precise)[0][0]/nb_precise 
                else :
                    eradication_drive[i] = np.where(s_loc > (a-s_loc)*y_axis_precise)[0][0]/nb_precise 
        ax.plot(abscisse_precise[index_eradication_drive]-0.5, eradication_drive*precision-0.5, color='#73c946ff', label="eradication drive", linewidth = 4)           
            
                    
        # Composite persistance line
        index_s1s2 = np.arange(np.where(x_axis_precise>=s_1)[0][0], np.where(x_axis_precise<=s_2)[0][-1]+1)
        if (conversion_timing == "zygote" and (1-h)*(1-c) > 0.5) or (conversion_timing == "germline" and h < 0.5) :                               
            if conversion_timing == "zygote" :    
                p_star = (x_axis_precise[index_s1s2]*(1-(1-c)*(1-h)) - c*(1-x_axis_precise[index_s1s2]))/(x_axis_precise[index_s1s2]*(1-2*(1-c)*(1-h)))  
                mean_fitness = (1-x_axis_precise[index_s1s2])*p_star**2+2*(c*(1-x_axis_precise[index_s1s2])+(1-c)*(1-x_axis_precise[index_s1s2]*h))*p_star*(1-p_star)+(1-p_star)**2 
            if conversion_timing == "germline" :           
                p_star = ((1-x_axis_precise[index_s1s2]*h)*(1+c)-1)/(x_axis_precise[index_s1s2]*(1-2*h))
                mean_fitness = (1-x_axis_precise[index_s1s2])*p_star**2+2*(1-x_axis_precise[index_s1s2]*h)*p_star*(1-p_star)+(1-p_star)**2   
            # index for which we have a composite persistance line in between y_min and y_max
            index_eradication_pop = index_s1s2[np.intersect1d(np.where(mean_fitness/(1-mean_fitness) <= y_max)[0], np.where(y_min <= mean_fitness/(1-mean_fitness))[0])]    
            # values of the composite persistance line
            eradication_pop = np.ones(len(index_eradication_pop))*(-1)      
            for i in range(len(eradication_pop)):
                # m_loc = local value of the mean fitness (inside the for loop)
                m_loc = x_axis_precise[index_eradication_pop[i]]  
                eradication_pop[i] = np.where((1-m_loc)/m_loc < y_axis_precise)[0][0]/nb_precise 
            ax.plot(abscisse_precise[index_eradication_pop]-0.5, eradication_pop*precision-0.5, color='#40720cff', linewidth = 4) 

    if x == "h" and y == "s":
        
        # S1 line
        # index for which s1 is in between y_min and y_max          
        index_s1 = np.intersect1d(np.where(c/(1-x_axis_precise*(1-c)) <= y_max)[0], np.where(c/(1-x_axis_precise*(1-c)) >= y_min)[0])         
        # values of s1
        s1_line = np.ones(len(index_s1))*(-1); s1_line[0] = c; s1_line[-1] = 1
        for i in range(1,len(s1_line)-1):
            # h_loc = local value of h (inside the for loop)
            h_loc = x_axis_precise[index_s1[i]]
            s1_line[i] = np.where(c/(1-h_loc*(1-c)) < y_axis_precise)[0][0]/nb_precise   
        ax.plot(abscisse_precise[index_s1]-0.5, s1_line*precision-0.5, label="s1", linewidth = 2)    
        
        # S2 line
        if conversion_timing == "zygote" :               
                # index for which s2 is in between y_min and y_max     
                index_s2 = np.intersect1d(np.where(c/(2*c + x_axis_precise*(1-c)) <= y_max)[0], np.where(c/(2*c + x_axis_precise*(1-c)) >= y_min)[0])  
                # values of s2
                s2_line = np.ones(len(index_s2))*(-1); s2_line[0] = 1/2; s2_line[-1] = c/(1+c)            
                for i in range(1,len(s2_line)-1):
                    # h_loc = local value of h (inside the for loop)
                    h_loc = x_axis_precise[index_s2[i]]
                    s2_line[i] = np.where(c/(2*c + h_loc*(1-c)) < y_axis_precise)[0][0]/nb_precise  
        if conversion_timing == "germline" :
                # index for which s2 is in between y_min and y_max (+1 because of x_axis_precise[1:], otherwise div. by 0)     
                index_s2 = 1+np.intersect1d(np.where(c/(x_axis_precise[1:]*2*c + x_axis_precise[1:]*(1-c)) <= y_max)[0], np.where(c/(x_axis_precise[1:]*2*c + x_axis_precise[1:]*(1-c)) >= y_min)[0])  
                # values of s2
                s2_line = np.ones(len(index_s2))*(-1); s2_line[-1] = c/(1+c)
                for i in range(len(s2_line)-1):
                    # h_loc = local value of h (inside the for loop)
                    h_loc = x_axis_precise[index_s2[i]]
                    s2_line[i] = np.where(c/(h_loc*2*c + h_loc*(1-c)) < y_axis_precise)[0][0]/nb_precise  
        ax.plot(abscisse_precise[index_s2]-0.5, s2_line*precision-0.5, label="s2", linewidth = 2)    
        
        

    # Axis and label sizes
    plt.gca().invert_yaxis()                 # inverse r axis (increasing values)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)
    ax.set(xlabel=x, ylabel=y)

    # Always save figure
    save_fig_or_data(f"heatmaps", fig, [], f"{conversion_timing}_c_{c}_r_{r}_{file_name}", bio_para, num_para)

    plt.show()




