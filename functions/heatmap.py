#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 15:32:47 2021

@author: lena
"""


# Libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 

# External functions 
from evolution import evolution
from tanaka import tanaka
from graph import save_fig_or_data


def heatmap(heatmap_type, heatmap_para, mod, bio_para, model_para, num_para, graph_para, what_to_do):
    # Parameters
    precision, smin, smax, rmin, rmax = heatmap_para 
    r,s,h,a,difW,difH,difD,c,homing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    T,L,M,N,mod,theta = num_para
    graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion, show_graph_ini = graph_para
    
    # Range for r and s
    delta_s = (smax-smin)/precision    # delta_s is the size of a simulation pixel (size mesured with the s scale)  
    s_range = np.arange(smin+delta_s/2,smax+delta_s/2,delta_s)       # s values for simulations (NB : smin and smax are not simulated, we simulate values centered on the simulation pixels)     
    delta_r = (rmax-rmin)/precision    # delta_r is the size of a simulation pixel (size mesured with the r scale)  
    r_range = np.arange(rmin+delta_r/2,rmax+delta_r/2,delta_r)       # r values for simulations (NB : rmin and rmax are not simulated, we simulate values centered on the simulation pixels)
    
    # Create a directory and save parameters.txt with r_range and s_range
    if save_fig :
        bio_para[0] = r_range; bio_para[1] = s_range
        save_fig_or_data(f"heatmap/{heatmap_type}/{bio_para[8]}/h_{h}_c_{c}/T_{T}_L_{L}_M_{M}", [], [], None, bio_para, num_para)
     
    # Print principal parameters
    print("\nwhat_to_do =", what_to_do); print("homing =", bio_para[8]);  print("\nr =", r_range); print("s =", s_range); print("h =", h); print("c =", c)
      
      
    # Initiate
    # matrix containing all the wave speeds for the different values of s and r 
    heatmap_values = np.zeros((precision,precision))      
    # separate the graph in two areas : positiv speed (right side) and negativ speed (left side)
    zero_line = np.matrix([[],[]])                        

    for r_index in range(0, precision) : 
        print("\n------- NEW r=", r_range[r_index])
        # Denote that we already meet (or not) the first pixel of the line with negative speed (to draw the zero line)
        zero_done = False
        for s_index in range(0, precision) :
            r = r_range[r_index] ; bio_para[0] = r
            s = s_range[s_index] ; bio_para[1] = s
            
            # Print main parameters
            print(f"\n s={np.round(s,2)},r={np.round(r,4)},c={np.round(c,2)},h={np.round(h,2)}")
            
            # Classical heatmap
            if heatmap_type == "classic" :
                # Speed value for evolution.py
                heatmap_values[r_index,s_index] = evolution(bio_para, model_para, num_para, graph_para, what_to_do)[0]
                # First pixel of the line with negative speed (to draw the zero line)
                if s_index != 0 and heatmap_values[r_index,s_index-1]*heatmap_values[r_index,s_index]<=0 and heatmap_values[r_index,s_index] != 0 and zero_done == False :
                    zero_line = np.append(zero_line,np.matrix([[s_index-0.5],[r_index]]), axis=1)
                    zero_done = True                     
            # Tanaka cubic speed value     
            if heatmap_type == "speed_cubic" : 
                heatmap_values[r_index,s_index] = np.abs(evolution(bio_para, model_para, num_para, graph_para, what_to_do)[0]-tanaka(s,"cubic",model_para,num_para,graph_para)[1])          
            # Tanaka fraction speed value     
            if heatmap_type == "speed_fraction" : 
                heatmap_values[r_index,s_index] = np.abs(evolution(bio_para, model_para, num_para, graph_para, what_to_do)[0]-tanaka(s,"fraction",model_para,num_para,graph_para)[1])            
            # Tanaka cubic r(1-n) value 
            if heatmap_type == "r_one_minus_n_cubic" : 
                speed_girardin, W, H, D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
                n = D+H+W
                heatmap_values[r_index,s_index] = np.max(abs(r*(1-n)))           
            # Tanaka fraction r(1-n) value 
            if heatmap_type == "r_one_minus_n_fraction" : 
                speed_girardin, W, H, D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
                p = D/(D+H+W); n = D+H+W
                heatmap_values[r_index,s_index] = np.max(abs(r*(1-n) - s*p*(2-p)/(1-s+s*(1-p)**2)))
                
            # Print each value of the heatmap once it is computed.    
            print(f"heatmap value ={heatmap_values[r_index,s_index]} \n") 
            
        # for each r, save the corresponding line
        if save_fig : 
            save_fig_or_data(f"heatmap/{heatmap_type}/{bio_para[8]}/h_{h}_c_{c}/T_{T}_L_{L}_M_{M}", [], heatmap_values[r_index,:], f"line_r_{np.round(r,2)}", bio_para, num_para)
           
    return(r_range, s_range, heatmap_values, zero_line) 
    
    
    
    
def print_heatmap(heatmap_values, zero_line, style, heatmap_type, heatmap_para, bio_para, num_para, save_fig) : 
    # Parameters
    precision, smin, smax, rmin, rmax = heatmap_para 
    r_range,s_range,h,a,difW,difH,difD,c,homing = bio_para
    T,L,M,N,mod,theta = num_para
        
    # Figure
    fig, ax = plt.subplots()     
    # Color choice
    if style == "simple" or style == "eradication" or style == "collapse" :
        # Negative scale
        colors1 = plt.cm.viridis(np.linspace(0, 0.8, 128))
        # Positive scale
        colors2 = plt.cm.plasma(np.flip(np.linspace(0., 1, 128)))
        # Merge the two
        colors = np.vstack((colors1, colors2))
        mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)
        # Plot heatmap values
        im = ax.imshow(heatmap_values,cmap=mymap, vmin=-4, vmax=4)
    else :
        # Plot heatmap values
        im = ax.imshow(heatmap_values,cmap='Blues',interpolation='bicubic')

    # Size of a simulation pixel, where the distance 'center of the first pixel' to 'center of the last pixel' is 1. 
    delta_square = 1/(precision-1)
    # Ticks positions : we want the ticks to start from the bordure, not the center
    ax.set_xticks(np.linspace(0-delta_square/2,1+delta_square/2,len(np.arange(smin,smax+0.1,0.1)))*(precision-1))     
    ax.set_yticks(np.linspace(0-delta_square/2,1+delta_square/2,len(np.arange(int(rmin),rmax+1,1)))*(precision-1))                    
    # Ticks labels
    ax.set_xticklabels(np.around(np.arange(smin,smax+0.1,0.1),2))                                    
    ax.set_yticklabels(np.around(np.arange(int(rmin),rmax+1,1),2))          
    
    # Colorbar
    ax.figure.colorbar(im, ax=ax)
    # Rotate the xtick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    # Problem if s = 1
    if np.isin(1, s_range) :                              # Avoid dividing by 0
        s_range = np.resize(s_range,len(s_range)-1); r_range = np.resize(r_range,len(r_range)-1)
        
    # Abscisse to plot lines    
    abscisse = np.arange(precision)       # (s_range-smin)*(precision-1)/(smax-smin)   
    
    # Lines always there
    if style == "simple" or style == "eradication" or style == "collapse" :
        # Zero line
        if np.shape(zero_line)[1] != 0 : 
           ax.plot(np.array(zero_line[0,:]).ravel(),np.array(zero_line[1,:]).ravel(),color="red",label="zero speed", linewidth = 1.5, linestyle='-.')
        # Eradication line
        eradication_drive = (s_range/(1-s_range)-rmin)*((precision-1)/(rmax-rmin))
        ax.plot(abscisse,eradication_drive, color='orangered',label="eradication drive", linewidth = 2)  
        ax.legend()
        
    # Eradication zone
    if style == "eradication":
        abscisse_for_zero_line = np.unique(np.array(zero_line[0,:]).ravel(),return_index=True)
        unique_zero_line = np.array(zero_line[1, abscisse_for_zero_line[1]]).ravel()
        abscisse_for_zero_line_into_s = abscisse_for_zero_line[0]*(smax-smin)/(precision-1)+smin
        eradication_drive_for_zero_line = (abscisse_for_zero_line_into_s/(1-abscisse_for_zero_line_into_s)-rmin)*((precision-1)/(rmax-rmin))
        # Fill part 1
        plt.fill_between(abscisse, eradication_drive, where=np.round(abscisse,2)<=zero_line[0,0], color='orangered')
        # Fill part 2
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
               
    # Set graph parameters
    ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
    ax.set_ylabel("r (growth rate)", fontsize=12)
    fig.tight_layout()
    plt.gca().invert_yaxis()                 # inverse r axis (increase)
    ax.xaxis.set_tick_params(labelsize=9)
    ax.yaxis.set_tick_params(labelsize=11)
    #ax.set_title(f"Speed of the wave (c={c} and h={h}, {homing})",fontsize = 13)
    #fig.suptitle(f"Heatmap : {heatmap_type}", fontsize=14)
    #plt.xlim([0,precision-1]);plt.ylim([0,precision-1])
    #plt.title("\n", f"Heatmap : {heatmap_type}")
    plt.show()
    
    # Save figures and datas 
    if save_fig :
        save_fig_or_data(f"heatmap/{heatmap_type}/{bio_para[8]}/h_{h}_c_{c}/T_{T}_L_{L}_M_{M}", fig, heatmap_values, f"{precision}_heatmap", bio_para, num_para)
        if heatmap_type == "classic" : 
            save_fig_or_data(f"heatmap/{heatmap_type}/{bio_para[8]}/h_{h}_c_{c}/T_{T}_L_{L}_M_{M}", [], zero_line, f"{precision}_zero_line", bio_para, num_para)


      
# To print loaded heatmaps :      
def upload_and_plot_heatmap(c, h, homing, style, heatmap_type, heatmap_para, bio_para, num_para,save_fig) : 
    # upload heatmap 
    heatmap_values = np.loadtxt(f'../outputs/heatmap/{heatmap_type}/{bio_para[8]}/h_{h}_c_{c}/heatmap_{heatmap_para[0]}.txt') 
    # upload zero_line
    if heatmap_type == "classic" :
        zero_line = np.loadtxt(f'../outputs/heatmap/{heatmap_type}/{bio_para[8]}/h_{h}_c_{c}/zero_line_{heatmap_para[0]}.txt')   
    # print heatmap  
    print_heatmap(heatmap_values, zero_line, style, heatmap_type, heatmap_para, bio_para, num_para, save_fig) 
    return(heatmap_values, zero_line)


# To load an already existing heatmap : 
load = False 
if load : 
    # style indicates which lines and zones to draw
    c=0.5; h=0.5; homing='germline'; style = 'simple'; save_fig = True;     
    # heatmap_type indicates if we plot the speed of the wave at r small or r very large.
    heatmap_type = 'classic'; precision = 50; smin = 0.1; smax = 0.9; rmin = 0 ; rmax = 12  
    # update new values in parameters vectors
    heatmap_para = [precision, smin, smax, rmin, rmax]; bio_para[8]=homing
    # calculate s and r ranges
    delta_s = (smax-smin)/precision; s_range = np.arange(smin+delta_s/2,smax+delta_s/2,delta_s); bio_para[0] = r_range         
    delta_r = (rmax-rmin)/precision; r_range = np.arange(rmin+delta_r/2,rmax+delta_r/2,delta_r); bio_para[1] = s_range  
    # load and plot the heatmap     
    heatmap_values, zero_line = upload_and_plot_heatmap(c, h, homing, style, heatmap_type, heatmap_para, bio_para, num_para, save_fig)

# The heatmap_values[r_index,s_index] correspond to the values s : s_range[s_index] and r : r_range[r_index]
#indice_r = np.where((5.3 < r_range) & (r_range < 5.5))[0] ; print("\nindice r :", indice_r)
#indice_s = np.where((0.87 < s_range) & (s_range < 0.89))[0] ; print("\nindice s :", indice_s)
#print(heatmap_values[indice_r,indice_s])







