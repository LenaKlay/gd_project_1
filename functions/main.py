#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""

# This is annoying : 
# for i in np.linspace(0.1,1,10) : 
# print(i)    
# different from :
# np.linspace(0.1,1,10)


######################### Librairies and functions ############################

# Load libraries
import numpy as np
import matplotlib.pyplot as plt

# Load functions
from evolution import evolution, evolution_2D
from tanaka import tanaka
from graph import save_fig_or_data
from pulled_pushed_wave import pulled_pushed_wave
from heatmap import load_heatmap, heatmap, print_heatmap        

# Change font to serif
plt.rcParams.update({'font.family':'serif'})


############################ What to do ? #####################################

# possible choices for what_to_do : 
# "evolution" "evolution 2D" "tanaka cubic" "tanaka fraction" "KPP" "pulled pushed" 
# "speed function of time" "speed function of s" "speed function of r" "heatmap"

what_to_do = "heatmap"


######################### General parameters ##################################


### General parameters

## Biological
r = 8                               # intrinsic growth rate
s = 0.2                               # fitness disadvantage for drive
h = 0.9                             # dominance coefficient
c = 0.85                            # conversion rate
conversion_timing = "germline"      # "zygote" or "germline"
# Eradication : r = 1, s = 0.52, h = 0.6, c = 0.85  (condition extinction drive only : s > r/(r+1))

# Particular growth/death terms
cas = "d_neg"
# detailled cases (a : allee effect coefficient)
if cas == "a" : growth_dynamic = "logistical"; death_dynamic = "constant"; a = -1  # a not taken into acount in case a (but I use a=-1 it in the persistent line, heatmap.py)
if cas == "b_pos": growth_dynamic = "allee_effect"; death_dynamic = "constant"; a = 0.2
if cas == "b_neg": growth_dynamic = "allee_effect"; death_dynamic = "constant"; a = -0.2
if cas == "c" : growth_dynamic = "constant"; death_dynamic = "logistical"; a = -1  # a not taken into acount in case c (but I use a=-1 it in the persistent line, heatmap.py)
if cas == "d_pos" : growth_dynamic = "constant"; death_dynamic = "allee_effect"; a = 0.2
if cas == "d_neg" : growth_dynamic = "constant"; death_dynamic = "allee_effect"; a = -0.2

# Diffusion
difWW = 1; difDW = 1; difDD = 1    # diffusion coefficient for resp. WW, WD or DD individuals

## Numerical
CI = "center"                      # Initial conditions : "center" for having the border in the center, "left" for having the border on the left
T = 500                            # final time
L = 2000                           # length of the spatial domain
M = T*12                           # number of time steps
N = L*2                            # number of spatial steps
theta = 0.5                        # discretization in space : theta = 0.5 for Crank Nicholson, theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  
    
## Save outputs
save_fig = True                    # Save the figures (.svg and .png) 



### Parameters specific for each what_to_do

## Evolution
graph_type = "Allele densities"                           # "Genotype densities", "Genotype frequencies", "Allele densities" or "Allele frequencies" (or None if we don't want any evolution graph fct of time)
show_graph_ini = True                                     # Show graph at time 0
show_graph_end = True                                     # Show graph at time T
wild = True; heterozygous = True; drive = True            # What to draw on the graph
grid = True                                               # A grid or not
semilogy = False                                          # semilogy = False : classical scale, semilogy = True : log scale for y
xlim = None                                               # x scale on the graph (xlim = None, means it's not specify)
mod = T//5                                               # Draw graph every ..mod.. time. Also used to know when tracking points in time graphics.
## Evolution 2D
CI_lenght = N//4
 
## Speed function of s
s_min = 0.32 ; s_max = 0.4   
s_nb_points = 10       # s values are taken in np.linspace(s_min,s_max,s_nb_points) 

## Speed function of r
r_min = 6 ; r_max = 10           
r_nb_points = 10       # r values are taken in np.linspace(r_min,r_max,r_nb_points) 

## Heatmap
x = "s"; y = "r"         # Heatmap axes
rlog = True          # r in log scale or not (when y = "r")
precision = 50       # Number of value on s and r scale for the heatmap
load = True          # do we load the datas (True) or create them (False)
migale = True        # if load == True, do the datas come from migale cluster, or from the folder "figures/heatmaps"
vmin = -8; vmax = 8      # Range min and max for the heatmap colour scale

### Small particularities
# Do not show graph for these tasks
if what_to_do in ["pulled pushed","speed function of time","speed function of s","speed function of r","heatmap"] :
    graph_type = None; show_graph_ini = False; show_graph_end = False    
# For r -> infinity, we need a big time precision (for the simulation to be stable)
if what_to_do in ["speed function of time","speed function of s","speed function of r"] :
    if r>=10 :  
        M = T*40; N = L*10
                     

#################### Group parameters for lisibility ##########################
        
bio_para = [r,s,h,difWW,difDW,difDD,c,conversion_timing,cas,a,growth_dynamic,death_dynamic]
num_para = [CI,T,L,M,N,theta,[vmin,vmax]]
graph_para = [wild, heterozygous, drive, mod, grid, semilogy, xlim, graph_type, show_graph_ini, show_graph_end, save_fig]

############################### Main program ##################################

## Evolution
if what_to_do == "evolution" :
    print("\nwhat_to_do =", what_to_do); print("conversion_timing =", conversion_timing); print("\nr =", r); print("s =", s); print("h =", h); print("c =", c)
    W, H, D, time, speed = evolution(bio_para, num_para, graph_para)
    print("\nspeed :", speed[-1])
    
    
## Evolution 2D
if what_to_do == "evolution 2D" :
    print("\nwhat_to_do =", what_to_do); print("conversion_timing =", conversion_timing); print("\nr =", r); print("s =", s); print("h =", h); print("c =", c)
    W, H, D, time, speed = evolution_2D(bio_para, num_para, graph_para, CI_lenght)
    print("\nspeed :", speed[-1])
    
    
## Pulled pushed wave (follow two neutral population in the wave; one at the front, the other elsewhere)
if what_to_do == "pulled pushed" : 
    pulled_pushed_wave(bio_para, num_para, graph_para)


## Tanaka cubic
if what_to_do == "tanaka cubic" :    
    print("\nwhat_to_do =", what_to_do); print("conversion_timing = zygote"); print("\ns =", s)    
    p_cubic, time_cubic, speed_cubic = tanaka(s,"cubic",num_para,graph_para)     
    
    
## Tanaka fraction  
if what_to_do == "tanaka fraction" :   
    print("\nwhat_to_do =", what_to_do); print("conversion_timing = zygote"); print("\ns =", s)   
    p_fraction, time_fraction, speed_fraction  = tanaka(s,"fraction",num_para,graph_para)
   
     
## KPP
if what_to_do == "KPP" :     
    print("\nwhat_to_do =", what_to_do)  
    p_KPP, time_KPP, speed_KPP = tanaka(s,"KPP",num_para,graph_para)



## Speed function of time 
if what_to_do == "speed function of time" :    
    print("\nwhat_to_do =", what_to_do); print("conversion_timing =", conversion_timing);  print("\nr =", r); print("s =", s); print("h =", h); print("c =", c)       
    # compute speed for Leo/Flo and Tanaka's models
    W, H, D, time_leoflo, speed_leoflo = evolution(bio_para, num_para, graph_para)
    p_cubic, time_cubic, speed_cubic = tanaka(s,"cubic",num_para,graph_para) 
    p_fraction, time_fraction, speed_fraction = tanaka(s,"fraction",num_para,graph_para)
    p_KPP, time_KPP, speed_KPP = tanaka(s,"KPP",num_para,graph_para)
    # figure
    fig, ax = plt.subplots()  
    ax.plot(time_cubic, speed_cubic, label="Tanaka cubic")
    ax.plot(range(0,T+1,mod), np.ones(int(T/mod+1))*(2-3*s)/np.sqrt(2*s), label="Tanaka cubic solution part.")
    ax.plot(time_fraction, speed_fraction, label="Tanaka fraction")
    if s < 0.5 : 
        ax.plot(time_KPP, speed_KPP, label="KPP numeric")
        ax.plot(np.arange(0,T+1,mod), np.ones(int(T/mod+1))*2*np.sqrt(1-2*s), label="KPP exact")         
    ax.plot(time_leoflo, speed_leoflo, label="Leo/Flo")
    ax.set(title = f"Speed function of time with s={s}")
    ax.set(xlabel='Time', ylabel='Speed')   
    ax.grid();plt.legend(); plt.show()      
    if save_fig : 
        save_fig_or_data(f"speed_fct_of_time/r_{r}_s_{s}_h_{h}_c_{c}", fig, [], "speeds_fct_of_time", bio_para, num_para)
        save_fig_or_data(f"speed_fct_of_time/r_{r}_s_{s}_h_{h}_c_{c}", [], speed_leoflo, "speed_flo_leo", bio_para, num_para)
        save_fig_or_data(f"speed_fct_of_time/r_{r}_s_{s}_h_{h}_c_{c}", [], speed_cubic, "speed_tanaka_cubic", bio_para, num_para)
        save_fig_or_data(f"speed_fct_of_time/r_{r}_s_{s}_h_{h}_c_{c}", [], speed_fraction, "speed_tanaka_fraction", bio_para, num_para)
        save_fig_or_data(f"speed_fct_of_time/r_{r}_s_{s}_h_{h}_c_{c}", [], speed_KPP, "speed_KPP", bio_para, num_para)
    
    
            
## Speed function of s
if what_to_do == "speed function of s" :
    # Comparison with the perfect conversion in the zygote model    
    r=0; conversion_timing = "zygote"; c=1; h=0    
    # Update parameters
    bio_para = [r,s,h,difWW,difDW,difDD,c,conversion_timing,cas,a,growth_dynamic,death_dynamic]
    # s values
    s_values = np.linspace(s_min,s_max,s_nb_points)
    fct_of_s = np.zeros((7,s_nb_points)) 
    fct_of_s[0,:] = s_values 
    # s values for s < 0.5
    list_s05 = []       
    # Print principal parameters
    print("\nwhat_to_do =", what_to_do); print("conversion_timing =", conversion_timing);  print("\nr =", r); print("s =", s_values); print("h =", h); print("c =", c, "\n")
    # Loop on s_values
    for s_index in range(len(s_values)) :
        # print and update s value
        bio_para[1] = np.round(s_values[s_index],3); print('\ns =',bio_para[1])     
        # second line = numerical speed Leo and Florence's model
        fct_of_s[1,s_index] = evolution(bio_para, num_para, graph_para)[4][-1]           
        # third line = exact bistable speed for Tanaka cubic 
        fct_of_s[2,s_index] = (2-3*s)/np.sqrt(2*s)     
        # fourth line = numerical speed of Tanaka fraction
        fct_of_s[3,s_index] = tanaka(s,"fraction",num_para,graph_para)[2][-1]  
        # fifth line = numerical speed of Tanaka cubic
        fct_of_s[4,s_index] = tanaka(s,"cubic",num_para,graph_para)[2][-1] 
        print(fct_of_s[3,s_index])
        if s <= 0.5 : 
            list_s05.append(s)
            # sixth line = numerical speed of KPP model (only defined when s < 0.5)
            fct_of_s[5,s_index] = tanaka(s,"KPP",num_para,graph_para)[2][-1]  
            print(fct_of_s[5,s_index])
            # seventh line = theorical speed of KPP model (only defined when s < 0.5)
            fct_of_s[6,s_index] = 2*np.sqrt(1-2*s)
            print(fct_of_s[6,s_index])
    # Figure
    fig, ax = plt.subplots()    
    ax.plot(fct_of_s[0,:],fct_of_s[2,:], label="(2-3*s)/sqrt(2*s)")    
    ax.plot(fct_of_s[0,:],fct_of_s[3,:], label="Model (12)")
    ax.plot(fct_of_s[0,:],fct_of_s[4,:], label="Barton")
    if r == 0 : ax.plot(list_s05,fct_of_s[1,0:len(list_s05)], label="Model (4)")
    else : print("r est different de 0.")
    ax.plot(list_s05,fct_of_s[5,0:len(list_s05)], label="Fisher KPP")
    ax.plot(list_s05,fct_of_s[6,0:len(list_s05)], label="2 sqrt(1-2s)")
    ax.grid(); ax.legend(); plt.show()
    if save_fig : 
        bio_para[1] = s_values 
        save_fig_or_data(f"speed_function_of_s/r_{r}_h_{h}_c_{c}", fig, fct_of_s, f"s_from_{s_min}_to_{s_max}", bio_para, num_para)
        
        

## Speed function of r
if what_to_do == "speed function of r" :
    # r values          
    r_values = np.linspace(r_min,r_max,r_nb_points)
    fct_of_r = np.zeros((2,len(r_values))) 
    fct_of_r[0,:] = r_values        
    # Print principal parameters
    print("\nwhat_to_do =", what_to_do); print("conversion_timing =", conversion_timing);  print("\nr =", r_values); print("s =", s); print("h =", h); print("c =", c)     
    # Loop on r_values 
    for r_index in range(len(r_values)) :
        # print and update s value
        bio_para[0] = np.round(r_values[r_index],3); print('\nr =',bio_para[0])# first line = numerical speed Leo and Florence's model
        fct_of_r[1,r_index] = evolution(bio_para, num_para, graph_para)[4][-1]     
    # Figure
    fig, ax = plt.subplots()    
    ax.plot(fct_of_r[0,:],fct_of_r[1,:], label="Leo/Flo")
    ax.plot(fct_of_r[0,:], 2*np.sqrt(1-2*s)*np.ones(len(r_values)), label="KPP exact")
    ax.grid(); ax.legend(); plt.show()
    if save_fig : 
        bio_para[0] = r_values         
        save_fig_or_data(f"speed_function_of_r/s_{s}_h_{h}_c_{c}", fig, fct_of_r, f"r_from_{r_min}_to_{r_max}", bio_para, num_para)


    
## Heatmaps
if what_to_do == "heatmap" :
        # Parameters
        T = 500; L = 2000; M = T*6; N = L; CI = 'center'
        # Update parameters
        num_para = [CI,T,L,M,N,theta,[vmin,vmax]] 
        bio_para = [r,s,h,difWW,difDW,difDD,c,conversion_timing,cas,a,growth_dynamic,death_dynamic]
        # Obtain the heatmap values                       
        if load : 
            heatmap_values, coex_values = load_heatmap(bio_para, rlog, precision, migale, x, y)
        else :
            heatmap_values = heatmap(bio_para, num_para, graph_para, rlog, precision, x, y)            
        # Print heatmap
        print_heatmap(heatmap_values, bio_para, num_para, rlog, precision, x, y, "fig_speed") 
        # Print coex       
        num_para[-1][0] = -1; num_para[-1][1] = 1
        print_heatmap(coex_values, bio_para, num_para, rlog, precision, x, y, "fig_coex") 
        # Superposition
        coex_values[np.where(coex_values == -1)] = -0.1
        num_para[-1][0] = -8; num_para[-1][1] = 8
        print_heatmap(heatmap_values+8*(coex_values+0.1), bio_para, num_para, rlog, precision, x, y, "fig_superposition") 

                





            
        

