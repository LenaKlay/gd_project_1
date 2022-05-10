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


############################### Load Stuff #######################################

# Load libraries
import numpy as np
import matplotlib.pyplot as plt


# Load functions
from evolution import evolution
from tanaka import tanaka
from graph import save_fig_or_data

# Change font to serif
plt.rcParams.update({'font.family':'serif'})


############################### Parameters ######################################

# Biological

# homing = "germline"
# s = 0.668
# h=0.5 et c=0.5
# r_range = array([ 0.12,  0.36,  0.6 ,  0.84,  (1.08),  (1.32),  1.56,  1.8 ,  2.04,
#        2.28,  2.52,  2.76,  3.  ,  3.24,  3.48,  3.72,  3.96,  4.2 ,
#        4.44,  4.68,  4.92,  5.16,  5.4 ,  5.64,  5.88,  6.12,  6.36,
#        6.6 ,  6.84,  7.08,  7.32,  7.56,  7.8 ,  8.04,  8.28,  8.52,
#        8.76,  9.  ,  9.24,  9.48,  9.72,  9.96, 10.2 , 10.44, 10.68,
#       10.92, 11.16, 11.4 , 11.64, 11.88])


# s = 0.876
# r_range = array([ 0.12,  0.36,  0.6 ,  0.84,  1.08,  1.32,  1.56,  1.8 ,  2.04,
#        2.28,  2.52,  2.76,  3.  ,  3.24,  3.48,  3.72,  3.96,  4.2 ,
#        4.44,  4.68,  4.92,  5.16,  (5.4 ,  5.64,  5.88,  6.12,  6.36),
#        6.6 ,  6.84,  7.08,  7.32,  7.56,  7.8 ,  8.04,  8.28,  8.52,
#        8.76,  9.  ,  9.24,  9.48,  9.72,  9.96, 10.2 , 10.44, 10.68,
#       10.92, 11.16, 11.4 , 11.64, 11.88])
    
r = 0    # growth rate
s = 0.38   # when selection only acts on survival with a fitness cost s (b=1 and d=1) 
h = 0.1     # and sh for heterozygous individuals
a = 0       # coefficient for allee effect (growth or death)

difW = 1   # diffusion coefficient for WW individuals
difH = 1   # diffusion coefficientrate for WD individuals
difD = 1   # diffusion coefficient rate for DD individuals

c = 0.75              # homing rate
homing = "zygote"   # "zygote" or "germline"

# Initialization

CI = "center"       # "center" for having the border in the center, "left" for having the border on the left

# How Population Grow and Decrease

growth_dynamic = "logistical"     # exponential or logistical (growth rate is r+1 for expon, r*coef+1 for logist)
death_dynamic = "exponential"      # exponential or logistical  (death rate is always 1)

# Linearization
  
linear_growth = False  
linear_mating = False          

# Numerical

T = 400       # final time
L = 800       # length of the spatial domain
M = T*8        # number of time steps
N = L*4          # number of spatial steps

theta = 0.5     # discretization in space : theta = 0.5 for Crank Nicholson
                # theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  
                   
# Graph

graph_type = "Allele proportions"                         # "Population densities", "Population proportions", "Allele densities" or "Allele proportions" (or None if we don't want any evolution graph fct of time)
show_graph_ini = True                                    # Show graph at time 0
show_graph_fin = True                                    # Show graph at time T
wild = True; heterozygous = True; drive = True            # What to draw on the graph
grid = True                                               # A grid or not
semilogy = False                                          # semilogy = False : classical scale, semilogy = True : log scale for y
xlim = None                                               # x scale on the graph (xlim = None, means it's not specify)
mod = int(T/4)                                            # Draw graph every ..mod.. time. Also used to know when tracking points in time graphics.
save_fig = True                                           # Save the figures (.pdf) 

# Speed calculus
speed_proportion = False            # True : use the wild-type number to compute the speed, False : use the wild-type proportion. 



# Group parameters for lisibility
bio_para = [r,s,h,a,difW,difH,difD,c,homing]
model_para = [CI,growth_dynamic,death_dynamic,linear_growth,linear_mating]
num_para = [T,L,M,N,mod,theta]
graph_para = [graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion,show_graph_ini,show_graph_fin]


############################ What to do ? #######################################

what_to_do = "heatmap"
# Bring the principal parameters together to make it easier.
# "evolution",  "tanaka cubic"  "tanaka fraction", "KPP" : simplest task, draw the propagation regarding the parameters above.
# "speed function of time" : idem + draw the speed as a function of time.
# "heatmap" : draw an heatmap function of s and r

############################### Main program #########################################

# Evolution
if what_to_do == "evolution" : 
        #c=np.linspace(0,1,10)
        #plt.plot(c, (2+3*c)**2-8*(1+c)*c)
        #det = (2+3*c)**2-8*(1+c)*c
        #plt.plot(c,((2+3*c)+np.sqrt((2+3*c)**2-8*(1+c)*c))/(4*(1+c)))
        #plt.plot(c,((2+3*c)-np.sqrt((2+3*c)**2-8*(1+c)*c))/(4*(1+c)))
        #plt.show()
        c=0.75
        h=0.1
        poly = (2*c + h*(1-c))*(2*(1-(1-c)*(1-h)))*s**2 +  s*(2*(1-c)*(1-h) -1 - 2*c - h*(1-c) -2*c*(1-(1-c)*(1-h))) + c 
        det = (2*(1-c)*(1-h) -1 - 2*c - h*(1-c) -2*c*(1-(1-c)*(1-h)))**2 - 4*(2*c + h*(1-c))*(2*(1-(1-c)*(1-h)))**c
        rac1 = (-(2*(1-c)*(1-h) -1 - 2*c - h*(1-c) -2*c*(1-(1-c)*(1-h))+np.sqrt(det))/(2*(2*c + h*(1-c))*(2*(1-(1-c)*(1-h)))))
        rac2 = (-(2*(1-c)*(1-h) -1 - 2*c - h*(1-c) -2*c*(1-(1-c)*(1-h))-np.sqrt(det))/(2*(2*c + h*(1-c))*(2*(1-(1-c)*(1-h)))))
        
        
        c=np.linspace(0.5,1,10)
        #plt.plot(c, (1-4*c-2*c**2)**2-16*c**3)
        #det = (2+3*c)**2-8*(1+c)*c
        plt.plot(c,(-(1-4*c-2*c**2)+np.sqrt( (1-4*c-2*c**2)**2-16*c**3 ))/(8*c**2))
        plt.plot(c,(-(1-4*c-2*c**2)-np.sqrt( (1-4*c-2*c**2)**2-16*c**3 ))/(8*c**2))
        plt.plot(c,np.ones(10))
        plt.show()
        
        c=0.75
        h=np.linspace(0,1,10)
        plt.plot(h,c/(1-h*(1-c)))
        plt.plot(h,c/(2*c + h*(1-c)))
        plt.show()
        s_1 = c/(1-h*(1-c))   
        if homing == "zygote" :
            s_3 = c/(2*c + h*(1-c))
            if s_1 < s_3 : 
                print("\nCoexistence", ":", np.round(s_1,3), "< s <", np.round(s_3,3))
            else :  
                print("\nBistability", ":", np.round(s_3,3), "< s <", np.round(s_1,3))
            if c*(1-2*s)-(1-c)*s*h > 0 : 
                print("Linear speed :", 2*np.sqrt(c*(1-2*s)-(1-c)*s*h))
                print("p_star :", (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h))))                
        if homing == "germline" : 
            s_2 = c/(2*c*h + h*(1-c))
            if s_1 < s_2 :
                print("\nCoexistence", ":", np.round(s_1,3) , "< s <", np.round(s_2,3))
            else :  
                print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
            if c*(1-2*s*h)-(1-c)*s*h > 0 : 
                print("Linear speed :", 2*np.sqrt(c*(1-2*s*h)-(1-c)*s*h))
                print("p_star :", ((1-s*h)*(1+c)-1)/(s*(1-2*h)))  

        print("\nwhat_to_do =", what_to_do); print("homing =", homing); print("\nr =", r); print("s =", s); print("h =", h); print("c =", c)
        W, H, D, time, speed = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
        # save the speed value
        #if save_fig :
        #    np.savetxt(f"../outputs/evolution/{homing}/s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}/r_{np.round(r,3)}/0_speed.txt", speed)        
        print("\nspeed :", speed[-1])
        
        # proportion
        if homing == "zygote" : 
            print("final prop :", ((D+0.5*H)/(W+H+D))[1000])
        if homing == "germline" : 
            print("final prop :", ((D+0.5*(1+c)*H)/(W+H+D))[1000], "\n")
        

        
        


# Tanaka cubic
if what_to_do == "tanaka cubic" : 
    
    print("\nwhat_to_do =", what_to_do); print("homing = zygote"); print("\ns =", s)
    
    p_cubic, time_cubic, speed_cubic = tanaka(s,"cubic",model_para,num_para,graph_para) 
    
    
# Tanaka fraction  
if what_to_do == "tanaka fraction" :
    
    print("\nwhat_to_do =", what_to_do); print("homing = zygote"); print("\ns =", s)
    
    p_fraction, time_fraction, speed_fraction  = tanaka(s,"fraction",model_para,num_para,graph_para)
    
    
# KPP
if what_to_do == "KPP" :  
    
    print("\nwhat_to_do =", what_to_do)
    
    p_KPP, time_KPP, speed_KPP = tanaka(s,"KPP",model_para,num_para,graph_para)



        
        
# Speed function of time 
if what_to_do == "speed function of time" :
    
    print("\nwhat_to_do =", what_to_do); print("homing =", homing);  print("\nr =", r); print("s =", s); print("h =", h); print("c =", c)
   
    # No graph at each time
    graph_para[0] = None; graph_para[-1] = False; graph_para[-2] = False  
    
    # Parameters
    if r<10 : 
        T = 100; L = 400; M = T*8; N = L*4; mod = int(T/100) 
    if r>=10 :  # for r -> infinity, we need a big time precision (for the simulation to be stable)
        T = 100; L = 400; M = T*40; N = L*4; mod = int(T/100) 
    
    # To compare speed, we have to use the WT proportion front (because Tanaka doesn't deal with number)
    speed_proportion = True 
       
    for s in [0.8] : 
        # update s value
        bio_para[1] = s; print(s)

        # compute speed for Leo/Flo and Tanaka's models
        W, H, D, time_leoflo, speed_leoflo = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
        p_cubic, time_cubic, speed_cubic = tanaka(s,"cubic",model_para,num_para,graph_para) 
        p_fraction, time_fraction, speed_fraction = tanaka(s,"fraction",model_para,num_para,graph_para)
        p_KPP, time_KPP, speed_KPP = tanaka(s,"KPP",model_para,num_para,graph_para)
    
        print('s=',s)
        fig, ax = plt.subplots()  
        ax.plot(time_cubic, speed_cubic, label="Tanaka cubic")
        ax.plot(range(0,T+1,mod), np.ones(int(T/mod+1))*(2-3*s)/np.sqrt(2*s), label="Tanaka cubic solution part.")
        ax.plot(time_fraction, speed_fraction, label="Tanaka fraction")
        if s < 0.5 : 
            ax.plot(time_KPP, speed_KPP, label="KPP numeric")
            #ax.plot(range(0,T+1,mod), np.ones(int(T/mod+1))*2*np.sqrt(1-2*s), label="KPP exact")         
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
        
    
            
# Speed function of s
if what_to_do == "speed function of s" :
        
    # No graph at each time
    graph_para[0] = None; graph_para[-1] = False; graph_para[-2] = False  
      
    # Parameters
    if r<10 : 
        T = 1600; L = 3200; M = T*6; N = L
    if r>=10 :  # for r -> infinity, we need a big time precision (for the simulation to be stable)
        T = 100; L = 400; M = T*40; N = L*4
    # Comparison with the perfect conversion in the zygote model    
    r=0; homing = "zygote"; c=1; h=0
    
    # Update parameters
    num_para = [T,L,M,N,mod,theta]
    bio_para = [r,s,h,a,difW,difH,difD,c,homing]
     
    # Set the x-axis (values of s)
    s_min = 0.32 ; s_max = 0.42
    nb_points = 100
    s_values = np.linspace(s_min,s_max,nb_points)
    fct_of_s = np.zeros((7,len(s_values))) 
    # first line = s values
    fct_of_s[0,:] = s_values 
    # s values for s < 0.5
    list_s05 = []   
    
    # Print principal parameters
    print("\nwhat_to_do =", what_to_do); print("homing =", homing);  print("\nr =", r); print("s =", s_values); print("h =", h); print("c =", c)
      
    # Loop on s_values
    for s_index in range(len(s_values)) :
        # print and update s value
        s = np.round(s_values[s_index],3); print(s); bio_para[1] = s 
                
        # line = numerical speed of Tanaka cubic
        #fct_of_s[1,s_index] = tanaka(s,"cubic",model_para,num_para,graph_para)[2][-1] 
        # third line = exact bistable speed for Tanaka cubic 
        #fct_of_s[2,s_index] = (2-3*s)/np.sqrt(2*s)     
        # fourth line = numerical speed of Tanaka fraction
        fct_of_s[3,s_index] = tanaka(s,"fraction",model_para,num_para,graph_para)[2][-1]  
        # fifth line = numerical speed Leo and Florence's model
        #fct_of_s[4,s_index] = evolution(bio_para, model_para, num_para, graph_para, what_to_do)[4][-1]
        if s <= 0.5 : 
            list_s05.append(s)
            # sixth line = numerical speed of KPP model (only defined when s < 0.5)
            fct_of_s[5,s_index] = tanaka(s,"KPP",model_para,num_para,graph_para)[2][-1]  
            # seventh line = theorical speed of KPP model (only defined when s < 0.5)
            fct_of_s[6,s_index] = 2*np.sqrt(1-2*s)
           
            
    
    # Figure
    fig, ax = plt.subplots()    
    #ax.plot(fct_of_s[0,:],fct_of_s[1,:], label="Tanaka Cubic")
    #ax.plot(fct_of_s[0,:],fct_of_s[2,:], label="(2-3*s)/np.sqrt(2*s)")    
    ax.plot(fct_of_s[0,:],fct_of_s[3,:], label="Tanaka Fraction")
    #if r == 0 : ax.plot(list_s05,fct_of_s[4,0:len(list_s05)], label="Leo/Flo")
    #else : print("r est different de 0.")
    ax.plot(list_s05,fct_of_s[5,0:len(list_s05)], label="KPP numerique")
    ax.plot(list_s05,fct_of_s[6,0:len(list_s05)], label="KPP theorique")
    ax.grid(); ax.legend(); plt.show()
    if save_fig : 
        bio_para[1] = s_values 
        save_fig_or_data(f"speed_function_of_s/r_{r}_h_{h}_c_{c}", fig, fct_of_s, f"s_from_{s_min}_to_{s_max}", bio_para, num_para)
        
        





# Speed function of s
if what_to_do == "speed function of r" :
    
    # No graph at each time
    graph_para[0] = None; graph_para[-1] = False; graph_para[-2] = False      
    
    # Parameters
    T = 1000; L = 4000; M = T*6; N = L; mod = int(T/100)          
    
    # Update parameters
    num_para = [T,L,M,N,mod,theta]
     
    # Set the x-axis (values of s)
    r_min = 6 ; r_max = 10           # if r > 12, set M = T*40
    nb_points = 5
    r_values = np.linspace(r_min,r_max,nb_points)
    fct_of_r = np.zeros((2,len(r_values))) 
    # first line = r values
    fct_of_r[0,:] = r_values      
    
    # Print principal parameters
    print("\nwhat_to_do =", what_to_do); print("homing =", homing);  print("\nr =", r_values); print("s =", s); print("h =", h); print("c =", c)
      
    # Loop on r_values 
    for r_index in range(len(r_values)) :
        # print and update s value
        r = np.round(r_values[r_index],3); print(r); bio_para[0] = r 
        # first line = numerical speed Leo and Florence's model
        fct_of_r[1,r_index] = evolution(bio_para, model_para, num_para, graph_para, what_to_do)[4][-1]
        print(fct_of_r)
       
    # Figure
    fig, ax = plt.subplots()    
    ax.plot(fct_of_r[0,:],fct_of_r[1,:], label="Leo/Flo")
    ax.plot(fct_of_r[0,:], 2*np.sqrt(1-2*s)*np.ones(len(r_values)), label="KPP exact")
    ax.grid(); ax.legend(); plt.show()
    if save_fig : 
        bio_para[0] = r_values         
        save_fig_or_data(f"speed_function_of_r/s_{s}_h_{h}_c_{c}", fig, fct_of_r, f"r_from_{r_min}_to_{r_max}", bio_para, num_para)

        
      

if what_to_do == "heatmap" :
        from heatmap import heatmap
        from heatmap import print_heatmap
    
        heatmap_type = "classic" ; print("heatmap_type =", heatmap_type, "\n")   # "classic"  "speed_cubic" "speed_fraction" "r_one_minus_n_cubic"  "r_one_minus_n_fraction"   
        bio_para[8] = "zygote"   # Homing : "zygote" or "germline"                                   

        graph_para[0] = None; graph_para[-1] = False; graph_para[-2] = False          # No evolution graph
        precision = 50              # Number of value on s and r scale (including 0 and 1) for the heatmap

        if heatmap_type == "classic" :
            T = 1000; L = 4000; M = T*6; N = L 
            smin = 0.1; smax = 0.9; rmin = 0; rmax = 12  
                
            # update parameters
            heatmap_para = [precision, smin, smax, rmin, rmax]  
            num_para = [T,L,M,N,mod,theta]
           
            bio_para[0], bio_para[1], heatmap_values, zero_line = heatmap(heatmap_type, heatmap_para, mod, bio_para, model_para, num_para, graph_para, what_to_do)
            r_range = bio_para[0]; s_range = bio_para[1]
            print_heatmap(heatmap_values, zero_line, "simple", heatmap_type, heatmap_para, bio_para, num_para, save_fig)
            #print_heatmap(heatmap_values, zero_line, "eradication", heatmap_type, heatmap_para, bio_para, num_para, save_fig)
            #print_heatmap(heatmap_values, zero_line, "collapse", heatmap_type, heatmap_para, bio_para, num_para, save_fig)

        else :   
            T = 1000; L = 4000; M = T*40; N = L 
            smin = 0.1; smax = 0.9; rmin = 50 ; rmax = 60 
            
            # update parameters
            heatmap_para = [precision, smin, smax, rmin, rmax] 
            num_para = [T,L,M,N,mod,theta]
            
            bio_para[0], bio_para[1], heatmap_values, zero_line = heatmap(heatmap_type, heatmap_para, mod, bio_para, model_para, num_para, graph_para, what_to_do)                          
            r_range = bio_para[0]; s_range = bio_para[1]
            print_heatmap(heatmap_values, zero_line, None, heatmap_type, heatmap_para, bio_para, num_para, save_fig)

                






################# Not re-tested ######################"



if what_to_do == "system+evolution" :
    
    # External function
    from which_system import which_system 
    from graph import graph
    
    # No graph at each time
    graph_para[0] = None          
  
    # run evolution
    W,H,D,time,speed = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
          
    # Determine the position of the border and print the number of the system involve
    epsilon = 0.0001              # numerical error accepted
    window = 'minimum'            # where to look : 'zoom border' or 'everywhere'  or 'minimum'   
    border, system, system_bis = which_system(W,H,D,epsilon,window,num_para)     
  
   
    # Draw (and save if save_fig==True) the wave at the end (normal and log scale).
    directory = "system_evolution"
    file = "wave"
    title_fig = f"t = {T}"
    if speed != None and border != None:  
        X = np.linspace(0,N,N+1)*(L/N)
        semilogy = False
        xlim = (border-75,border+75)
        graph(X,W,H,D,T,graph_para,bio_para,num_para)
        semilogy = True 
        xlim = (0,L)       
        graph(X,W,H,D,T,graph_para,bio_para,num_para)
        
        

            
if what_to_do == "roots and speed function of s" : 
    
    s_min = 0.3 ; s_max = 0.8
    s_values = np.linspace(s_min,s_max,31)
    fct_of_s = np.zeros((7,len(s_values))) 
    fct_of_s[0,:] = s_values      # first line = s values
         
    for s_index in range(len(s_values)) :
            s = np.round(s_values[s_index],3); print(s)
            
            W, H, D, time, speed = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
            # Speed
            fct_of_s[1,s_index] = speed
            # Alpha
            alpha = (-(1-2*s)*(1-r) + np.sqrt((1-2*s)**2 * (1-r)**2 + 8 * r**2 * (1-s)))/(4 * r * (1-s))
            fct_of_s[2,s_index] = alpha
            # C_sty 
            fct_of_s[3,s_index] = (-speed-np.sqrt(speed**2 + 4*r*(2*(1-s)*alpha+1)))/2
            # D_sty 
            fct_of_s[4,s_index] = (-speed-np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2    
            # D_sty_plus 
            fct_of_s[5,s_index] = (-speed+np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2   
            # Condition 2 (a ou b ?)
            fct_of_s[6,s_index] = speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1) 
     
    np.savetxt(f'vitesse_en_fct_de_s.txt', fct_of_s)    
    fig, ax = plt.subplots()    
    ax.plot(fct_of_s[0,:], fct_of_s[1,:], label="vitesse", linestyle='--', color='black',linewidth = 1)
    ax.plot(fct_of_s[0,:], fct_of_s[3,:], label="C", color='blue')
    ax.plot(fct_of_s[0,:], fct_of_s[4,:], label="D-", color='deepskyblue')
    ax.plot(fct_of_s[0,:], fct_of_s[5,:], label="D+", color='deeppink')
    ax.plot(fct_of_s[0,:], np.zeros(len(s_values)), color='black', linewidth = 0.5)
    ax.set(xlabel='s (fitness disadvantage)')   
    ax.grid()
    ax.legend()
    plt.show()
    
    fig.savefig(f'vitesse_en_fct_de_s_model_{model}.pdf')
    
    fig, ax = plt.subplots()    
    #ax.plot(fct_of_s[0,:],fct_of_s[2,:], label="alpha", color="gold")
    ax.plot(fct_of_s[0,:],fct_of_s[6,:], label="condition 2", color="orange")
    ax.plot(fct_of_s[0,:],np.zeros(len(s_values)), color='black',linewidth = 0.5)
    ax.grid()
    ax.legend()
    plt.show()
        
       
            
        
if what_to_do == "limite r infini" :  
    
    T = 1000; L = 4000; M = T*40; N = L; mod=int(T/10)         #T = 150; L = 900; M = 2000; N = 3000  
    CI = "center"
    graph_type = "Population density" 
    linear_growth = False ; linear_mating = False              
    homing = "zygote"        # "zygote" or "germline" or "minimum"
    c = 1; h = 1
    s = 0.5

    # update parameters
    heatmap_para = [precision, smin, smax, rmin, rmax]  
    model_para = [CI,growth_dynamic,death_dynamic,linear_growth,linear_mating]
    num_para = [T,L,M,N,mod,theta]
    graph_para = [graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion]

    #nb_points = 15          # ATTENTION MODIF SALE POUR r=0
    #rmin = 50
    #rmax = 60
    nb_points=1; rmin=0; rmax=0; M = T*6
    
    other_case = 0
    
    if other_case == 1 : 
        growth_dynamic = "allee_effect"  
        death_dynamic = "exponential"    
    if other_case == 2 : 
        growth_dynamic = "exponential"  
        death_dynamic = "logistical"
    if other_case == 3 : 
        growth_dynamic = "exponential"  
        death_dynamic = "allee_effect" 
       
    cv_r1moinsn_fct_of_r = np.zeros((3,nb_points)) 
    cv_r1moinsn_fct_of_r[0,:] = np.linspace(rmin,rmax,nb_points)         # r values
    
    cv_speed_fct_of_r = np.zeros((5,nb_points)) 
    cv_speed_fct_of_r[0,:] = np.linspace(rmin,rmax,nb_points)            # r values
    cv_speed_fct_of_r[1,:] = np.ones(nb_points)*(2-3*s)/np.sqrt(2*s)     # theoritical speed for Tanaka cubic
    
    for i in range(0,nb_points) : 
        print(i)
        
        r = np.round(cv_r1moinsn_fct_of_r[0,i],6)
        W, H, D, time_leoflo, speed_leoflo = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
        p = D/(D+H+W); n = D+H+W
        
        # When r tends to infinity : does r(1-n) cv to 0 ? to  s*p*(2-p)/(1-s+s*(1-p)**2) ? 
        cv_r1moinsn_fct_of_r[1,i] = np.max(abs(r*(1-n) - s*p*(2-p)/(1-s+s*(1-p)**2)))
        cv_r1moinsn_fct_of_r[2,i] = np.max(abs(r*(1-n)))
        
        # When r tends to infinity : does the speed of the system Girardin 2021 tends to the speed Tanaka cubic ? Tanaka fraction ?
        cv_speed_fct_of_r[2,i] = tanaka(s,"cubic",model_para,num_para,graph_para)[2][-1] 
        cv_speed_fct_of_r[3,i] = tanaka(s,"fraction",model_para,num_para,graph_para)[2][-1]
        cv_speed_fct_of_r[4,i] = speed_leoflo
             
    np.savetxt(f'cv_r_un_moins_n_s_{s}.txt', cv_r1moinsn_fct_of_r) 
    np.savetxt(f'cv_speed_s_{s}.txt', cv_speed_fct_of_r) 
    
    fig, ax = plt.subplots()  
    ax.plot(cv_r1moinsn_fct_of_r[0,:], cv_r1moinsn_fct_of_r[1,:], label="cv Tanaka fraction ? ", color="orange")
    ax.plot(cv_r1moinsn_fct_of_r[0,:], cv_r1moinsn_fct_of_r[2,:], label="cv Tanaka cubic ?", color="green" )
    ax.set(xlabel='r (growth rate)', ylabel='tends to 0 if cv')   
    ax.grid();plt.legend();plt.show()
    fig.savefig(f'cv_r_un_moins_n_s_{s}.pdf')
    fig, ax = plt.subplots() 
    ax.plot(cv_speed_fct_of_r[0,:], cv_speed_fct_of_r[1,:], label="speed Tanaka cubic (theo)", color="blue")
    ax.plot(cv_speed_fct_of_r[0,:], cv_speed_fct_of_r[2,:], label="speed Tanaka cubic", color="green")
    ax.plot(cv_speed_fct_of_r[0,:], cv_speed_fct_of_r[3,:], label="speed Tanaka fraction", color="orange")
    ax.plot(cv_speed_fct_of_r[0,:], cv_speed_fct_of_r[4,:], label="speed Leo Florence", color="red")
    ax.set(xlabel='r (growth rate)', ylabel='speed')   
    ax.grid();plt.legend();plt.show()
    fig.savefig(f'cv_speed_s_{s}.pdf')



