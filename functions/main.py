#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 18:13:08 2021

@author: lena
"""

# This is annoying : 
#for i in np.linspace(0.1,1,10) : 
#  print(i)    
# different from :
#np.linspace(0.1,1,10)




############################### Load Stuff #######################################

# Load libraries
import numpy as np
import matplotlib.pyplot as plt


# Load functions
from evolution import evolution

# Change font to serif
plt.rcParams.update({'font.family':'serif'})


############################### Parameters ######################################

# Biological
    
r = 0.7   # growth rate
s = 0.49  # when selection only acts on survival with a fitness cost s (b=1 and d=1) 
h = 1     # and sh for heterozygous individuals
a = 0     # coefficient for allee effect (growth or death)

difW = 1   # diffusion coefficient for WW individuals
difH = 1   # diffusion coefficientrate for WD individuals
difD = 1   # diffusion coefficient rate for DD individuals

c = 0.7             # homing rate
homing = "zygote"   # "zygote" or "germline"

# Initialization

CI = "center"       # "center" for having the border in the center, "left" for having the border on the left

# How Population Grow and Decrease

growth_dynamic = "logistical"     # exponential or logistical (growth rate is r+1 for expon, r*coef+1 for logist)
death_dynamic = "exponential"      # exponential or logistical  (death rate is always 1)

max_capacity = 1                   # for logistical growth or death


# Linearization
  
linear_growth = False    # With linear_growth = True, max_capacity = 1 is mandatory for a good linear approx. 
linear_mating = False          

# Numerical

T = 100         # final time
L = 400         # length of the spatial domain
M = T*6        # number of time steps
N = L           # number of spatial steps


theta = 0.5     # discretization in space : theta = 0.5 for Crank Nicholson
                # theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  
                   
# Graph

wild =True; heterozygous = True; drive = True              # What to draw on the graph
theorical_wave = False                                     # Draw the theorical wave or not

origin = 190                                               # where to center the theorical wave(made to be changed manually)   
grid = True                                                # A grid or not
semilogy = False                                           # semilogy = False : classical scale, semilogy = True : log scale for y
ylim = False                                               # y scale on the graph (ylim = False, means it's not specify)
xlim = False                                               # x scale on the graph (ylim = False, means it's not specify)
mod = T                                                    # Draw graph every ..mod.. time
graph_type = "Individuals"                                 # "Individuals" or "Proportions" (or None if we don't want any graph)
save_figure = True                                         # Save the figures (.pdf) 

# Speed calculus
speed_proportion = False            # True : use the wild-type number to compute the speed, False : use the wild-type proportion. 



# Group parameters for lisibility
bio_para = [r,s,h,a,difW,difH,difD,c,homing]
model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
num_para = [T,L,M,N,theta]
graph_para = [wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion]



############################ What to do ? #######################################

what_to_do = "evolution" 

# Bring the principal parameters together to make it easier.
# "evolution" : simplest task, draw the propagation regarding the parameters above.
# "system+evolution" : draw the propagation for the linear system and add the theoritical values at the end.
# "calcul de la vitesse en fonction du temps" : idem + draw the speed as a function of time.
# "heatmap" : draw an heatmap



############################### Main program #########################################

if what_to_do == "evolution" :        
    speed, W, H, D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
    
    
    
if what_to_do == "system+evolution" :
    from which_system import which_system 
    from graph import graph
    
    CI = "center"
    graph_type = "Individuals"
    linear_growth = False ; linear_mating = False              
  
    # update parameters
    bio_para = [r,s,h,a,difW,difH,difD,c,homing]
    model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
    num_para = [T,L,M,N,theta]
    graph_para = [wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion]
    
    # run evolution
    speed,position,W,H,D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
          
    # Determine the position of the border and print the number of the system involve
    epsilon = 0.0001              # numerical error accepted
    window = 'minimum'            # where to look : 'zoom border' or 'everywhere'  or 'minimum'   
    border, system, system_bis = which_system(W,H,D,epsilon,window,num_para)     
  
   
    # Draw the wave at the end (normal and log scale), and the theoritical values.
    if speed != None and border!=None:  
        X = np.linspace(0,N,N+1)*(L/N)
        semilogy = False
        xlim = (border-75,border+75)
        ylim = (-0.3,1.3)
        graph(X,W,H,D,T,speed,graph_para,r,s,L)
        semilogy = True 
        xlim = (0,L)                
        ylim = (10**(-60),10)
        graph(X,W,H,D,T,speed,graph_para,r,s,L) 
        
        
     
                
        
        
# Calcul de la vitesse en fonction du temps 
if what_to_do == "vitesse en fonction du temps" :
    from tanaka import tanaka
    
    CI = "center"
    graph_type = None
    linear_growth = False ; linear_mating = False   
    
    # update parameters
    model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
    graph_para = [wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion]
    
    
    speed_proportion = True   # To compare speed, we have to use the WT proportion front (because Tanaka doesn't deal with number)
    
    for s in [0.3,0.8] : 

        speed, W, H, D, speed_fct_of_time_girardin = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
        p_cubic, speed_cubic, speed_fct_of_time_cubic = tanaka(s,"cubic",model_para,num_para,graph_para) 
        p_fraction, speed_fraction, speed_fct_of_time_fraction = tanaka(s,"fraction",model_para,num_para,graph_para)
    
        print('s=',s)
        fig, ax = plt.subplots()  
        ax.plot(speed_fct_of_time_cubic[0,:], speed_fct_of_time_cubic[1,:], label="Tanaka cubic")
        ax.plot(range(0,T+1,mod), np.ones(int(T/mod+1))*(2-3*s)/np.sqrt(2*s), label="Tanaka cubic solution part.")
        if s < 0.5 : 
            ax.plot(range(0,T+1,mod), np.ones(int(T/mod+1))*2*np.sqrt(1-2*s), label="KKP r=1-2*s")
        ax.plot(speed_fct_of_time_fraction[0,:], speed_fct_of_time_fraction[1,:], label="Tanaka fraction") 
        ax.plot(speed_fct_of_time_girardin[0,:],speed_fct_of_time_girardin[1,:], label="Girardin")
        ax.set(title = f"Vitesse en fonction du tps avec s={s}")
        ax.set(xlabel='Time', ylabel='Speed')   
        ax.grid();plt.legend(); plt.show()        
    
        
        if  linear_growth == True and linear_mating == True :    
            # Determine the position of the border
            epsilon = 0.0001             
            window = 'minimum'        
            border, system, system_bis = which_system(W,H,D,epsilon,window) 
             
            # Draw the wave at the end (normal and log scale), and the theoritical values.
            if speed != None :                                                         
                X = np.linspace(0,N,N+1)*(L/N)
                semilogy = False
                xlim = (border-75,border+75)
                ylim = (-0.3,1.3)
                graph(X,W,H,D,T,"Individuals",speed,border)  
                semilogy = True 
                xlim = (0,L)                
                ylim = (10**(-60),10)
                graph(X,W,H,D,T,"Individuals",speed,border)  
            
            
   
if what_to_do == "vitesse en fonction de s" :  
    from tanaka import tanaka
    
    T = 1000; L = 4000; M = T*6; N = L; mod=20            
    CI = "left"
    graph_type = None
    linear_growth = False ; linear_mating = False   
    model = 'comparaison'    
    
    # update parameters
    model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
    num_para = [T,L,M,N,theta]
    graph_para = [wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion]
    
   
    s_min = 0.3 ; s_max = 0.8
    s_values = np.linspace(s_min,s_max,31)
    fct_of_s = np.zeros((7,len(s_values))) 
    fct_of_s[0,:] = s_values      # first line = s values
                
    if model == 'girardin' : 
        for s_index in range(len(s_values)) :
            s = np.round(s_values[s_index],3); print(s)
            
            speed, W, H, D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
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
        ax.plot(fct_of_s[0,:],fct_of_s[1,:], label="vitesse", linestyle='--', color='black',linewidth = 1)
        ax.plot(fct_of_s[0,:],fct_of_s[3,:], label="C", color='blue')
        ax.plot(fct_of_s[0,:],fct_of_s[4,:], label="D-", color='deepskyblue')
        ax.plot(fct_of_s[0,:],fct_of_s[5,:], label="D+", color='deeppink')
        ax.plot(fct_of_s[0,:],np.zeros(len(s_values)), color='black',linewidth = 0.5)
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
        
            
    if model == 'comparaison' :
        
        for s_index in range(len(s_values)) :
            s = np.round(s_values[s_index],3); print(s) 
            
            fct_of_s[1,s_index] = (2-3*s)/np.sqrt(2*s)     # theorical speed for one solution of Tanaka cubic  
            fct_of_s[2,s_index] = tanaka(s,"cubic",model_para,num_para,graph_para)[1]   # numerical speed for Tanaka cubic
            fct_of_s[3,s_index] = tanaka(s,"fraction",model_para,num_para,graph_para)[1]    # numerical speed for Tanaka fraction
            fct_of_s[4,s_index], W, H, D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)  # speed for Leo and Florence's model
            if s < 0.5 : 
                fct_of_s[5,s_index] = 2*np.sqrt(1-2*s)          # speed of KKP model (with r=1-2*s)
    
         
        np.savetxt(f'vitesse_en_fct_de_s_model_{model}.txt', fct_of_s) 
        fig, ax = plt.subplots()    
        #ax.plot(fct_of_s[0,:],fct_of_s[1,:], label="Cubic sol.part.")
        ax.plot(fct_of_s[0,:],fct_of_s[2,:], label="Cubic num.")
        ax.plot(fct_of_s[0,:],fct_of_s[3,:], label="Fraction num.")
        #ax.plot(fct_of_s[0,:],fct_of_s[4,:], label="Leo/Flo num.")
        ax.plot(fct_of_s[0,:],fct_of_s[5,:], label="KKP r=1-2s")
        ax.grid(); ax.legend(); plt.show()
        fig.savefig(f'vitesse_en_fct_de_s_model_{model}.pdf')

            
            
            
            

if what_to_do == "heatmap" :
        from heatmap import heatmap
        from heatmap import print_heatmap
    
        heatmap_type = "classic"    #   "classic"  "speed_cubic" "speed_fraction" "r_one_minus_n_cubic"  "r_one_minus_n_fraction"                                          

        CI = "center"
        graph_type = None
        precision = 30  # number of value on s and r scale (including 0 and 1) for the heatmap
        
        # update parameters
        graph_para = [wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion]
   
        
        if heatmap_type == "classic" :
            for h in [0.9] :
                T = 500; L = 2000; M = T*6; N = L 
                smin = 0.3; smax = 0.9; rmin = 0 ; rmax = 12 
                linear_growth = False ; linear_mating = False  
                
                # update parameters
                heatmap_para = [precision, smin, smax, rmin, rmax]  
                model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
                num_para = [T,L,M,N,theta]
                             
                s_range, r_range, heatmap_values, zero_line = heatmap(heatmap_type, heatmap_para, mod, bio_para, model_para, num_para, graph_para, what_to_do)
                print_heatmap(heatmap_values, zero_line, "simple", heatmap_type, heatmap_para, bio_para, save_figure)
                print_heatmap(heatmap_values, zero_line, "eradication", heatmap_type, heatmap_para, bio_para, save_figure)
                print_heatmap(heatmap_values, zero_line, "collapse", heatmap_type, heatmap_para, bio_para, save_figure)
        else :   
            T = 1000; L = 4000; M = T*40; N = L 
            smin = 0.3; smax = 0.9; rmin = 50 ; rmax = 60 
            linear_growth = False ; linear_mating = False              
            homing = "zygote"       
            c = 1; h = 1
            
            # update parameters
            heatmap_para = [precision, smin, smax, rmin, rmax] 
            bio_para = [r,s,h,a,difW,difH,difD,c,homing]
            model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
            num_para = [T,L,M,N,theta]
            
            s_range, r_range, heatmap_values, zero_line = heatmap(heatmap_type, heatmap_para, mod, bio_para, model_para, num_para, graph_para, what_to_do)                          
            print_heatmap(heatmap_values, zero_line, None, heatmap_type, heatmap_para, bio_para, save_figure)
            

        
               
        
if what_to_do == "limite r infini" :  
    
    T = 1000; L = 4000; M = T*40; N = L; mod=int(T/10)         #T = 150; L = 900; M = 2000; N = 3000  
    CI = "center"
    graph_type = "Individuals"
    linear_growth = False ; linear_mating = False              
    homing = "zygote"        # "zygote" or "germline" or "minimum"
    c = 1; h = 1
    s = 0.5

    # update parameters
    heatmap_para = [precision, smin, smax, rmin, rmax]  
    model_para = [CI,growth_dynamic,death_dynamic,max_capacity,linear_growth,linear_mating]
    num_para = [T,L,M,N,theta]
    graph_para = [wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion]
   
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
        speed_girardin, W, H, D = evolution(bio_para, model_para, num_para, graph_para, what_to_do)
        p = D/(D+H+W); n = D+H+W
        
        # When r tends to infinity : does r(1-n) cv to 0 ? to  s*p*(2-p)/(1-s+s*(1-p)**2) ? 
        cv_r1moinsn_fct_of_r[1,i] = np.max(abs(r*(1-n) - s*p*(2-p)/(1-s+s*(1-p)**2)))
        cv_r1moinsn_fct_of_r[2,i] = np.max(abs(r*(1-n)))
        
        # When r tends to infinity : does the speed of the system Girardin 2021 tends to the speed Tanaka cubic ? Tanaka fraction ?
        cv_speed_fct_of_r[2,i] = tanaka(s,"cubic",model_para,num_para,graph_para)[1]   
        cv_speed_fct_of_r[3,i] = tanaka(s,"fraction",model_para,num_para,graph_para)[1]
        cv_speed_fct_of_r[4,i] = speed_girardin
             
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



