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




############################ What to do ? #######################################

what_to_do = "evolution"     

# Bring the principal parameters together to make it easier.
# "evolution" : simplest task, draw the propagation regarding the parameters above.
# "system+evolution" : draw the propagation for the linear system and add the theoritical values at the end.
# "calcul de la vitesse en fonction du temps" : idem + draw the speed as a function of time.
# "heatmap" : draw an heatmap


############################## Libraries ########################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors 
import scipy.sparse as sp
import scipy.sparse.linalg  as la
import time

plt.rcParams.update({'font.family':'serif'})
print(plt.rcParams['font.family'])



############################### Parameters ######################################

# Biological
    
r = 0     # growth rate
s = 0.45  # when selection only acts on survival with a fitness cost s (b=1 and d=1) 
h = 1     # and sh for heterozygous individuals
a = 0     # coefficient for allee effect (growth or death)

difW = 1   # diffusion coefficient for WW individuals
difH = 1   # diffusion coefficientrate for WD individuals
difD = 1   # diffusion coefficient rate for DD individuals

c = 1               # homing rate
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
M = T*40        # number of time steps
N = L           # number of spatial steps


theta = 0.5     # discretization in space : theta = 0.5 for Crank Nicholson
                # theta = 0 for Euler Explicit, theta = 1 for Euler Implicit  
                   
# Graph

wild =True; heterozygous = True; drive = True; N_alpha = False        # What to draw on the graph
theorical_wave = False                                     # Draw the theorical wave or not

origin = 190                                               # where to center the theorical wave(made to be changed manually)   
grid = True                                                # A grid or not
semilogy = False                                           # semilogy = False : classical scale, semilogy = True : log scale for y
ylim = False                                               # y scale on the graph (ylim = False, means it's not specify)
xlim = False                                               # x scale on the graph (ylim = False, means it's not specify)
mod = 20                                                   # Draw graph every ..mod.. time
graph_type = "Proportions"                                 # "Individuals" or "Proportions" (or None if we don't want any graph)

# Speed calculus
speed_proportion = False            # True : use the wild-type number to compute the speed, False : use the wild-type proportion. 

############################### Results #########################################



    
    
    
if what_to_do == "evolution" :
    
    from functions.evolution import evolution 
        
    speed, W, H, D = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)
           
    
if what_to_do == "system+evolution" :
    
    from functions.evolution import evolution 
    from functions.which_system import which_system 

    # Penser à changer la configuration graphique !!! (inline vs automatic)
    T = 2000; L = 1000; M = T*6; N = L; mod=250            #T = 150; L = 900; M = 2000; N = 3000  
    CI = "center"
    graph_type = "Individuals"
    linear_growth = False ; linear_mating = False              
    homing = "germline"        # "zygote" or "germline" or "minimum"
    c = 0.7; h = 0.7
    s = 0.7; r = 7
    theta = 0.5        
            

    # Penser à changer la configuration graphique !!! (inline vs automatic)
    #T = 40; L =120; M = T*6; N = L; mod=20             #T = 150; L = 900; M = 2000; N = 3000  
    #CI = "center"
    #graph_type = "Individuals"
    #grid=False
    #linear_growth = False ; linear_mating = False             
    #homing = "zygote"        # "zygote" or "germline" or "minimum"
    #r = 2; s = 0.4
    #theta = 1                
            
    speed,position,W,H,D = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)
    #print('\nVitesse numérique :', speed)
    #print('Vitesse théorique :', np.sqrt((8*(r-0.5)**2)/(3*(1+r))), '\n')
    
    #speed_vector = np.zeros(6)
    #for s_index in range(6) :
    #    s = np.linspace(0.47,0.52,6)[s_index]
    #    speed, position, W, H, D = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)
    #    speed_vector[s_index] = speed   
        
    # Determine the position of the border and print the number of the system involve
    epsilon = 0.0001              # numerical error accepted
    window = 'minimum'            # where to look : 'zoom border' or 'everywhere'  or 'minimum'   
    border, system, system_bis = which_system(W,H,D,epsilon,window)     
    
    #derivative = True
    #X = np.linspace(0,N,N+1)*(L/N) 
    #graph(X,W,H,D,T,"Individuals",speed,border)  
    #derivative = False
    
    
   
    # Draw the wave at the end (normal and log scale), and the theoritical values.
    if speed != None and border!=None:  
        X = np.linspace(0,N,N+1)*(L/N)
        semilogy = False
        xlim = (border-75,border+75)
        ylim = (-0.3,1.3)
        graph(X,W,H,D,T,"Individuals",speed,border)  
        semilogy = True 
        xlim = (0,L)                
        ylim = (10**(-60),10)
        graph(X,W,H,D,T,"Individuals",speed,border)  
        
        
     
                
        
        
# Calcul de la vitesse en fonction du temps 
if what_to_do == "vitesse en fonction du temps" :
    
    # Penser à changer la configuration graphique !!! (inline vs automatic)
    T = 1000; L = 4000; M = T*6; N = L; mod=100         #T = 150; L = 900; M = 2000; N = 3000  
    CI = "center"
    graph_type = None
    linear_growth = False ; linear_mating = False              
    homing = "zygote"        # "zygote" or "germline" or "minimum"
    c = 1; h = 1
    r = 0
    
    speed_proportion = True   # To compare speed, we have to use the WT proportion front (because Tanaka doesn't deal with number)
    
    for s in [0.3,0.8] : 

        speed, W, H, D, speed_fct_of_time_girardin = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)  
        p_cubic, speed_cubic, speed_fct_of_time_cubic = tanaka(s,T,L,M,N,theta,mod,"cubic")
        p_fraction, speed_fraction, speed_fct_of_time_fraction = tanaka(s,T,L,M,N,theta,mod,"fraction")
    
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
    
    # Penser à changer la configuration graphique !!! (inline vs automatic)
    T = 1000; L = 4000; M = T*6; N = L; mod=20            
    CI = "left"
    graph_type = None
    linear_growth = False ; linear_mating = False                
    homing = "zygote"       # "zygote" or "germline" or "minimum"
    r = 0
    theta = 0.5    
    model = 'comparaison'    
    
    
    s_min = 0.3 ; s_max = 0.8
    s_values = np.linspace(s_min,s_max,31)
    fct_of_s = np.zeros((7,len(s_values))) 
    fct_of_s[0,:] = s_values      # first line = s values
                
    if model == 'girardin' : 
        for s_index in range(len(s_values)) :
            s = np.round(s_values[s_index],3); print(s)
            
            speed, W, H, D = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)
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
            fct_of_s[2,s_index] = tanaka(s,T,L,M,N,theta,mod,"cubic")[1]   # numerical speed for Tanaka cubic
            fct_of_s[3,s_index] = tanaka(s,T,L,M,N,theta,mod,"fraction")[1]    # numerical speed for Tanaka fraction
            fct_of_s[4,s_index], W, H, D = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)  # speed for Leo and Florence's model
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
    
        heatmap_type = "classic"    #   "classic"  "speed_cubic" "speed_fraction" "r_one_minus_n_cubic"  "r_one_minus_n_fraction"                                          

        CI = "center"
        graph_type = None
        theta = 0.5                
        precision = 30  # number of value on s and r scale (including 0 and 1) for the heatmap
              
        tic = time.clock()
        
        if heatmap_type == "classic" :
            for h in [0.9] :
                T = 500; L = 2000; M = T*6; N = L 
                smin = 0.3; smax = 0.9; rmin = 0 ; rmax = 12  
                linear_growth = False ; linear_mating = False              
                homing = "zygote"        # "zygote" or "germline" 
                c = 0.9 #; h = 0.7
                s_range, r_range, heatmap_values, zero_line = heatmap(precision,smin,smax,rmin,rmax, "classic")
                print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, "simple")  
                print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, "eradication")  
                print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, "collapse") 
                np.savetxt(f'zero_line_{homing}_c_{c}_h_{h}.txt', zero_line) 
                np.savetxt(f'heatmap_{homing}_c_{c}_h_{h}.txt', heatmap_values)               # heatmap_values = np.loadtxt('heatmap_{homing}_c_{c}_h_{h}.txt') to get it back
        else :   
            T = 1000; L = 4000; M = T*40; N = L 
            smin = 0.3; smax = 0.9; rmin = 50 ; rmax = 60 
            linear_growth = False ; linear_mating = False              
            homing = "zygote"        # "zygote" or "germline" 
            c = 1; h = 1
            s_range, r_range, heatmap_values, zero_line = heatmap(precision,smin,smax,rmin,rmax,heatmap_type)                                  
            print_heatmap(homing, c, h, s_range, r_range, heatmap_values, zero_line, None)  
            np.savetxt(f'heatmap_{heatmap}.txt', heatmap_values)  

        toc = time.clock()
        print("Temps de calcul :", toc - tic)
        
    
               
        
if what_to_do == "limite r infini" :  
    
    T = 1000; L = 4000; M = T*40; N = L; mod=int(T/10)         #T = 150; L = 900; M = 2000; N = 3000  
    CI = "center"
    graph_type = "Individuals"
    linear_growth = False ; linear_mating = False              
    homing = "zygote"        # "zygote" or "germline" or "minimum"
    c = 1; h = 1
    s = 0.5
    theta = 0.5    
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
        speed_girardin, W, H, D = evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod)
        p = D/(D+H+W); n = D+H+W
        
        # When r tends to infinity : does r(1-n) cv to 0 ? to  s*p*(2-p)/(1-s+s*(1-p)**2) ? 
        cv_r1moinsn_fct_of_r[1,i] = np.max(abs(r*(1-n) - s*p*(2-p)/(1-s+s*(1-p)**2)))
        cv_r1moinsn_fct_of_r[2,i] = np.max(abs(r*(1-n)))
        
        # When r tends to infinity : does the speed of the system Girardin 2021 tends to the speed Tanaka cubic ? Tanaka fraction ?
        cv_speed_fct_of_r[2,i] = tanaka(s,T,L,M,N,theta,mod,"cubic")[1]
        cv_speed_fct_of_r[3,i] = tanaka(s,T,L,M,N,theta,mod,"fraction")[1]
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
    

#s=0.8    # Attention à modifier s en haut aussi !!!!
#cv_r1moinsn_fct_of_r = np.loadtxt(f'cv_r_un_moins_n_s_{s}.txt')
#cv_speed_fct_of_r = np.loadtxt(f'cv_speed_s_{s}.txt')

      
# Brouillon          
#      

#r = 0.7; s = 0.51
#alpha = (-(1-2*s)*(1-r) + np.sqrt((1-2*s)**2 * (1-r)**2 + 8 * r**2 * (1-s)))/(4 * r * (1-s))
#A_scr = (2*(1-s)*r)/(2*(1-s)*(1-r+2*alpha*r)-(1-r))              

# dt = T/M    # time
# dx = L/N    # spatial
# X = np.linspace(0,N,N+1)*dx  
        
 
#linear_speed_range = speed_range       
#precision = 20  
#smin = 0.3; smax = 0.9; rmin = 0.1 ; rmax = 5   
#s_range = np.linspace(smin,smax,precision)      
#r_range = np.linspace(rmin,rmax,precision) 
#linear_speed_range = np.loadtxt('heatmap_zygote_c_1_h_1.txt')
#i=18;j=18
#print('r_range =', r_range[i]); print('s_range =', s_range[j])
#print('vitesse', linear_speed_range[i,j])




