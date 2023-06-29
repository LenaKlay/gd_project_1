#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 17:15:16 2022

@author: lena
"""

### Librairies

import numpy as np
import matplotlib.pyplot as plt
import os


### Load functions

from L_speed_eigen_values_ger import num, continu, discrete, L
from gw_space import gw    

### Numerical positions and distances  
def dist(time, index_time, wt, drive, K, ref_values, dir_save):
    
    # Position of a reference value (ref_values), the last individuals and a density define as nb_drive
    posi_pic = np.zeros((2,len(index_time)))    
    posi_1 = np.zeros((2,len(index_time)))  
    posi_nb = np.zeros((2,len(index_time)))  
    for i in range(len(index_time)) :   
        j = index_time[i]
        # Percent pic (first line : wt, second line : drive)
        posi_pic[0,i] = np.where(wt[j,:]>ref_values[0])[0][0]*dx 
        posi_pic[1,i] = np.where(drive[j,:]>ref_values[1])[0][0]*dx 
        # Distance from the last individual to the pic (first line : wt, second line : drive)
        posi_1[0,i] = np.where(wt[j,:]!=0)[0][0]*dx
        posi_1[1,i] = np.where(drive[j,:]!=0)[0][0]*dx 
        # Distance from the pic to last individual (first line : wt, second line : drive)
        posi_nb[0,i] = np.where(wt[j,:]>nb_drive)[0][0]*dx
        posi_nb[1,i] = np.where(drive[j,:]>nb_drive)[0][0]*dx 

    #fig, ax = plt.subplots()
    #ax.plot(time[index_time], (posi_pic[0,:]-posi_1[0,:])-np.mean(posi_pic[0,:]-posi_1[0,:]), label = "wild-type")
    #ax.plot(time[index_time], (posi_pic[0,:]-posi_nb[1,:])-np.mean(posi_pic[0,:]-posi_nb[1,:]), label = f"drive")
    #ax.set(xlabel='Time', ylabel='Distance')
    #ax.set_title(f"Space variations of the last wild-type individual \n and the last site with a drive density above {nb_drive}.")
    #plt.legend()
    #fig.savefig(f"{dir_save}/distance_s_{s}_nb_{nb_drive}.png", format='png')
    #plt.show()
                
    dist_1 = posi_pic - posi_1
    dist_nb = posi_pic - posi_nb
    return(posi_pic, posi_1, posi_nb, dist_1, dist_nb)
    
    
  

# Plot the wave of drive and wild-type individuals
def plot_wave(time, index_time, dx, nb_graph, wt, drive, ref_values, dir_save):  
    for i in index_time[np.linspace(0,len(index_time)-1,nb_graph).astype(int)]:   
        fig, ax = plt.subplots()
        ax.semilogy(np.arange(nb_sites)*dx, wt[i,:], label = "wt", color = "cornflowerblue")
        ax.semilogy(np.arange(nb_sites)*dx, drive[i,:], label = "drive", color = "crimson")
        ax.hlines(ref_values[0], xmin = 0, xmax = nb_sites*dx, color = "black", linewidth = 1, linestyle ="--")
        ax.set(xlabel='Space', ylabel='Number of individuals') #, ylim = [0,1.1*K*dx])
        ax.set_title(f"t = {time[i]}")
        plt.legend()
        fig.savefig(f"{dir_save}/t_{time[i]}.png", format='png')
        plt.show()
        
               
    
# to test the function : lambda_back =[lWbd, lDbd]  
    
def plot_end_wave(wt, drive, nb_drive, ref_values, index_time, lambda_back, dx, s, dir_save):  
    # Left margin after the last drive individual (blanc space)
    left_margin = int(np.round(5/dx,0))
    # Number of drive wave superposed in the graphic
    nb_values = 400
    # Vector to calculate the intercept mean, for drive and wt
    density = np.ones((2,len(index_time[-nb_values:]))).astype('int')*(-1)  
    # Loop to shape the figure on the datas
    min_abscisse = np.zeros(nb_values)
    for i in index_time[-nb_values:] :
        last_index = np.where(wt[i,:]>1000)[0][0]  #max(np.where(drive[i,:]> ref_values[1])[0][0], 
        first_individual = min([max([np.where(drive[i,:]>0)[0][0], nb_sites//2]), np.where(wt[i,:]>0)[0][0]])
        min_abscisse[i-index_time[-nb_values]] = (first_individual-left_margin-last_index+1)*dx  
    # Figure             
    fig, ax = plt.subplots(figsize=[-min(min_abscisse)/10, 5.25])
    # Superposition of the numerical waves
    for i in index_time[-nb_values:] :
         #index_wt_ref_value = np.where(wt[i,:]> ref_values[0])[0][0]
         index_drive_ref_value = np.where(drive[i,:]> ref_values[1])[0][0]         
         last_index = max(index_drive_ref_value, np.where(wt[i,:]>1000)[0][0])  # max : to see some wt if there are far away
         first_individual = min([max([np.where(drive[i,:]>0)[0][0], nb_sites//2]), np.where(wt[i,:]>0)[0][0]])
         abscisse = np.arange(first_individual-left_margin-last_index+1,1)*dx
         ax.semilogy(abscisse, drive[i, first_individual-left_margin:last_index], color="crimson", alpha=0.3) 
         ax.semilogy(abscisse, wt[i, first_individual-left_margin:last_index], color="cornflowerblue", alpha=0.3)  
         density[:, i-index_time[-nb_values]] = [wt[i, last_index-1], drive[i, last_index-1]]
    # Legend
    ax.semilogy(0,0, color="crimson", label="drive") 
    ax.semilogy(0,0, color="cornflowerblue", label="wild-type")
    ax.set_ylim([None,ref_values[1]]); ax.set_xlim([min(min_abscisse),0])
    # Theoritical exponential function
    beta = np.mean(density[1,:]); exp_drive = beta*np.exp(lambda_back[1]*abscisse)
    gamma = np.mean(density[0,:]); exp_wt = gamma*np.exp(lambda_back[0]*abscisse)
    exp_D1 = np.where(exp_drive>1)[0][0]; exp_WT1 = np.where(exp_wt>1)[0][0]
    ax.semilogy(abscisse[exp_D1:], exp_drive[exp_D1:], color="black")
    ax.semilogy(abscisse[exp_WT1:], exp_wt[exp_WT1:], color="black")           
    #ax.set_title(f"Wild-type and drive densities at the back of the wave.")
    ax.set_ylabel('Densities'); ax.set_xlabel('Space') 
    plt.legend()
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/end_wave_s_{s}.png", format='png')    
    plt.show()
    
 
# Plot the distance from the last individual (density=1) to the first density = nb_drive
def eradication_time_mvt(time, index_time, matrix, nb_drive):
    
    # Find the max to work only on the increasing section (important in the drive case)
    # FIRST TIME : on the right of this point, the density is always > nb_drive (before the pic for drive) 
    space_index_1 = np.where(matrix[index_time[0],:]<nb_drive)[0][-1]+1
    # LAST TIME : on the left of this point, the density = 0 everywhere 
    space_index_2 = np.where(matrix[index_time[-1],:]!=0)[0][0]-1
    # We consider the space in between (space_index_1 -> space_index_2), i.e. spatial point where the density went from < nb_drive to 0.
    # For each spatial site in this area, we record the time until extinction.
    erad_time = np.ones(space_index_2-space_index_1)*(-1) 
        
    for i in range(space_index_2-space_index_1):
        x = np.arange(space_index_1,space_index_2+1)[i]
        # first index of time for which the density is > nb_drive forever      
        index_nb = np.where(matrix[:,x]>nb_drive)[0][-1] 
        # first index of time for which the density is != 0
        index_0 = np.where(matrix[:,x]!=0)[0][-1]
        # Save extinction time
        erad_time[i] = time[index_0] - time[index_nb]

    return(erad_time) 
            
            
    
# Distance function of time, and histogram
def plot_distances(index_time, dist_1, dist_nb, ref_values, nb_drive, dir_save): 
    for density in ["1_to_2000000",f"{nb_drive}_to_2000000",f"1_to_{nb_drive}"] :
        # In between two densities
        if density == "1_to_2000000": 
            dist = dist_1; first_dens = 1; last_dens = ref_values[0]
        if density == f"{nb_drive}_to_2000000": 
            dist = dist_nb; first_dens = nb_drive; last_dens = ref_values[0]
        if density == f"1_to_{nb_drive}": 
            dist = dist_1 - dist_nb; first_dens = 1; last_dens = nb_drive
        # distances function of time
        fig, ax = plt.subplots()
        ax.plot(time[index_time], dist[0,:], label = "wt")
        ax.plot(time[index_time], dist[1,:], label = "drive")
        ax.set(xlabel='Time', ylabel='Distance')
        ax.set_title(f"Distance from density {first_dens} to density {last_dens}.")
        plt.legend()
        fig.savefig(f"{dir_save}/distance_time_{density}.png", format='png')
        plt.show()   
        # histogram of these values
        fig, ax = plt.subplots()
        bins = [x - 0.5 for x in range(int(max(dist[1,:]))+1)]
        ax.hist([dist[0,:], dist[1,:]], bins = bins, histtype = 'bar', label = ['wt','drive'], density=True)
        ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
        ax.set_title(f"Distance from density {first_dens} to density {last_dens}.")
        #plt.legend()
        fig.savefig(f"{dir_save}/distance_histogram_{density}.png", format='png')
        plt.show()
        
        
def comparaison_distance_time(v, s, dist_1, dist_nb, erad, gw, allele, dir_save):
    
    dis = (dist_1 - dist_nb)[allele]
    title = ["Wt","Drive"][allele]
    col = [["yellowgreen","green","blue"],["red", "orange", "gold"]][allele]
    # histogram of these value
    fig, ax = plt.subplots()
    bins = [x - 0.5 for x in range(int(np.max(dis))+10)]
    ax.hist([dis, erad*v, gw*v], bins = bins, histtype = 'bar', label = ['distance last indi.','(time wave)*speed', '(time gw)*speed'], density=True, color = col)
    for i in range(3) :
        vect = [dis, erad*v, gw*v]
        ax.vlines(np.mean(vect[i]), 0, 0.3, color=col[i], linewidth=2, linestyle="-.")
    ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
    #ax.set_title(f"{title} : distance VS time of eradication * speed (from {nb_drive} -> 0).")
    plt.legend()
    fig.savefig(f"{dir_save}/distance_time_s_{s}.png", format='png')
    plt.show()
    
    col = [["yellowgreen","green"],["red", "orange"]][allele]
    # histogram of these value
    fig, ax = plt.subplots()
    bins = [x - 0.5 for x in range(int(np.max(dis))+10)]
    ax.hist([dis, erad*v], bins = bins, histtype = 'bar', label = ['distance last indi.','(time wave)*speed'], density=True, color = col)
    for i in range(2) :
        vect = [dis, erad*v]
        ax.vlines(np.mean(vect[i]), 0, 0.3, color=col[i], linewidth=2, linestyle="-.")
    ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
    #ax.set_title(f"{title} : distance VS time of eradication * speed (from {nb_drive} -> 0).")
    plt.legend()
    fig.savefig(f"{dir_save}/distance_time_s_{s}_bis.png", format='png')
    plt.show()
    
    
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(dis)), dis, s=5)
    ax.set(xlabel='Time', ylabel=f'Distance (from {nb_drive} -> 0)')
    ax.set_title(f"{title} : Distance (from {nb_drive} -> 0) in the wave, at each time.") 
    fig.savefig(f"{dir_save}/{title}_dist.png", format='png')
    plt.show()
    
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(erad)), erad, s=5)
    ax.set(xlabel='Sites', ylabel='Eradication time')
    ax.set_title(f"{title} : Eradication time, at each spatial site.") 
    fig.savefig(f"{dir_save}/{title}_erad.png", format='png')
    plt.show()
    
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(gw)), gw, s=5)
    ax.set(xlabel='Sites', ylabel='Eradication time')
    ax.set_title(f"{title} : Galton-Watson eradication time.") 
    fig.savefig(f"{dir_save}/{title}_gw.png", format='png')
    plt.show()
    
            
            
        
def plot_histo_on_wave(wt, drive, index_time, nb_drive, dist_1, dist_nb, s, dir_save):   
    # Left margin on the left of the graphic (blanc space)
    left_margin = int(np.round(20/dx,0))
    # Number of drive wave superposed in the graphic
    nb_values = 400
    # Distance between sites of density 1 and nb_drive
    dist = dist_1 - dist_nb + 1    
    # Loop to shape the figure on the datas
    min_abscisse = np.zeros(nb_values)
    for i in index_time[-nb_values:] :
        last_index = np.where(wt[i,:]>nb_drive)[0][0]
        first_individual = min([max([np.where(drive[i,:]>0)[0][0], nb_sites//2]), np.where(wt[i,:]>0)[0][0]])
        min_abscisse[i-index_time[-nb_values]] = (first_individual-left_margin-last_index+1)*dx  
    # Figure             
    fig, ax1 = plt.subplots(figsize=[-min(min_abscisse)/10, 5.25])
    ax2 = ax1.twinx()
    # Drive   
    for i in index_time[-nb_values:]:
        last_index = np.where(wt[i,:]>nb_drive)[0][0]
        first_individual = min([max([np.where(drive[i,:]>0)[0][0], nb_sites//2]), np.where(wt[i,:]>0)[0][0]]) # min([np.where(drive[i,:]>10)[0][0], np.where(wt[i,:]>0)[0][0]])
        abscisse = np.arange(first_individual-left_margin-last_index+1,1)*dx       
        ax2.semilogy(abscisse, drive[i, first_individual-left_margin:last_index], color="crimson", alpha = 0.3)   
    ax2.semilogy(0, 0, label="Drive", color="crimson", linewidth=2) 
    # Histogramme
    bins = [x - 0.5 for x in range(int(abscisse[0]),0)]
    ax1.hist(-dist[0,:], bins = bins, histtype = 'bar', label = ['dist erad'], density = True, color="yellowgreen", alpha = 0.4)
    # Wild-type
    ax2.semilogy(abscisse, wt[index_time[-1], first_individual-left_margin:last_index], label="Wild-type", color="cornflowerblue", linewidth=2)
    # Wild-type last individual
    index_first_ind = np.where(wt[index_time[-10], first_individual-left_margin:last_index]>0)[0][0]
    ax2.semilogy(abscisse[index_first_ind-1:index_first_ind+1], wt[index_time[-10], first_individual-left_margin+index_first_ind-1:first_individual-left_margin+index_first_ind+1], color="limegreen", linewidth=2, alpha = 0.8, label = "Last wild-type indi.")
    #ax2.hlines(nb_drive, (first_individual-30-last_index)*dx, 0, color="black", linestyle="-.") 
    #ax1.set_title(f"Extinction histogramme wt")
    ax2.set_ylim([None,ref_values[1]]); ax2.set_xlim([min(min_abscisse),0])
    plt.legend()
    ax1.set_xlabel('Distances')
    ax2.set_ylabel('Densities')
    ax1.set_ylabel('Histogramme proportions')
    ax1.set_ylim([0,0.3])
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/histogram_s_{s}.png", format='png')
    plt.show()
    

#def heatmap_distance(v):
#    ref_values = 
#    difference = 
#    L(ref_values, lambda_pos)
#    fig, ax = plt.subplots()
#    im = ax.imshow(harvest)
#    fig.tight_layout()
#    plt.show()
            
    

### Parameters and datas ###
    
# Choose some parameters 
nb_drive = 100              # Threshold for drive ("enough" drive for the chasing not to happen)
conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
dx = 1                      # spatial step
K = 10**8                    # Carrying capacity on one space unit
s = 0.3                      # Disadvantage for drive
r = 0.1                     # Intrasic growth rate


# Load the other parameters
dir_load = f"{conv_timing}_K_{int(np.log10(K))}_dx_{dx}_s_{s}_r_{r}"
file = open(f"{dir_load}/parameters.txt", "r")
para = file.read()
para_list = para.replace(' ', '').split("\n")[:-1]
print(para_list)
K = int(para_list[1].replace('K=', ''))
nb_sites = int(para_list[2].replace('nb_sites=', ''))
T = float(para_list[4].replace('T=', ''))
dt = float(para_list[5].replace('dt=', ''))
m = float(para_list[6].replace('m=', ''))
c = float(para_list[10].replace('c=', ''))
h = float(para_list[11].replace('h=', ''))
file.close()


# Load datas
time = np.loadtxt(f"{dir_load}/time.txt")              # Vector time
wt = np.loadtxt(f"{dir_load}/nW_matrix.txt")           # Matrix wt : columns = times, row = space
drive = np.loadtxt(f"{dir_load}/nD_matrix.txt")        # Matrix drive : columns = times, row = space


# Determine parameters from data
start = int(3*time[-1]/5)   # Time sequence to look at the exponential section (not to early because of CI)
end = int(time[-1]-10)      # Time sequence to look at the exponential section (not to late because of the section being partly outside the window)
index_time = np.intersect1d(np.where(time>start)[0], np.where(time<end)[0])
    
    
# Parameters for figures
ref_values = [0.001*K*dx]*2   # Reference position
nb_ind = [10, 10]
nb_graph = 2       # Number of graphs shown
nb_sinus = 30      # Number of time values used in the sinus graph
show_graph = False



# Save

### Results and save
dir_save = f"{dir_load}/outputs"
if not os.path.exists(dir_save): os.mkdir(dir_save)

# Preliminar
posi_pic, posi_1, posi_nb, dist_1, dist_nb = dist(time, index_time, wt, drive, K, ref_values, dir_save)
#plot_wave(time, index_time, dx, nb_graph, wt, drive, ref_values, dir_save)
v_num, [lWbn, lDbn], L_num = num(index_time, time, dx, wt, drive, nb_ind, nb_drive, posi_pic, ref_values, dist_1)
v_cont, lDfc, [lWbc, lDbc] = continu(conv_timing,s,h,c,r); L_cont = L(ref_values, [lWbc, lDbc]) 
v_dis, lDfd, [lWbd, lDbd] = discrete(conv_timing,s,h,c,r,m,dx,dt); L_dis = L(ref_values, [lWbd, lDbd]) 
   
print(lDfc, lDfd)


# Plot the end of the wave (exp section )
plot_end_wave(wt, drive, nb_drive, ref_values, index_time, [lWbd, lDbd], dx, s, dir_save)

# Eradication time in the wave
erad_wt = eradication_time_mvt(time, index_time, wt, nb_drive)

# Eradication time galton-watson
nb_i = 1000; max_density = 10**5; snapshot = False; left_margin=30
gw_wt = gw(time, start, end, wt, 0, max_density, left_margin, nb_drive, m, nb_i, dx, T, dt, r, s, c, h, lWbd, snapshot, dir_load, dir_save)
  
#Comparaison
comparaison_distance_time(v_num, s, dist_1, dist_nb, erad_wt, gw_wt, 0, dir_save)

#plot_distances(index_time, dist_1, dist_nb, directory, ref_values, nb_drive)
plot_histo_on_wave(wt, drive, index_time, nb_drive, dist_1, dist_nb, s, dir_save)


# Save the results in the file "figures_values.txt"
file = open(f"{dir_load}/speed_lambdas_L.txt", "w") 
file.write(f"Speed num : {v_num}")
file.write(f"\nSpeed cont : {v_cont}")
file.write(f"\nSpeed dis : {v_dis}")
file.write(f"\nLambdas front cont : {lDfc}")
file.write(f"\nLambdas front dis : {lDfd}")
file.write(f"\nLambdas back num : {[lWbn,lDbn]}")
file.write(f"\nLambdas back cont : {[lWbc,lDbc]}")
file.write(f"\nLambdas back dis : {[lWbd,lDbd]}")
file.write(f"\nL num : {L_num}")
file.write(f"\nL cont : {L_cont}")
file.write(f"\nL dis : {L_dis}")
file.close() 




