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

from L_speed_eigen_values import num, continu, discrete, L1, L2

    

### Numerical positions and distances  
def dist(time, index_time, wt, drive):
    
    # Percent pic values
    pic_values = [2000000, 2000000] 
    # pic_values = [0.02*max(wt[index_time[-1],:]), 0.2*max(drive[index_time[-1],:])]
    
    # Position of a reference value (pic_values), the last individuals and a density define as nb_drive
    posi_pic = np.zeros((2,len(index_time)))    
    posi_1 = np.zeros((2,len(index_time)))  
    posi_nb = np.zeros((2,len(index_time)))  
    for i in range(len(index_time)) :   
        j = index_time[i]
        # Percent pic (first line : wt, second line : drive)
        posi_pic[0,i] = np.where(wt[j,:]>pic_values[0])[0][0]*dx 
        posi_pic[1,i] = np.where(drive[j,:]>pic_values[1])[0][0]*dx 
        # Distance from the last individual to the pic (first line : wt, second line : drive)
        posi_1[0,i] = np.where(wt[j,:]!=0)[0][0]*dx
        posi_1[1,i] = np.where(drive[j,:]!=0)[0][0]*dx 
        # Distance from the pic to last individual (first line : wt, second line : drive)
        posi_nb[0,i] = np.where(wt[j,:]>nb_drive)[0][0]*dx
        posi_nb[1,i] = np.where(drive[j,:]>nb_drive)[0][0]*dx 

    fig, ax = plt.subplots()
    ax.plot(time[index_time], (posi_pic[0,:]-posi_1[0,:])-np.mean(posi_pic[0,:]-posi_1[0,:]), label = "wild-type")
    ax.plot(time[index_time], (posi_pic[0,:]-posi_nb[1,:])-np.mean(posi_pic[0,:]-posi_nb[1,:]), label = f"drive")
    ax.set(xlabel='Time', ylabel='Distance')
    ax.set_title(f"Space variation of the last wild-type individual \n and the last site with a drive density above {nb_drive}.")
    plt.legend()
    fig.savefig(f"{directory}/distance_s_{s}_nb_{nb_drive}.png", format='png')
    plt.show()
                
    dist_1 = posi_pic - posi_1
    dist_nb = posi_pic - posi_nb
    return(posi_pic, posi_1, posi_nb, dist_1, dist_nb, pic_values)
    
    
  

# Plot the wave of drive and wild-type individuals
def plot_wave(time, index_time, dx, nb_graph, wt, drive, pic_values, directory):  
    for i in index_time[np.linspace(0,len(index_time)-1,nb_graph).astype(int)]:   
        fig, ax = plt.subplots()
        ax.plot(np.arange(nb_sites)*dx, wt[i,:], label = "wt", color = "cornflowerblue")
        ax.plot(np.arange(nb_sites)*dx, drive[i,:], label = "drive", color = "orange")
        ax.hlines(pic_values[0], xmin = 0, xmax = nb_sites*dx, color = "cornflowerblue"); ax.hlines(pic_values[1], xmin = 0, xmax = nb_sites*dx, color = "orange");
        ax.set(xlabel='Space', ylabel='Number of individuals', ylim = [0,1.1*K*dx])
        ax.set_title(f"t = {time[i]}")
        plt.legend()
        fig.savefig(f"{directory}/t_{time[i]}.png", format='png')
        plt.show()
        
        
# Plot the individual densities at the end of the wave (where it decreases exponentially), for t = index_time[-1]
def plot_end_wave(wt, drive, index_time, pic_values, dx):
    for allele in range(2):
        matrix = [wt, drive][allele]
        title = ["wt", "drive"][allele]
        vect = matrix[index_time[-1],:]
        last_index = np.where(vect>pic_values[allele])[0][0]   
        first_individual = np.where(vect>0)[0][0]      
        abscisse = np.arange(first_individual-20-last_index,0)*dx
        fig, ax = plt.subplots()
        ax.semilogy(abscisse, vect[first_individual-20:last_index], label="data")
        ax.set_title(f"dx={dx}, {title}")
        plt.legend()
        fig.savefig(f"{directory}/dx_{dx}_{title}_end_wave.png", format='png')
        plt.show()
    
        
# Plot the distance from the last individual (density=1) to the first density = nb_drive
def eradication_time_mvt(time, index_time, matrix, nb_drive):
    
    #dist_nb = np.zeros((2,len(index_time))) 
    #for i in range(len(index_time)) :   
    #    j = index_time[i]
    #    # Distance from the last individual (density=1), to density = nb_drive
    #    dist_nb[0,i] = (np.where(wt[j,:]>nb_drive)[0][0] - np.where(wt[j,:]!=0)[0][0]) * dx
    #    dist_nb[1,i] = (np.where(drive[j,:]>nb_drive)[0][0] - np.where(drive[j,:]!=0)[0][0]) * dx 
    
    #fig, ax = plt.subplots()
    #bins = [x - 0.5 for x in range(60)]
    #ax.hist([dist_1_nb[0,:]/speed, dist_1_nb[1,:]/speed], bins = bins, histtype = 'bar', label = ['wt','drive'], density=True)
    #ax.vlines(np.mean(dist_1_nb[0,:]/speed), 0, 0.5, color="black", linewidth=2, linestyle="-.")
    #ax.vlines(np.mean(dist_1_nb[1,:]/speed), 0, 0.5, color="black", linewidth=2, linestyle="-.")
    #ax.set(xlabel='Time', xlim = [0,60], ylim = [0,0.35])
    #ax.set_title(f"Distance from density=1 to density={nb_drive} divided by the speed.\n dt = {np.round(dt,10)}, runs = {len(dist_nb[0,:])},  mean wt = {np.round(np.mean(dist_1_nb[0,:]),3)} and drive = {np.round(np.mean(dist_1_nb[1,:]),3)}")
    #plt.legend()
    #fig.savefig(f"{directory}/with_space_extinction_lenght_0_{nb_drive}_divided_by_speed.png", format='png')
    #plt.show() 
    
    #max_density_in_space = np.where(matrix[index_time[0],:]==np.max(matrix[index_time[0],:]))[0][0]

    
    # Find the max to work only on the increasing section (important in the drive case)
    # FIRST TIME : on the right of this point, the density is always > nb_drive (before the pic for drive) 
    space_index_1 = np.where(matrix[index_time[0],:]<nb_drive)[0][-1]+1
    # LAST TIME : on the left of this point, the density = 0 everywhere 
    space_index_2 = np.where(matrix[index_time[-1],:]!=0)[0][0]-1
    # We consider the space in between (space_index_1 -> space_index_2), i.e. spatial point where the density went from < nb_drive to 0.
    # For each spatial site in this area, we record the time until extinction.
    erad_time = np.ones(space_index_2-space_index_1)*(-1) 
        
        
        # Draft
        #for t in np.arange(index_time[0], index_time[0]+15): 
        #    print(t)
        #    vect = wt[t, space_index_1-30:space_index_1+30]
        #    fig, ax = plt.subplots()
        #    bins = [x - 0.5 for x in range(len(vect)+1)]
        #    plt.hist(np.arange(len(vect)), bins = bins, weights = vect,
        #             histtype = 'barstacked', color = 'cornflowerblue')  
        #    ax.set(xlabel=f'Spatial sites', ylabel='Number of individuals')
        #    ax.set_title(f"Density at time {t}")
        #    ax.set(ylim = [0,15])
        #    fig.savefig(f"test_ext_time/{t}.png", format='png')  
        #    plt.show() 
            

    for i in range(space_index_2-space_index_1):
        x = np.arange(space_index_1,space_index_2+1)[i]
        # Find the max to work only on the decreasing section regarding to time (important in the drive case)
        #max_density_in_time = np.where(matrix[:,x]==np.max(matrix[:,x]))[0][0]
        # first time after which density < nb_drive forever (in the time considered)
        
        index_nb = np.where(matrix[:,x]>nb_drive)[0][-1]         # mean = 10.35     
        #index_nb = np.where(matrix[max_density_in_time:,x]<nb_drive)[0][0]+max_density_in_time           # mean = 9.33
                                  
        # first index of time with density == 0 (after the pic)
        #index0 = np.where(matrix[max_density_in_time:,x]==0)[0][0]+max_density_in_time
        # first time after which density == 0 forever (in the time considered)
        index_0 = np.where(matrix[:,x]!=0)[0][-1]
        # Save extinction time
        erad_time[i] = time[index_0] - time[index_nb]
        # Density going from nb_drive to 0, function of time
        #if i == 0 or i == space_index_2-space_index_1-1 or i == (space_index_2-space_index_1-1)//2 :
        #    fig, ax = plt.subplots()
        #    ax.plot(time[index_nb:index_0], matrix[index_nb:index_0,x])
        #    ax.set(xlabel='Time', ylabel='Density')
        #    ax.set_title(f"Density for site {x}, allele {allele} (0=wt,1=drive)")
        #    fig.savefig(f"{directory}/density_for_site_{x}_allele_{allele}.png", format='png')
        #    plt.show()   
    
    
    #fig, ax = plt.subplots()
    #bins = [x - 0.5 for x in range(60)]
    #ax.hist([erad_wt, erad_drive], bins = bins, histtype = 'bar', label = ['wt','drive'], density=True)
    #ax.vlines(np.mean(erad_wt), 0, 0.5, color="black", linewidth=2, linestyle="-.")
    #ax.vlines(np.mean(erad_drive), 0, 0.5, color="black", linewidth=2, linestyle="-.")
    #ax.set(xlabel='Time', xlim = [0,60], ylim = [0,0.35])
    #ax.set_title(f"Extinction time for each spatial site at the end of the wave.\n dt = {np.round(dt,10)}, runs = {len(erad_wt)} and {len(erad_drive)}, mean wt = {np.round(np.mean(erad_wt),3)} and drive = {np.round(np.mean(erad_drive),3)}")
    #plt.legend()
    #fig.savefig(f"{directory}/with_space_extinction_time_for_each_spatial_site.png", format='png')
    #plt.show() 
    
    return(erad_time) #, erad_drive)
            
            
    
# Distance function of time, and histogram
def plot_distances(index_time, dist_1, dist_nb, directory, pic_values, nb_drive): 
    for density in ["1_to_2000000",f"{nb_drive}_to_2000000",f"1_to_{nb_drive}"] :
        # In between two densities
        if density == "1_to_2000000": 
            dist = dist_1; first_dens = 1; last_dens = pic_values[0]
        if density == f"{nb_drive}_to_2000000": 
            dist = dist_nb; first_dens = nb_drive; last_dens = pic_values[0]
        if density == f"1_to_{nb_drive}": 
            dist = dist_1 - dist_nb; first_dens = 1; last_dens = nb_drive
        # distances function of time
        fig, ax = plt.subplots()
        ax.plot(time[index_time], dist[0,:], label = "wt")
        ax.plot(time[index_time], dist[1,:], label = "drive")
        ax.set(xlabel='Time', ylabel='Distance')
        ax.set_title(f"Distance from density {first_dens} to density {last_dens}.")
        plt.legend()
        fig.savefig(f"{directory}/distance_time_{density}.png", format='png')
        plt.show()   
        # histogram of these values
        fig, ax = plt.subplots()
        bins = [x - 0.5 for x in range(int(max(dist[1,:]))+1)]
        ax.hist([dist[0,:], dist[1,:]], bins = bins, histtype = 'bar', label = ['wt','drive'], density=True)
        ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
        ax.set_title(f"Distance from density {first_dens} to density {last_dens}.")
        #plt.legend()
        fig.savefig(f"{directory}/distance_histogram_{density}.png", format='png')
        plt.show()
        
        
def comparaison_distance_time(v, dist_1, dist_nb, erad, gw, gw_int, allele):
    
    dis = (dist_1 - dist_nb)[allele]
    title = ["Wt","Drive"][allele]
    col = [["blue","cornflowerblue", "green", "yellowgreen"],["red", "orange", "gold", "yellowgreen"]][allele]
    # histogram of these value
    fig, ax = plt.subplots()
    bins = [x - 0.5 for x in range(int(np.max(dis))+10)]
    ax.hist([dis, erad*v, gw*v, gw_int*v], bins = bins, histtype = 'bar', label = ['dist','time*speed', 'gw*speed','gw_int*speed'], density=True, color = col)
    for i in range(4) :
        vect = [dis, erad*v, gw*v, gw_int*v]
        ax.vlines(np.mean(vect[i]), 0, 0.3, color=col[i], linewidth=2, linestyle="-.")
    ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
    ax.set_title(f"{title} : distance VS time of eradication * speed (from {nb_drive} -> 0).")
    plt.legend()
    fig.savefig(f"{directory}/distance_time_0_{nb_drive}_{title}.png", format='png')
    plt.show()
    
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(dis)), dis, s=5)
    ax.set(xlabel='Time', ylabel=f'Distance (from {nb_drive} -> 0)')
    ax.set_title(f"{title} : Distance (from {nb_drive} -> 0) in the wave, at each time.") 
    fig.savefig(f"{directory}/{title}_dist.png", format='png')
    plt.show()
    
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(erad)), erad, s=5)
    ax.set(xlabel='Sites', ylabel='Eradication time')
    ax.set_title(f"{title} : Eradication time, at each spatial site.") 
    fig.savefig(f"{directory}/{title}_erad.png", format='png')
    plt.show()
    
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(gw)), gw, s=5)
    ax.set(xlabel='Sites', ylabel='Eradication time')
    ax.set_title(f"{title} : Galton-Watson eradication time.") 
    fig.savefig(f"{directory}/{title}_gw.png", format='png')
    plt.show()
    
            
            
        
def plot_histo_on_wave(wt, drive, index_time, nb_drive, dist_1, dist_nb):     
    dist = dist_1 - dist_nb + 1     
    last_index = np.where(wt[index_time[-1],:]>nb_drive)[0][0]
    first_individual = min([np.where(drive[index_time[-1],:]>10)[0][0], np.where(wt[index_time[-1],:]>0)[0][0]])
    abscisse = np.arange(first_individual-20-last_index,0)*dx
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    bins = [x - 0.5 for x in range(int(abscisse[0]),0)]
    ax1.hist(-dist[0,:], bins = bins, histtype = 'bar', label = ['dist erad'], density=True, color = "lightblue")
    ax2.semilogy(abscisse, wt[index_time[-1], first_individual-20:last_index], label="wt")
    ax2.semilogy(abscisse, drive[index_time[-1], first_individual-20:last_index], label="drive")
    ax2.hlines(nb_drive, (first_individual-20-last_index)*dx, 0, color="black", linestyle="-.") 
    ax1.set_title(f"Extinction histogramme wt (black line = {nb_drive})")
    plt.legend()
    ax1.set_xlabel('Distances')
    ax2.set_ylabel('Densities')
    ax1.set_ylabel('Histogramme proportions')
    fig.savefig(f"{directory}/histogramme_end_wave_{nb_drive}.png", format='png')
    plt.show()
    

    
    

### Parameters and datas ###
    
# Choose K and dx   
dx = 1     
s = 0.4  
m = 0.2 
nb_drive = 100

# Load the other parameters
directory = f"dx_{dx}_s_{s}_m_{m}"
file = open(f"{directory}/parameters.txt", "r")
para = file.read()
para_list = para.replace(' ', '').split("\n")[:-1]
print(para_list)
K = float(para_list[1].replace('K=', ''))
dt = float(para_list[5].replace('dt=', ''))
#m = float(para_list[6].replace('m=', ''))
file.close()

 
# Load datas
time = np.loadtxt(f"{directory}/time.txt")              # Vector time
wt = np.loadtxt(f"{directory}/nW_matrix.txt")           # Matrix wt : columns = times, row = space
drive = np.loadtxt(f"{directory}/nD_matrix.txt")        # Matrix drive : columns = times, row = space
gw = np.loadtxt(f"{directory}/wt_exp_m_0.2_extinction.txt")
gw_int = np.loadtxt(f"{directory}/wt_exp_m_0.2_extinction_tps_entier.txt")


# Determine parameters from data
start = int(3*time[-1]/5)   # Time sequence to look at the exponential section (not to early because of CI)
end = int(time[-1]-10)      # Time sequence to look at the exponential section (not to late because of the section being partly outside the window)
index_time = np.intersect1d(np.where(time>start)[0], np.where(time<end)[0])
    
    
# Parameters for figures
nb_ind = [10, 10]
nb_graph = 2       # Number of graphs shown
nb_sinus = 30      # Number of time values used in the sinus graph
show_graph = False
nb_sites = np.shape(wt)[1]



# Save

### Results and save
directory = f"3_space_dx_{dx}_s_{s}_m_{m}"
if not os.path.exists(directory): os.mkdir(directory)

# Preliminar
posi_pic, posi_1, posi_nb, dist_1, dist_nb, pic_values = dist(time, index_time, wt, drive)
v_num = num(index_time, time, dx, wt, drive, nb_ind, nb_drive, posi_pic, pic_values, dist_1)[0]
   

# Plot the wave and its speed, and distances of the last individuals to 20% of the pic.
plot_wave(time, index_time, dx, nb_graph, wt, drive, pic_values, directory)
#plot_distances(index_time, dist_1, dist_nb, directory, pic_values, nb_drive)
#erad_wt = eradication_time_mvt(time, index_time, wt, nb_drive)
#plot_end_wave(wt, drive, index_time, pic_values, dx)
#plot_histo_on_wave(wt, drive, index_time, nb_drive, dist_1, dist_nb)

# , erad_drive
#comparaison_distance_time(v_num, dist_1, dist_nb, erad_wt, gw, gw_int, 0)


# Plot the renormalized graphs.
#plot_sinus(index_time, dx, lambda_pos_dis, lambda_neg_dis, L_dis1, "Nothing", "Nothing", nb_ind, pic_values, directory)
#plot_sinus(index_time, dx, lambda_pos_dis, lambda_neg_dis, L_dis2, A_dis2, B_dis2, nb_ind, pic_values, directory)


speed_eigen_value = False

if speed_eigen_value : 
    # Speed and eigen values
    v_num, lambda_pos_num, L_num_mean, L_var = num(index_time, time, dx, wt, drive, nb_ind, nb_drive, posi_pic, pic_values, dist_1)
    v_cont, lambda_pos_cont, lambda_neg_cont = continu(s)
    v_dis, lambda_pos_dis, lambda_neg_dis = discrete(s,m,dx,dt)
    L_cont1 = L1(pic_values, lambda_pos_cont) 
    L_dis1 = L1(pic_values, lambda_pos_dis)  
    L_cont2, A_cont2, B_cont2 = L2(pic_values, lambda_pos_cont, lambda_neg_cont, dx) 
    L_dis2, A_dis2, B_dis2 = L2(pic_values, lambda_pos_dis, lambda_neg_dis, dx)    
    
    # Print results
    def p(text, num): print(text, np.round(num,3))
    print("\n- speed"); p("num :", v_num); p("con :", v_cont); p("dis :",v_dis)
    print("\n- lamb pos drive");p("num :",lambda_pos_num[0]); p("con :",lambda_pos_cont[0]); p("dis :",lambda_pos_dis[0])
    print("\n- lamb pos wt");p("num :",lambda_pos_num[1]); p("con :",lambda_pos_cont[1]); p("dis :",lambda_pos_dis[1])
    print("\n- distance wt"); p("num  :",L_num_mean[0]); p("con1 :",L_cont1[0]); p("dis1 :",L_dis1[0]); p("con2 :",L_cont2[0]); p("dis2 :",L_dis2[0])
    print("\n- distance drive"); p("num  :",L_num_mean[1]); p("con1 :",L_cont1[1]); p("dis1 :",L_dis1[1]); p("con2 :",L_cont2[1]); p("dis2 :",L_dis2[1])
    
    # Save the results in the file "figures_values.txt"
    file = open(f"{directory}/dx_{dx}.txt", "w") 
    file.write(f"Speed cont : {v_cont}")
    file.write(f"\nSpeed num  : {v_num}")
    file.write(f"\nSpeed dis : {v_dis}")
    file.write(f"\n\nWT L num : {L_num_mean[0]}")
    file.write(f"\nWT L dis1 : {L_dis1[0]}")
    file.write(f"\nWT L dis2 : {L_dis2[0]}")
    file.write(f"\nWT L con1  : {L_cont1[0]}")
    file.write(f"\nWT L con2  : {L_cont2[0]}")
    file.write(f"\n\nDrive L num : {L_num_mean[1]}")
    file.write(f"\nDrive L dis1 : {L_dis1[1]}")
    file.write(f"\nDrive L dis2 : {L_dis2[1]}")
    file.write(f"\nDrive L con1  : {L_cont1[1]}")
    file.write(f"\nDrive L con2  : {L_cont2[1]}")
    file.close() 

