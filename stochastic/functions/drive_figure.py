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

### Graphs parameters 

plt.rcParams.update({'font.family':'serif'})
label_size = 12
legend_size = 10

### Load functions

from L_speed_eigen_values import num, continu, discrete, L
from gw_space import gw, ini_exp

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
    # Difference of the positions
    dist_1 = posi_pic - posi_1
    dist_nb = posi_pic - posi_nb
    return(posi_pic, posi_1, posi_nb, dist_1, dist_nb)
    
    


# to test the function : lambda_back =[lWbd, lDbd]  
    
def plot_end_wave(wt, drive, nb_drive, ref_values, index_time, lambda_back, dx, s, x_graph_values, dir_save):  
    # Number of drive wave superposed in the graphic
    nb_values = 2000
    # Vector to calculate the intercept mean, for drive and wt
    ref_inside_window_drive = False; absc_drive_for_mean = np.ones(len(index_time[-nb_values:])).astype('int')*(-1)  
    density_for_mean = np.ones((2,len(index_time[-nb_values:]))).astype('int')*(-1)  
    # Figure    
    abscisse = np.arange(-x_graph_values,1)*dx         
    fig, ax = plt.subplots(figsize=[-int(abscisse[0])/10, 5.25])
    # Superposition of the numerical waves
    for i in index_time[-nb_values:] :
        after_chasing = np.where(wt[i,:]==0)[0][-1]
        last_index = after_chasing + np.where(wt[i,after_chasing:]>ref_values[0])[0][0]  
        first_index = last_index - int(x_graph_values) 
        # Wild-type
        ax.semilogy(abscisse, wt[i, first_index:last_index+1], color="cornflowerblue", alpha=0.3) 
        # Drive
        ax.semilogy(abscisse, drive[i, first_index:last_index+1], color="crimson", alpha=0.3)        
        # Save values for mean
        
        # Wild-type
        value_wt = wt[i, last_index]
        
        # Drive
        if ref_values[1] < drive[i, last_index] :            
            index_drive = after_chasing + np.where(drive[i, after_chasing:] > ref_values[1])[0][0] - first_index
            value_drive = drive[i, index_drive+first_index]; ref_inside_window_drive = True
            
        else: value_drive = drive[i, last_index]
        
        # Save
        
        density_for_mean[:, i-index_time[-nb_values]] = [value_wt, value_drive]  
        if ref_inside_window_drive : absc_drive_for_mean[i-index_time[-nb_values]] = abscisse[index_drive]
        
    # Legend
    ax.semilogy(0,0, color="crimson", label="drive") 
    ax.semilogy(0,0, color="cornflowerblue", label="wild-type")
    ax.set_ylim([None,10**5]); ax.set_xlim([int(abscisse[0]),0])
    # Theoritical exponential function
    
    # Wild-type
    gamma = np.mean(density_for_mean[0,:])
    exp_wt = gamma*np.exp(lambda_back[0]*abscisse)
    exp_WT1 = np.where(exp_wt>1)[0][0]
    ax.semilogy(abscisse[exp_WT1:], exp_wt[exp_WT1:], color="black")   
    
    # Drive
    
    if not ref_inside_window_drive: beta = np.mean(density_for_mean[1,:])
    else : beta = np.mean(density_for_mean[1,:])/np.exp(np.mean(absc_drive_for_mean)*lambda_back[1])
    
    exp_drive = beta*np.exp(lambda_back[1]*abscisse)
    exp_D1 = np.where(exp_drive>1)[0][0]
    
    if ref_inside_window_drive: expDref = np.where(exp_drive>ref_values[1])[0][0]
    else : expDref = -1   
    
    ax.semilogy(abscisse[exp_D1:expDref+1], exp_drive[exp_D1:expDref+1], color="black")
   
    # Graph labels, legend, save
    ax.set_ylabel('Number of individuals'); ax.set_xlabel('Space') 
    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)   
    plt.rc('legend', fontsize=legend_size)      
    plt.legend()
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/end_wave_s_{s}_{int(np.log10(K))}.png", format='png')    
    plt.show()
    
    return(np.where(wt[-10,:]>0)[0][0]<np.where(drive[-10,:]>0)[0][0])

    
# to test the function : lambda_back =[lWbd, lDbd]  
    
def plot_end_wave_wt(wt, index_time, lambda_back, dx, s, x_graph_values, limit_exp, dir_save):  
    # Number of wt wave superposed in the graphic
    nb_values = 2000
    # Vector to calculate the intercept mean, for drive and wt
    absc_wt_for_mean = np.ones(len(index_time[-nb_values:])).astype('int')*(-1)  
    density_wt_for_mean = np.ones(len(index_time[-nb_values:])).astype('int')*(-1)  
    # Figure    
    abscisse = np.arange(-x_graph_values,1)*dx         
    fig, ax = plt.subplots()
    # Superposition of the numerical waves
    for i in index_time[-nb_values:] :
        last_index = np.where(wt[i,:]>=10**7)[0][0]  
        first_index = last_index - int(x_graph_values) 
        # Plot
        ax.semilogy(abscisse, wt[i, first_index:last_index+1], color="cornflowerblue", alpha=0.3)         
        # Wild-type
        index_wt = np.where(wt[i,:] > limit_exp)[0][0] - first_index
        value_wt = wt[i, index_wt+first_index]
        # Save        
        density_wt_for_mean[i-index_time[-nb_values]] = value_wt
        absc_wt_for_mean[i-index_time[-nb_values]] = abscisse[index_wt]
        
    # Legend
    ax.semilogy(0,0, color="cornflowerblue", label="wild-type")
    ax.set_ylim([None,10**7]); ax.set_xlim([int(abscisse[0]),0])
    # Theoritical exponential function
    
    # Exponential function
    gamma = np.mean(density_wt_for_mean)/np.exp(np.mean(absc_wt_for_mean)*lambda_back[0])    
    exp_wt = gamma*np.exp(lambda_back[0]*abscisse)
    exp_W1 = np.where(exp_wt>1)[0][0]
    expWref = np.where(exp_wt>limit_exp)[0][0]    
    ax.semilogy(abscisse[exp_W1:expWref+1], exp_wt[exp_W1:expWref+1], color="black")
   
    # Graph labels, legend, save
    ax.set_ylabel('Number of individuals'); ax.set_xlabel('Space') 
    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)   
    plt.rc('legend', fontsize=legend_size)      
    plt.legend()
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/end_wave_wt_s_{s}.png", format='png')    
    plt.show() 
    
    return(np.mean(absc_wt_for_mean))
    
    
    
    

# Plot the distance from the last individual (density=1) to the first density = nb_drive
def eradication_time_mvt(time, index_time, matrix, nb_drive):
    
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
            
 
        
        
# Galton-Watson : Histogram for different values of max densities  
    
def data_histogram_gw(nb_sites, lWbd, nb_drive, allele, nb_i, T, dt, dx, r, s, c, h, m, x_graph_values, decal, dir_load, dir_save):
    exposants = np.arange(3, int(np.log10(K))-1)
    extinction_list = np.zeros((len(exposants), nb_i))
    for j in exposants:
        max_density = 10**j
        exp, exp_nb = ini_exp(nb_sites, lWbd, nb_drive, max_density, dx, x_graph_values)         
        extinction_list[j-3], exp = gw(exp, exp_nb, allele, nb_i, T, dt, dx, r, s, c, h, m, dir_load, dir_save)
    # adjust x position
    exp_decal = np.zeros(len(exp))    
    exp_decal[np.where(exp!=0)[0]+int(np.round(decal/dx,0))] = exp[np.where(exp!=0)[0]]
    return(exposants, extinction_list, exp_decal)
    
def histogram_gw(exposants, extinction_list, exp, x_graph_values):
    col = ["darkturquoise", "deepskyblue", "dodgerblue", "royalblue", "blue"][5-len(exposants):]
    al = [0.3,0.4,0.5,0.5,0.8]
    abscisse = np.arange(-x_graph_values,1)
    # Initial conditions
    fig, ax = plt.subplots()
    for j in exposants:
        part_exp = np.zeros(len(exp))
        # Densities
        if j == exposants[0]:
            part_exp[:np.where(exp>=10**(j))[0][0]] = exp[:np.where(exp>=10**(j))[0][0]] 
        else : 
            part_exp[np.where(exp>=10**(j-1))[0][0]:np.where(exp>=10**(j))[0][0]] = exp[np.where(exp>=10**(j-1))[0][0]:np.where(exp>=10**(j))[0][0]]      
        bins = [x for x in range(len(exp))]
        kwargs = dict(histtype='stepfilled', alpha=al[j-3], bins=bins, weights = part_exp, color=col[j-3], log=True)
        plt.hist(np.arange(len(part_exp)), **kwargs)  
        # Limit on the right   
        ax.vlines(np.where(exp>=10**j)[0][0], 0, 10**j, color=col[j-3], linewidth=2, label = np.arange(3, int(np.log10(K)))[j-3])
        ax.hlines(10**j, 0, np.where(exp>=10**j)[0][0], color=col[j-3], linewidth=2, label = np.arange(3, int(np.log10(K)))[j-3])
    # Site initially of density close to nb_drive
    ax.scatter(np.where(exp>=nb_drive)[0][0], exp[np.where(exp>=nb_drive)[0][0]], alpha=0.5, color="black", s=100, zorder=10)
    ax.vlines(np.where(exp>=nb_drive)[0][0], 0, exp[np.where(exp>=nb_drive)[0][0]], alpha=0.5, color="black", linewidth=4)        
    ax.set(xlabel='Space', ylabel='Number of individuals')
    ax.set_xticks(np.linspace((len(exp)-1)%10, len(exp)-1, len(np.arange(abscisse[0], abscisse[-1]+1, 10))))    
    ax.set_xticklabels(np.arange(abscisse[(len(exp)-1)%10], abscisse[-1]+1, 10)) 
    ax.set_ylim([None,10**7]); ax.set_xlim([0,-int(abscisse[0])])
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/gw_ini_s_{s}.png", format='png')
    plt.show()
    
    #  Histogram of extinction time values
    al = [0.3,0.4,0.5,0.5,0.8]
    fig, ax = plt.subplots()
    for j in exposants:
        bins = [x - 0.5 for x in range(int(max(extinction_list[j-3,:]))+11)]
        kwargs = dict(histtype='stepfilled', alpha=al[j-3], bins=bins, color=col[j-3], label = f"max density = $10^{j}$", density=True)
        ax.hist(extinction_list[j-3,:]*v_num, **kwargs)
        ax.vlines(np.mean(extinction_list[j-3,:]*v_num), 0, 0.4, color=col[j-3], linewidth=2, linestyle="-.")
    ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
    #ax.set_title(f"{title} : distance VS time of eradication * speed (from {nb_drive} -> 0).")
    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)   
    plt.rc('legend', fontsize=legend_size)  
    plt.legend()
    ax.set(xlabel='Time', xlim = [0,30], ylim = [0,0.35])
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/gw_histogram_s_{s}.png", format='png')
    plt.show()
    
       
     
        
def comparaison_distance_time_gw(v, s, dist_1, dist_nb, erad, gw, allele, dir_save):
 
    # Histogram
    dis = (dist_1 - dist_nb)[allele]
    title = ["time_gw","","distance_time"]
    col = [["blue","green","yellowgreen"],["red", "orange", "gold"]][allele]
    leg = ['(time gw)*speed', '(time wave)*speed','distance last indi.']
    vect = [gw*v, erad*v, dis]
    for ir in [[1,0], [1,2]] : 
        # histogram of these value
        fig, ax = plt.subplots()
        bins = [x - 0.5 for x in range(int(np.max(dis))+10)]
        for i in ir :      
            ax.hist(vect[i], bins = bins, histtype = 'stepfilled', density=True, color = col[i], alpha = (i==ir[0])*0.8+(i==ir[-1])*0.5, linewidth = 2, label = leg[i])        
        for i in ir :     
            ax.hist(vect[i], bins = bins, histtype = 'step', density=True, color = col[i], alpha = 1, linewidth = 2)
            ax.vlines(np.mean(vect[i]), 0, 0.3, color=col[i], linewidth=2, linestyle="-.")
        ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
        #ax.set_title(f"{title} : distance VS time of eradication * speed (from {nb_drive} -> 0).")
        ax.xaxis.label.set_size(label_size)
        ax.yaxis.label.set_size(label_size)   
        plt.rc('legend', fontsize=legend_size)  
        plt.legend()
        plt.tight_layout() 
        fig.savefig(f"{dir_save}/{title[i]}_s_{s}.png", format='png')
        plt.show()
    
    # Scatter plot for the position of the last individual
    fig, ax = plt.subplots()
    plt.scatter(np.arange(len(dis)), dis, s=5)
    ax.set(xlabel='Time', ylabel=f'Distance (from {nb_drive} -> 0)')
    ax.set_title(f"{title} : Distance (from {nb_drive} -> 0) in the wave, at each time.") 
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/{title}_dist.png", format='png')
    plt.show()   
            
        
def plot_histo_on_wave(wt, drive, index_time, nb_drive, dist_1, dist_nb, s, x_graph_values, dir_save):   
    # Number of drive wave superposed in the graphic
    nb_values = 400
    # Distance between sites of density 1 and nb_drive
    dist = dist_1 - dist_nb + 1   
    # Figure             
    abscisse = np.arange(-x_graph_values,1)*dx       
    fig, ax1 = plt.subplots(figsize=[-int(abscisse[0])/10, 5.25])
    ax2 = ax1.twinx()
    # Drive   
    for i in index_time[-nb_values:]:
        last_index = np.where(wt[i,:]>nb_drive)[0][0]
        first_index = last_index - int(x_graph_values) 
        ax2.semilogy(abscisse, drive[i, first_index:last_index+1], color="crimson", alpha = 0.3)   
    ax2.semilogy(0, 0, label="Drive", color="crimson", linewidth=2) 
    # Histogramme
    bins = [x - 0.5 for x in range(int(abscisse[0]),0)]
    ax1.hist(-dist[0,:], bins = bins, histtype = 'bar', label = ['dist erad'], density = True, color="yellowgreen", alpha = 0.6)
    # Wild-type
    ax2.semilogy(abscisse, wt[index_time[-1], first_index:last_index+1], label="Wild-type", color="cornflowerblue", linewidth=3)
    # Wild-type last individual
    index_first_ind = np.where(wt[index_time[-10], first_index:last_index+1]>0)[0][0]
    ax2.semilogy(abscisse[index_first_ind-1:index_first_ind+1], wt[index_time[-10], first_index+index_first_ind-1:first_index+index_first_ind+1], color="cornflowerblue", linewidth=3)#, label = "Last wild-type indi.")
    #ax2.hlines(nb_drive, (first_individual-30-last_index)*dx, 0, color="black", linestyle="-.") 
    #ax1.set_title(f"Extinction histogramme wt")
    ax2.set_ylim([None,ref_values[1]]); ax2.set_xlim([abscisse[0],0])
    ax1.xaxis.label.set_size(label_size)
    ax1.yaxis.label.set_size(label_size)   
    ax2.yaxis.label.set_size(label_size)   
    plt.rc('legend', fontsize=legend_size)  
    plt.legend()
    ax1.set_xlabel('Distances')
    ax2.set_ylabel('Number of individuals')
    ax1.set_ylabel('Histogramme proportions')
    ax1.set_ylim([0,0.3])
    plt.tight_layout() 
    fig.savefig(f"{dir_save}/histogram_s_{s}.png", format='png')
    fig.savefig(f"{dir_save}/histogram_s_{s}.svg", format='svg')
    plt.show()
    


    
     
    

### Parameters and datas ###
    
# Choose some parameters 
nb_drive = 100              # Threshold for drive ("enough" drive for the chasing not to happen)
conv_timing = "ger"         # Conversion timing : "ger" or "zyg"
dx = 1                      # spatial step
K = 10**8                    # Carrying capacity on one space unit
s = 0.3                     # Disadvantage for drive
r = 0.1                     # Intrasic growth rate


# Load the other parameters
dir_load = f"../../../stoch_not_save/datas/{conv_timing}_K_{int(np.log10(K))}_dx_{dx}_s_{s}_r_{r}"
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

load = True
if load: 
    time = np.loadtxt(f"{dir_load}/time.txt")              # Vector time
    wt = np.loadtxt(f"{dir_load}/nW_matrix.txt")           # Matrix wt : columns = times, row = space
    drive = np.loadtxt(f"{dir_load}/nD_matrix.txt")        # Matrix drive : columns = times, row = space


# Determine parameters from data
start = int(3*time[-1]/5)   # Time sequence to look at the exponential section (not to early because of CI)
end = int(time[-1]-10)      # Time sequence to look at the exponential section (not to late because of the section being partly outside the window)
index_time = np.intersect1d(np.where(time>start)[0], np.where(time<end)[0])
    
    
# Parameters for figures
if K > 10**5: ref_values = [0.001*K*dx]*2   #  Reference position for large K
if K <= 10**5: ref_values = [0.01*K*dx]*2   #  Reference position for small K
nb_ind = [10, 10]
nb_graph = 2       # Number of graphs shown


# Where to save results
dir_save = f"../outputs/K_{int(np.log10(K))}_s_{s}"
if not os.path.exists(dir_save): os.mkdir(dir_save)


# Preliminar
posi_pic, posi_1, posi_nb, dist_1, dist_nb = dist(time, index_time, wt, drive, K, ref_values, dir_save)


#plot_wave(time, index_time, dx, nb_graph, wt, drive, ref_values, dir_save)
v_num, [lWbn, lDbn], L_num = num(index_time, time, dx, wt, drive, nb_ind, nb_drive, posi_pic, ref_values, dist_1)
v_cont, lDfc, [lWbc, lDbc] = continu(conv_timing,s,h,c,r); L_cont = L(ref_values, [lWbc, lDbc]) 
v_dis, lDfd, [lWbd, lDbd] = discrete(conv_timing,s,h,c,r,m,dx,dt); L_dis = L(ref_values, [lWbd, lDbd]) 


# Plot the end of the wave (exp section)
x_graph_values = 160/dx
chasing = plot_end_wave(wt, drive, nb_drive, ref_values, index_time, [lWbd, lDbd], dx, s, x_graph_values, dir_save)

if not chasing: 
    # Eradication time in the wave
    erad_wt = eradication_time_mvt(time, index_time, wt, nb_drive)
    
    # Eradication time galton-watson
    nb_i = 1000
    x_graph_values = 50/dx
    if K == 10**8: decal = plot_end_wave_wt(wt, index_time, [lWbd, lDbd], dx, s, x_graph_values, 10**6, dir_save)
    exposants, extinction_list, exp = data_histogram_gw(nb_sites, lWbd, nb_drive, 0, nb_i, T, dt, dx, r, s, c, h, m, x_graph_values, decal, dir_load, dir_save)
    histogram_gw(exposants, extinction_list, exp, x_graph_values)
        
    #Comparaison
    gw_wt = extinction_list[-1,:]
    comparaison_distance_time_gw(v_num, s, dist_1, dist_nb, erad_wt, gw_wt, 0, dir_save)
    
    #plot_distances(index_time, dist_1, dist_nb, directory, ref_values, nb_drive)
    x_graph_values = 160/dx
    plot_histo_on_wave(wt, drive, index_time, nb_drive, dist_1, dist_nb, s, x_graph_values, dir_save)
    
    
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



#fig, ax = plt.subplots()
#ax.plot(np.exp(np.linspace(-4,0,100)), color = "black")
#ax.set_ylim([0,1]); ax.set_xlim([0,100])
#plt.tight_layout()     
#fig.savefig(f"exp.svg", format='svg')
#plt.show()



#p_star = ((1-h)*(1+c)-1)/((1-2*h))
#mean_fitness = 2*(1-h)*p_star*(1-p_star)+(1-p_star)**2   
#(1-mean_fitness)/mean_fitness      