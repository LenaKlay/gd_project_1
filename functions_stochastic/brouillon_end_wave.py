#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 18:01:36 2023

@author: lena
"""

def plot_end_wave_bis(wt, drive, nb_drive, pic_values, dx, i, al, ax, legend):
    last_index = max(np.where(drive[i,:]> pic_values[1])[0][0], np.where(wt[i,:]>100)[0][0])
    first_individual = min([max([np.where(drive[i,:]>10)[0][0], nb_sites//2]), np.where(wt[i,:]>0)[0][0]])
    abscisse = np.arange(first_individual-20-last_index,0)*dx
    if legend:       
        ax.semilogy(abscisse, drive[i, first_individual-20:last_index], color="crimson", alpha=al, label="Drive")
        ax.semilogy(abscisse, wt[i, first_individual-20:last_index], color="cornflowerblue", alpha=al, label="Wild-type")
        return([wt[i, last_index-1], drive[i, last_index-1]], abscisse) 
    else: 
        ax.semilogy(abscisse, drive[i, first_individual-20:last_index], color="crimson", alpha=al) 
        ax.semilogy(abscisse, wt[i, first_individual-20:last_index], color="cornflowerblue", alpha=al)        
        return([wt[i, last_index-1], drive[i, last_index-1]]) 


def plot_end_wave(wt, drive, nb_drive, pic_values, index_time, lambda_back, dx, s, dir_save):  
    fig, ax = plt.subplots()
    density = np.zeros((2,len(index_time))).astype('int')
    # Superposition of the numerical waves
    for i in index_time[:-1] :
        density[:, i-index_time[0]] = plot_end_wave_bis(wt, drive, nb_drive, pic_values, dx, i, 0.3, ax, False)
    density[:,-1], abscisse = plot_end_wave_bis(wt, drive, nb_drive, pic_values, dx, index_time[-1], 1, ax, True)
    # Theoritical exponential function
    beta = np.mean(density[1,:])/np.exp(lambda_back[1]*abscisse[-1]); exp_drive = beta*np.exp(lambda_back[1]*abscisse)
    gamma = np.mean(density[0,:])/np.exp(lambda_back[0]*abscisse[-1]); exp_wt = gamma*np.exp(lambda_back[0]*abscisse)
    exp_D1 = np.where(exp_drive>1)[0][0]; exp_WT1 = np.where(exp_wt>1)[0][0]
    ax.semilogy(abscisse[exp_D1:], exp_drive[exp_D1:], color="black")
    ax.semilogy(abscisse[exp_WT1:], exp_wt[exp_WT1:], color="black")           
    ax.set_title(f"Wild-type and drive densities at the back of the wave.")
    ax.set_ylabel('Densities'); ax.set_xlabel('Space') 
    plt.legend()
    fig.savefig(f"{dir_save}/end_wave_s_{s}.png", format='png')
    plt.show()