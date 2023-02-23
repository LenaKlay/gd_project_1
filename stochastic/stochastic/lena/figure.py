#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 17:15:16 2022

@author: lena
"""

# Librairies
import numpy as np
import matplotlib.pyplot as plt


# Fonctions 

def plot_wave(index_time, time, dx, nb_graph, wt, drive):  
    for i in index_time[np.linspace(0,len(index_time)-1,nb_graph).astype(int)]:  
        
        # drive frequency
        drive_freq = -1*np.ones(np.shape(wt)[1])
        index_pop = np.where(wt[i,:]+drive[i,:]!=0)[0]
        drive_freq[index_pop] = drive[i,index_pop]/(wt[i,index_pop]+drive[i,index_pop])
        # plot
        fig, ax = plt.subplots()
        ax.plot(np.arange(nb_sites)*dx, wt[i,:], label = "wt")
        ax.plot(np.arange(nb_sites)*dx, drive[i,:], label = "drive")
        ax.set(xlabel='Space', ylabel='Number of individuals', ylim = [0,1.1*K*dx])
        ax.set_title(f"t = {time[i]}")
        plt.legend()
        fig.savefig(f"t_{time[i]}_{K}.png", format='png')
        plt.show()
        
        fig, ax = plt.subplots()
        first_wt, last_wt = section_exp(wt[i,:])
        ax.plot(np.arange(nb_sites)[first_wt:last_wt]*dx, drive_freq[first_wt:last_wt], label = "p")
        ax.set_title(f"wt, t = {time[i]}")
        plt.legend()
        plt.show()
        
        fig, ax = plt.subplots()
        first_drive, last_drive = section_exp(drive[i,:])
        ax.plot(np.arange(nb_sites)[first_drive:last_drive]*dx, drive_freq[first_drive:last_drive], label = "p")
        ax.set_title(f"drive, t = {time[i]}")
        plt.legend()
        plt.show()
        

def plot_distance_last_ind(index_time, time, dx, wt, drive, s, K):
    last_ind_list = np.zeros((2,len(time)))  # first line : position of the last wt individual, second line : position of the last drive individual
    pic_list = np.zeros(len(time))           # position of 20% of the drive pic
   
    for i in range(len(time)) :        
        pic_list[i] = np.where(drive[i,:]>pic_percent*max(drive[i,:]))[0][0]*dx 
        last_ind_list[0,i] = np.where(wt[i,:]!=0)[0][0]*dx
        last_ind_list[1,i] = np.where(drive[i,:]!=0)[0][0]*dx 

    # distances function of time
    fig, ax = plt.subplots()
    ax.plot(time[index_time], pic_list[index_time]-last_ind_list[0,index_time], label = "wt")
    ax.plot(time[index_time], pic_list[index_time]-last_ind_list[1,index_time], label = "drive")
    ax.set(xlabel='Time', ylabel='Distance')
    ax.set_title("Distance from the last individual to 20 percent of the pic")
    plt.legend()
    fig.savefig(f"distance_time_{K}.png", format='png')
    plt.show()

    # histogram of those same values
    wt_hist = pic_list[index_time]-last_ind_list[0,index_time]
    drive_hist = pic_list[index_time]-last_ind_list[1,index_time]
    fig, ax = plt.subplots()
    ax.hist([wt_hist, drive_hist], bins = 40, histtype = 'bar', label = ['wt','drive'])
    ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
    ax.set_title("Distance from the last individual to 20 percent of the pic")
    plt.legend()
    fig.savefig(f"distance_histogram_{K}.png", format='png')
    plt.show()
    
    # print mean and variances
    print("mean wt :", np.round(np.mean(wt_hist),3)); print("mean drive :", np.round(np.mean(drive_hist),3))
    print("variance wt :", np.round(np.var(wt_hist),3)); print("variance drive :", np.round(np.var(drive_hist),3))
    
    # save speed value
    speed = np.round(np.mean(np.diff(pic_list[index_time]/np.diff(time)[-1])),3)

    return(speed)
 

# return the section increasing exponentially (index and values). In this section we want strict positivity and croissance.  
def section_exp(vect):   
    # last_exp : where we last the exponential section
    last_index = np.where(vect>pic_percent*max(vect))[0][0]
    # first_exp : where we start the exponential section
    posit_index = np.where(vect[:last_index]==0)[0][-1]+1    # index for strict positivity
    first_index = posit_index   #(posit_index + last_index)//2                # we take only the last half of the value so that we have a quite stable section
    return(first_index, last_index)

# return the approximate exponential coefficient for the exponentially increasing section
def find_q_value(index_time, time, matrix):
    q_list = np.zeros(len(index_time))
    for i in index_time :
        vect = matrix[i,:]
        first_exp, last_exp = section_exp(vect)
        if last_exp-first_exp > 1 : 
            q_list[i-index_time] = np.exp(np.polyfit(np.arange(last_exp-first_exp), np.log(vect[first_exp:last_exp]), 1)[0])
    q_list = np.delete(q_list, np.where(q_list==0))
    q = np.mean(q_list)
    return(q, np.log(q))  
    
    
# return coefficients A and B for a wave of the form A*exp(lambda+) + B*exp(lambda-)  (behind the front, in the exponential section)
def coef_exp(speed, K, order_0, pic20_value, L):    
    lambda_pos = (-speed + np.sqrt(speed**2+4*order_0))/2
    lambda_neg = (-speed - np.sqrt(speed**2+4*order_0))/2  
    #L = np.log(K)/lambda_pos     # on ne peut pas diviser par 0 ! 
    A = pic20_value * np.exp(-L*lambda_pos)/( np.exp(-L*lambda_pos) - np.exp(-L*lambda_neg))
    B = pic20_value - A
    return(A, lambda_neg, B, lambda_pos)
    
    
    
# plot the exponentially increasing section but renormalized
def plot_sinus(index_time, time, dx, matrix, q, allele, speed, s):
    
    if allele == "wt" : order_0 = 1
    if allele == "drive" : order_0 = s
    
    fig, ax = plt.subplots()
    for i in index_time[np.linspace(0,len(index_time)-1,nb_sinus).astype(int)]: 
        # vect = allele repartition (wt or drive) at one precise time
        vect = matrix[i,:]
        # the exponentional section being in between first_exp and last_exp 
        first_exp, last_exp = section_exp(vect)   
        # renormalization for the exponential section
        suite_geo = [q**n for n in range(last_exp-first_exp)]
        # abscisse for plots
        abscisse = np.arange(first_exp-last_exp,0)*dx
        # plot exponential section
        ax.plot(abscisse, vect[first_exp:last_exp])
        # plot renormalized exponential section
        ax.plot(abscisse, vect[first_exp:last_exp]*suite_geo[::-1])
        # fit an exponential curve
        A, lambda_neg, B, lambda_pos = coef_exp(speed, K, order_0, vect[last_exp], (last_exp-first_exp)*dx)        
        ax.plot(abscisse, A*np.exp(abscisse*lambda_neg) + B*np.exp(abscisse*lambda_pos), color="black")
        ax.plot(abscisse, A*np.exp(abscisse*(lambda_neg-lambda_pos)) + B, color="gray")
    ax.set(xlabel='Space', ylabel='Densities')
    ax.set_title(f"Density and renormalized density in the exponential section. \n-{allele}-")
    fig.savefig(f"{allele}_renormalized_{K}.png", format='png')
    plt.show()
    
    

    

# Parameters
K = 1000000000    
s = 0.4   
dx = 0.1         
pic_percent = 0.2
nb_graph = 2      # Number of graphs shown
nb_sinus = 30     # Number of time values used in the sinus graph
start = 60        # Time sequence to look at the exponential section (not to early because of CI)
last = 100         # Time sequence to look at the exponential section (not to late because of the section being partly outside the window)


# Load data Florence  
#datacsv = open(f"flo_{K}/result_{K}.csv")
#datamatrix = np.loadtxt(datacsv, delimiter=",")
#datamatrix = datamatrix[:-3,:]   # suppress last values (where the pic goes behond the window)
#time = datamatrix[:,0]
#wt = datamatrix[:,np.arange(1,nb_sites*2,2)]
#drive = datamatrix[:,np.arange(1,nb_sites*2,2)+1]

# Load data Léna
time = np.loadtxt(f"lena_{K}_dx_{dx}/time.txt")
wt = np.loadtxt(f"lena_{K}_dx_{dx}/nW_matrix.txt")
drive = np.loadtxt(f"lena_{K}_dx_{dx}/nD_matrix.txt")
# Determine some parameters from data
index_time = np.intersect1d(np.where(time>start)[0], np.where(time<last)[0])
nb_sites = np.shape(wt)[1]


file = open(f"figures_values.txt", "w") 
speed = plot_distance_last_ind(index_time, time, dx, wt, drive, s, K); print("\nspeed :", speed); file.write(f"Speed : {speed}")  
plot_wave(index_time, time, dx, nb_graph, wt, drive)

q_wt, q_wt_log = find_q_value(index_time, time, wt); print(f"wt : q = {np.round(q_wt,3)} and log(q) = {np.round(q_wt_log,3)} and lamba_pos : {np.round((-speed + np.sqrt(speed**2+4))/2,3)}"); file.write(f"\nwt : q = {np.round(q_wt,3)} and log(q) = {np.round(q_wt_log,3)}") 
q_drive, q_drive_log = find_q_value(index_time, time, drive); print(f"drive : q = {np.round(q_drive,3)} and log(q) = {np.round(q_drive_log,3)} and lamba_pos : {np.round((-speed + np.sqrt(speed**2+4*s))/2,3)}"); file.write(f"\ndrive : q = {np.round(q_drive,3)} and log(q) = {np.round(q_drive_log,3)}") 

plot_sinus(index_time, time, dx, wt, q_wt, "wt", speed, s)
plot_sinus(index_time, time, dx, drive, q_drive, "drive", speed, s)

file.close() 




