#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 17:15:16 2022

@author: lena
"""

### Librairies

import numpy as np
import matplotlib.pyplot as plt


### Functions 
 
### Preliminary statements       
def pic_dist(index_time, dx, wt, drive):
    # position of 20% of the drive pic
    pic_list = np.zeros(len(index_time))    
    # first line : distance of the last wt individual, second line : distance of the last drive individual   
    dist_list = np.zeros((2,len(index_time)))   
    for i in range(len(index_time)) :   
        j = index_time[i]
        pic_list[i] = np.where(drive[j,:]>pic_percent*max(drive[j,:]))[0][0]*dx 
        dist_list[0,i] = pic_list[i] - np.where(wt[j,:]!=0)[0][0]*dx
        dist_list[1,i] = pic_list[i] - np.where(drive[j,:]!=0)[0][0]*dx 
    # Mean and variances
    mean = np.zeros(2); var = np.zeros(2)
    for i in range(2):
        mean[i] = np.round(np.mean(dist_list[i,:]),3)
        var[i] = np.round(np.var(dist_list[i,:]),3)
    return(pic_list, dist_list, mean, var)

### Speeds and distances
def speed_num(time, pic_list):
    v_num = np.mean(np.diff(pic_list))/(time[-1]-time[-2])
    return(v_num)
    
def speed_dis1(s,m,dx,dt):
    lambda_pos = np.linspace(0.001,2,2001)
    v_dis1 = np.min(np.log((1+(1-2*s)*dt)*(1-m+m*np.cosh(lambda_pos*dx)))/(lambda_pos*dt))
    return(v_dis1)


def distance_continu(allele, matrix, speed, K):
    if allele == "wt" : order_0 = 1
    if allele == "drive" : order_0 = s
    lambda_pos = (-speed + np.sqrt(speed**2+4*order_0))/2
    L_cont = np.log(pic_percent*K)/lambda_pos
    return(L_cont)
    
def distance_semi_num(allele, matrix, index_time, speed, K, dx):
    if allele == "wt" : order_0 = 1
    if allele == "drive" : order_0 = s
    vect = matrix[index_time[-1],:]
    first_exp, last_exp = section_exp(vect)  
    L_semi_num = coef_exp(speed, K, order_0, vect[last_exp], (last_exp-first_exp)*dx)[0]
    return(L_semi_num)

def distance_dis1(allele, v_dis1, s, m, dt, dx):
    lambda_vect = np.linspace(0,2,2001)
    if allele == "drive" : eq = np.exp(-lambda_vect*v_dis1*dt) - (1-s*dt)*(1-m+m*np.cosh(lambda_vect*dx))
    if allele == "wt" : eq = np.exp(-lambda_vect*v_dis1*dt) - (1-dt)*(1-m+m*np.cosh(lambda_vect*dx))
    lambda_dis1 = lambda_vect[np.where(abs(eq)==np.min(abs(eq)))[0][0]]
    L_dis1 = np.log(pic_percent*K)/lambda_dis1
    return(L_dis1)


### Plots 
    
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
        
def plot_distances(index_time, pic_list, dist_list):    
    # distances function of time
    fig, ax = plt.subplots()
    ax.plot(time[index_time], dist_list[0,:], label = "wt")
    ax.plot(time[index_time], dist_list[1,:], label = "drive")
    ax.set(xlabel='Time', ylabel='Distance')
    ax.set_title("Distance from the last individual to {pic_percent}*pic")
    plt.legend()
    fig.savefig(f"distance_time_{K}.png", format='png')
    plt.show()   
    # histogram of those same values
    wt_hist = dist_list[0,:]
    drive_hist = dist_list[1,:]
    fig, ax = plt.subplots()
    ax.hist([wt_hist, drive_hist], bins = 40, histtype = 'bar', label = ['wt','drive'])
    ax.set(xlabel='Distances', ylabel='Frequencies of each distances')
    ax.set_title("Distance from the last individual to 20 percent of the pic")
    plt.legend()
    fig.savefig(f"distance_histogram_{K}.png", format='png')
    plt.show()   


# return the section increasing exponentially (index and values). In this section we want strict positivity and croissance.  
def section_exp(vect): 
   # last_exp : where we last the exponential section
    last_index = np.where(vect>pic_percent*max(vect))[0][0]
    # first_exp : where we start the exponential section
    if np.any(vect[:last_index]==0):
        posit_index = np.where(vect[:last_index]==0)[0][-1]+1    # index for strict positivity
        first_index = posit_index   #(posit_index + last_index)//2                # we take only the last half of the value so that we have a quite stable section
    else : first_index = "outside_window"
    return(first_index, last_index)

# return the approximate exponential coefficient for the exponentially increasing section
def find_q_value(index_time, time, matrix):
    q_list = np.zeros(len(index_time))
    for i in index_time :
        # vect = allele repartition (wt or drive) at time[i]
        vect = matrix[i,:]
        # index for begin and end of the exponential section considered
        first_exp, last_exp = section_exp(vect)
        if type(first_exp)==str : 
            break
        elif last_exp-first_exp > 1 : 
            # compute the exponential coefficient of vect[first_exp:last_exp]
            q_list[i-index_time] = np.exp(np.polyfit(np.arange(last_exp-first_exp), np.log(vect[first_exp:last_exp]), 1)[0])
    q_list = np.delete(q_list, np.where(q_list==0))
    q = np.mean(q_list)
    return(q, np.log(q))  
    
    
# return coefficients A and B for a wave of the form A*exp(lambda+) + B*exp(lambda-)  (behind the front, in the exponential section)
def coef_exp(speed, K, order_0, pic20_value, L):    
    lambda_pos = (-speed + np.sqrt(speed**2+4*order_0))/2
    lambda_neg = (-speed - np.sqrt(speed**2+4*order_0))/2  
    A = pic20_value * np.exp(-L*lambda_pos)/( np.exp(-L*lambda_pos) - np.exp(-L*lambda_neg))
    B = pic20_value - A
    # L : Theoritical lenght in between the last individual and the pic of 20%
    if abs(A) < 0.001 : L = np.log(pic_percent*K)/lambda_pos    
    if abs(B) < 0.001 : L = np.log(pic_percent*K)/lambda_neg   
    return(L, A, lambda_neg, B, lambda_pos)
    
      
# plot the exponentially increasing section renormalized
def plot_sinus(index_time, time, dx, matrix, q, allele, speed, s):
    
    if allele == "wt" : order_0 = 1
    if allele == "drive" : order_0 = s
    
    fig, ax = plt.subplots()
    for i in index_time[np.linspace(0,len(index_time)-1,nb_sinus).astype(int)]: 
        # vect = allele repartition (wt or drive) at time[i]
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
        L, A, lambda_neg, B, lambda_pos = coef_exp(speed, K, order_0, vect[last_exp], (last_exp-first_exp)*dx)
        # -> classic exponential curve (fit)
        ax.plot(abscisse, A*np.exp(abscisse*lambda_neg) + B*np.exp(abscisse*lambda_pos), color="black")   
        # -> normalized exponential curve (fit)
        ax.plot(abscisse, A*np.exp(abscisse*(lambda_neg-lambda_pos)) + B, color="gray")
    ax.set(xlabel='Space', ylabel='Densities')
    ax.set_title(f"Density and renormalized density in the exponential section. \n-{allele}-")
    fig.savefig(f"{allele}_renormalized_{K}.png", format='png')
    plt.show()

    return(L)
    

    

### Parameters and datas
    
# Choose K and dx
K = 1000000000    
dx = 0.1         

# Load the other parameters
file = open(f"lena_{K}_dx_{dx}/parameters.txt", "r")
para = file.read()
para_list = para.replace(' ', '').split("\n")[:-1]
print(para_list)
dt = float(para_list[5].replace('dt=', ''))
m = float(para_list[6].replace('m=', ''))
s = float(para_list[8].replace('s=', ''))
file.close()

# Load datas
time = np.loadtxt(f"lena_{K}_dx_{dx}/time.txt")              # Vector time
wt = np.loadtxt(f"lena_{K}_dx_{dx}/nW_matrix.txt")           # Matrix wt : columns = times, row = space
drive = np.loadtxt(f"lena_{K}_dx_{dx}/nD_matrix.txt")        # Matrix drive : columns = times, row = space

# Determine parameters from data
start = 300         # Time sequence to look at the exponential section (not to early because of CI)
last = 1000         # Time sequence to look at the exponential section (not to late because of the section being partly outside the window)
index_time = np.intersect1d(np.where(time>start)[0], np.where(time<last)[0])
nb_sites = np.shape(wt)[1]

# Parameters for figures
pic_percent = 0.2
nb_graph = 4       # Number of graphs shown
nb_sinus = 30      # Number of time values used in the sinus graph
show_graph = False



### Results
def p(text,num): print(text,np.round(num,5))

# Preliminaries
pic, dist, mean, var = pic_dist(index_time, dx, wt, drive)

# Plot the wave and its speed, and distances of the last individuals to 20% of the pic.

if show_graph : 
    plot_wave(index_time, time, dx, nb_graph, wt, drive)
    plot_distances(index_time, pic, dist)
    
# Speeds 
v_cont = np.round(2*np.sqrt(1-2*s),3); p("\nSpeed cont :", v_cont)
v_num = speed_num(time, pic); p("Speed num :", v_num)
v_dis1 = speed_dis1(s,m,dx,dt); p("Speed dis1 :", v_dis1)
# Distances
allele = "wt"; matrix = wt
L_wt_cont = distance_continu(allele, matrix, v_num, K); p("\nL wt cont :", L_wt_cont)
p("L wt num :", mean[0])
L_wt_semi_num = distance_semi_num(allele, matrix, index_time, v_num, K, dx); p("L wt semi num :", np.round(L_wt_semi_num,3))
L_wt_dis1 = distance_dis1(allele, v_dis1, s, m, dt, dx); p("L wt dis1 :", np.round(L_wt_dis1,3))
allele = "drive"; matrix = drive
L_drive_cont = distance_continu(allele, matrix, v_num, K); p("\nL drive cont :", L_drive_cont)
p("L drive num :",mean[1])
L_drive_semi_num = distance_semi_num(allele, matrix, index_time, v_num, K, dx); p("L drive semi num :", L_drive_semi_num)
L_drive_dis1 = distance_dis1(allele, v_dis1, s, m, dt, dx); p("L drive dis1 :", L_drive_dis1)


    

# Find numerically the exponential coefficient.
if show_graph :
    q_wt, q_wt_log = find_q_value(index_time, time, wt); print(f"wt : q = {np.round(q_wt,3)} and log(q) = {np.round(q_wt_log,3)} and lamba_pos : {np.round((-v_num + np.sqrt(v_num**2+4))/2,3)}")
    q_drive, q_drive_log = find_q_value(index_time, time, drive); print(f"drive : q = {np.round(q_drive,3)} and log(q) = {np.round(q_drive_log,3)} and lamba_pos : {np.round((-v_num + np.sqrt(v_num**2+4*s))/2,3)}")

# Plot the renormalized graphs.
if show_graph :
    plot_sinus(index_time, time, dx, wt, q_wt, "wt", v_num, s)
    plot_sinus(index_time, time, dx, drive, q_drive, "drive", v_num, s)






# Save the results in the file "figures_values.txt"
#file = open(f"figures_values.txt", "w") 
#file.write(f"Speed num : {v_num}")
#file.write(f"\nwt : q = {np.round(q_wt,3)} and log(q) = {np.round(q_wt_log,3)}") 
#file.write(f"\ndrive : q = {np.round(q_drive,3)} and log(q) = {np.round(q_drive_log,3)}") 
#file.close() 




# BROUILLON

# CONTINU
#la_d = - np.sqrt(1-2*s)+np.sqrt(1-s)
#la_wt = - np.sqrt(1-2*s)+np.sqrt(2*(1-s))

# DISCRET
#speed_la = np.log((1+(1-2*s)*dt)*(1-m+m*np.cosh(la*dx)))/(la*dt)
#fig, ax = plt.subplots()
#ax.plot(la, speed_la)  
#plt.show() 




