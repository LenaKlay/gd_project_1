#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 22:37:47 2023

@author: lena
"""


### Librairies

import numpy as np

import matplotlib.pyplot as plt


# return coefficients A and B for a wave of the form A*exp(lambda+) + B*exp(lambda-)  (behind the front, in the exponential section)   
def coef_exp(wt, drive, index_time, first_exp, last_exp, vect, lambda_pos, lambda_neg, dx, allele, nb_ind, pic_values): 
    A = np.zeros(2); B = np.zeros(2)
    for allele in range(2): 
        matrix = [wt, drive][allele]
        # Density at the last time value
        vect = matrix[index_time[-1],:]
        # Exponential section at the last time value
        first_exp, last_exp = section_exp(vect, allele, nb_ind, pic_values)  
        # Lenght of the exponential section
        lenght_exp = (last_exp - first_exp)*dx 
        # Exponential coefficents
        B[allele] = (vect[first_exp] - vect[last_exp]*np.exp(-lenght_exp*lambda_pos[allele]))/( np.exp(-lenght_exp*lambda_neg[allele]) - np.exp(-lenght_exp*lambda_pos[allele]))
        A[allele] = vect[last_exp] - B[allele]
    return(A, B)
  

### Preliminar: return the section increasing exponentially (index and values). In this section we want strict positivity and croissance.  
def section_exp(vect, allele, nb_ind, pic_values): 
    # last_exp : where we last the exponential section
    last_index = np.where(vect>pic_values[allele])[0][0]
    # first_exp : where we start the exponential section
    if np.any(vect[:last_index]==0):
        # we ensure strict positivity 
        posit_index = np.where(vect[:last_index]==0)[0][-1]+1
        # and nb_ind on the first site or before
        nb_ind_index = np.where(vect[:last_index]>nb_ind[allele])[0][0]    
        # max to have both conditions       
        first_index = max(posit_index, nb_ind_index)
    else : first_index = "outside_window"
    return(first_index, last_index)



# Speed and lambdas (numeric)
def num(index_time, time, dx, wt, drive, nb_ind, nb_drive, posi_pic, pic_values, dist_1):  
    # Mean and variances of the distance of the last individual
    L_mean = np.zeros(2); L_var = np.zeros(2)
    for allele in range(2):
        L_mean[allele] = np.round(np.mean(dist_1[allele,:]),3)
        L_var[allele] = np.round(np.var(dist_1[allele,:]),3)       
    # Numerical speed  
    v = np.mean(np.diff(posi_pic[0,int(3*np.shape(posi_pic)[1]/4):]))/(time[-1]-time[-2])
    # Numerical positive lambda (approximate exponential coefficient for the exponentially increasing section)
    lambda_pos_list = np.zeros((2,len(index_time))); lambda_pos = np.zeros(2)
    for allele in range(2): 
        matrix = [wt, drive][allele]
        for i in index_time: 
            # vect = allele repartition (wt and drive) at time[i]
            vect = matrix[i,:]
            # index for begin and end of the exponential section considered
            first_exp, last_exp = section_exp(vect, allele, nb_ind, pic_values)
            # compute the exponential coefficient of vect[first_exp:last_exp]
            if first_exp != "outside_window":
                lambda_pos_list[allele,i-index_time] = np.polyfit(np.arange(first_exp-last_exp,0)*dx, np.log(vect[first_exp:last_exp]), 1)[0]
        # Mean of the values found
        lambda_pos[allele] = np.mean(np.delete(lambda_pos_list[allele,:], np.where(lambda_pos_list[allele,:]==0)))
    return(v, lambda_pos, L_mean, L_var)  
    

# Speed and lambdas (continuous)
def continu(s):
    # Speed 
    v = 2*np.sqrt(1-2*s)
    # Lamdba pos and neg
    lambda_pos = np.zeros(2); lambda_neg = np.zeros(2)
    for allele in range(2) : 
        order_0 = [1,s][allele]
        lambda_pos[allele] = (-v + np.sqrt(v**2+4*order_0))/2
        lambda_neg[allele] = (-v - np.sqrt(v**2+4*order_0))/2  
    return(v, lambda_pos, lambda_neg)
    
# Speed and lambdas (discrete)
def discrete(s,m,dx,dt):
    # Speed 
    lambda_front = np.linspace(0.001,2,20001)
    v = np.min(np.log((1+(1-2*s)*dt)*(1-m+m*np.cosh(lambda_front*dx)))/(lambda_front*dt))
    # Lamdba pos and neg
    lambda_vect = np.linspace(-2,2,40002)
    lambda_pos = np.zeros(2); lambda_neg = np.zeros(2)
    for allele in range(2) : 
        order_0 = [1,s][allele]
        eq = np.exp(-lambda_vect*v*dt) - (1-order_0*dt)*(1-m+m*np.cosh(lambda_vect*dx))
        lambda_neg[allele] = lambda_vect[np.where(abs(eq[:20001])==np.min(abs(eq[:20001])))[0][0]]
        lambda_pos[allele] = lambda_vect[20001 + np.where(abs(eq[20001:])==np.min(abs(eq[20001:])))[0][0]]
    return(v, lambda_pos, lambda_neg) 
        

# Return coefficients A and B for a wave of the form A*exp(lambda+) + B*exp(lambda-)  
# (behind the front, in the exponential section)   
def profil_A_B(lambda_pos, lambda_neg, pic_values) : 
    B = np.zeros(2)
    for allele in range(2):
        B = pic_values[allele] * np.exp(lambda_pos[allele])/(np.exp(lambda_pos[allele]) - np.exp(lambda_neg[allele]))
    A = pic_values - B
    return(A,B)
    
# Return coefficients A and B, and L (with a third equation)
def profil_L_A_B(lambda_pos, lambda_neg, pic_values, allele, dx) : 
    L = np.zeros(2); B = np.zeros(2)
    for allele in range(2):
        L_vect = np.linspace(0,200,20001)
        B_vect = 1/(np.exp(-L_vect*lambda_neg[allele])*(np.exp(dx*lambda_neg[allele])-np.exp(dx*lambda_pos[allele])))
        eq = np.log(-B_vect/(pic_values[allele]-B_vect))/(lambda_neg[allele]-lambda_pos[allele]) - L_vect
        L = L_vect[np.where(abs(eq)==np.min(abs(eq)))[0][0]]
        B = B_vect[np.where(abs(eq)==np.min(abs(eq)))[0][0]]
    A = pic_values - B
    return(L, A, B)  

# Distance of the last individual (density striclty positive)
def L1(pic_values, lambda_pos):    
    L = [np.log(pic_values[0])/lambda_pos[0], np.log(pic_values[1])/lambda_pos[1]]
    return(L)
 
# Distance of the last individual (density beign zero at one point)
def L2(pic_values, lambda_pos, lambda_neg, dx):  
    B = np.zeros(2); L = np.zeros(2)
    for allele in range(2) :
        pic = pic_values[allele]; l_pos = lambda_pos[allele]; l_neg = lambda_neg[allele]
        L_vect = np.linspace(0.01,200,20001)
        B_vect = (-pic*np.exp(-L_vect*l_pos)) / (np.exp(-L_vect*l_neg)-np.exp(-L_vect*l_pos)) 
        eq = B_vect*(np.exp((dx-L_vect)*l_neg) - np.exp((dx-L_vect)*l_pos)) - (1-pic*np.exp((dx-L_vect)*l_pos))         
        L[allele] = L_vect[np.where(abs(eq)==np.min(abs(eq)))[0][0]]
        B[allele] = B_vect[np.where(abs(eq)==np.min(abs(eq)))[0][0]]
    A = pic_values - B
    return(L, A, B)

    
# plot the exponentially increasing section renormalized
def plot_sinus(index_time, wt, drive, dx, lambda_pos, lambda_neg, L, A, B, nb_ind, pic_values, directory, nb_sinus):

    for allele in range(2): 
        matrix = [wt, drive][allele]
        title = ["wt", "drive"][allele]        
        fig, ax = plt.subplots()
        
        for i in index_time[np.linspace(0,len(index_time)-1,nb_sinus).astype(int)]: 
            # vect = allele repartition (wt or drive) at time[i]
            vect = matrix[i,:]
            # the exponentional section being in between first_exp and last_exp 
            first_exp, last_exp = section_exp(vect, allele, nb_ind, pic_values)          
            # abscisse for plots
            abscisse = np.arange(first_exp-last_exp,0)*dx
            # plot exponential section
            ax.plot(abscisse, vect[first_exp:last_exp])
            # renormalization for the exponential section
            suite_geo = [np.exp(lambda_pos[allele])**(n*dx) for n in range(last_exp-first_exp)]
            # plot renormalized exponential section
            ax.plot(abscisse, vect[first_exp:last_exp]*suite_geo[::-1])
            
        # Compute A and B if not given
        if type(A) == str : 
            A, B = coef_exp(first_exp, last_exp, vect, lambda_pos, lambda_neg, dx, allele, nb_ind, pic_values)
        # Exponential curve
        ax.plot(abscisse, A[allele]*np.exp(abscisse*lambda_pos[allele]) + B[allele]*np.exp(abscisse*lambda_neg[allele]), color="black", label = "exp curve")   
        # Normalized exponential curve
        ax.plot(abscisse, A[allele] + B[allele]*np.exp(abscisse*(lambda_neg[allele]-lambda_pos[allele])), color="gray", label = "renormalized exp curve")
        ax.set(xlabel='Space', ylabel='Densities', ylim = [0,pic_values[allele]*1.2])
        ax.set_title(f"Density and renormalized density in the exponential section. \n-{title}-")
        plt.legend()
        fig.savefig(f"{directory}/{title}_renormalized.png", format='png')
        plt.show()
        
        fig, ax = plt.subplots()
        ax.semilogy(abscisse, vect[first_exp:last_exp], label="data")
        ax.semilogy(abscisse, A[allele]*np.exp(abscisse*lambda_pos[allele]) + B[allele]*np.exp(abscisse*lambda_neg[allele]), color="black", label="exp curve")   
        ax.set_title(f"Density and fit in a log scale. \n-{title}-")
        plt.legend()
        fig.savefig(f"{directory}/{title}_fit.png", format='png')
        plt.show()
        


























    
    
    
    
    
    
    
    
    
    
    
    
    

    
          
