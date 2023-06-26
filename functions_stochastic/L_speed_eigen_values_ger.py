#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 18 22:37:47 2023

@author: lena
"""


### Librairies

import numpy as np
  



# Return the section increasing exponentially (index and values). In this section we want strict positivity and croissance.  
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
    L_mean = np.zeros(2)
    for allele in range(2):
        L_mean[allele] = np.round(np.mean(dist_1[allele,:]),3)               
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
    return(v, lambda_pos, L_mean)  
    


# Speed and lambdas (continuous)
def continu(conv_timing,s,h,c,r):
    # Conversion timing
    if conv_timing == "ger" :
        term = (1-s*h)*(1+c)
    if conv_timing == "zyg" :
        term = 2*c*(1-s) + (1-c)*(1-s*h) 
    # Speed 
    v = 2*np.sqrt(term-1)
    # Lambdas
    lDfc = -v/2
    lDbc = -v/2 + np.sqrt(term - (r+1)*(1-s))
    lWbc = -v/2 + np.sqrt(term - (r+1)*(1-s*h)*(1-c))
    return(v, lDfc, lDbc, lWbc)

    
# Speed and lambdas (discrete)
def discrete(conv_timing,s,h,c,r,m,dx,dt):
    # Conversion timing
    if conv_timing == "ger" :
        term = (1-s*h)*(1+c)
    if conv_timing == "zyg" :
        term = 2*c*(1-s) + (1-c)*(1-s*h) 
    # Lambda vector
    lambda_vect = np.linspace(0.001,2,20001)
    # Speed 
    v = np.min(np.log((1+(term-1)*dt)*(1-m+m*np.cosh(lambda_vect*dx)))/(lambda_vect*dt))
    # Drive lambda at the front
    eqDfd = np.exp(lambda_vect*v*dt) - (1-(term-1)*dt)*(1-m+m*np.cosh(lambda_vect*dx))
    lDfd = lambda_vect[np.where(abs(eqDfd)==np.min(abs(eqDfd)))[0][0]]
    # Drive and Wild-type lamdba at the back
    eqDbd = np.exp(-lambda_vect*v*dt) - (1+((r+1)*(1-s)-1)*dt)*(1-m+m*np.cosh(lambda_vect*dx))
    eqWbd = np.exp(-lambda_vect*v*dt) - (1+((r+1)*(1-s*h)*(1-c)-1)*dt)*(1-m+m*np.cosh(lambda_vect*dx))
    lDbd = lambda_vect[np.where(abs(eqDbd)==np.min(abs(eqDbd)))[0][0]]
    lWbd = lambda_vect[np.where(abs(eqWbd)==np.min(abs(eqWbd)))[0][0]]
    return(v, lDfd, lDbd, lWbd) 

    
# Distance of the last individual (density striclty positive)
def L(pic_values, lambda_pos):    
    L = [np.log(pic_values[0])/lambda_pos[0], np.log(pic_values[1])/lambda_pos[1]]
    return(L)
 




























    
    
    
    
    
    
    
    
    
    
    
    
    

    
          
