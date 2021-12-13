#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:28:55 2021

@author: lena
"""


# Libraries
import numpy as np


def which_system(W,H,D,epsilon,window,num_para) :
    # Parameters
    T,L,M,N,mod,theta = num_para

    system = [[] for _ in range(len(W))]
    for i in range(0,len(W)) :
        if (1-W-D)[i] >= W[i] :
            if (1-W-D)[i] >= D[i] :
                if W[i]>=D[i] : 
                    system[i].append(1)
                if W[i]<=D[i] : 
                    system[i].append(2)
            if (1-W-D)[i] <= D[i] : 
                if W[i]>=D[i] : 
                    system[i].append(3)
                if W[i]<=D[i] : 
                    system[i].append(4)
        if (1-W-D)[i] <= W[i] :
            if (1-W-D)[i] >= D[i] :
                if W[i]>=D[i] : 
                    system[i].append(5)
                if W[i]<=D[i] : 
                    system[i].append(6)
            if (1-W-D)[i] <= D[i] : 
                if W[i]>=D[i] : 
                    system[i].append(7)
                if W[i]<=D[i] : 
                    system[i].append(8)  
    system_bis = [[] for _ in range(len(W))]
    for i in range(0,len(W)) :
        if (1-W-D)[i]+epsilon >= W[i] :
            if (1-W-D)[i]+epsilon >= D[i] :
                if W[i]+epsilon>=D[i] : 
                    system_bis[i].append(1)
                if W[i]<=D[i]+epsilon : 
                    system_bis[i].append(2)
            if (1-W-D)[i] <= D[i]+epsilon : 
                if W[i]+epsilon>=D[i] : 
                    system_bis[i].append(3)
                if W[i]<=D[i]+epsilon : 
                    system_bis[i].append(4)
        if (1-W-D)[i] <= W[i]+epsilon :
            if np.minimum(1-W-D,D)[i]+epsilon >= D[i] :
                if W[i]+epsilon>=D[i] : 
                    system_bis[i].append(5)
                if W[i]<=D[i]+epsilon : 
                    system_bis[i].append(6)
            if (1-W-D)[i] <= D[i]+epsilon : 
                if W[i]+epsilon>=D[i] : 
                    system_bis[i].append(7)
                if W[i]<=D[i]+epsilon : 
                    system_bis[i].append(8) 
                    
    if ([2] in system) and ([7] in system) : 
        end_2 = N-system[::-1].index([2])
        start_7 = system.index([7])
        border = int((end_2+start_7)/2)

        if window == 'everywhere' :
            print(system,"\n")
            print(system_bis,"\n")
        if window == 'zoom border' : 
            print(system[border-75*int(N/L):border+75*int(N/L)],"\n") 
            print(system_bis[border-75*int(N/L):border+75*int(N/L)],"\n")
        if window == 'minimum' : 
            print(system[border-4*int(N/L):border+4*int(N/L)],"\n") 
     
        if window != 'minimum' :     
            print("End of system (2) :", end_2/N*L,"\n")
            print("Start of system (7) :", start_7/N*L,"\n")
    
        border = int(border/N*L)
        
        return(border,system,system_bis)
        
    else :
        print(system[0:100])
        return(None,system,system_bis)