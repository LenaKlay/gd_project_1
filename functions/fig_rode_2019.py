#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 4 15:32:47 2021

@author: lena
"""


import numpy as np
import matplotlib.pyplot as plt


for c in np.linspace(0.1,0.9,9) :
    h=np.linspace(0.1,1,9)
    fig, ax = plt.subplots()
    ax.plot(h,c/(1-h*(1-c)), label='s1')
    #ax.plot(h, c/(2*c + h*(1-c)), label='s2 zyg')
    ax.plot(h, c/(2*c*h + h*(1-c)), label='s2 ger')
    ax.vlines(0.4,0,1)
    ax.vlines(0.6,0,1)
    #ax.plot(h,2*np.sqrt(c*(1-2*s*h)- s*h*(1-c)))
    ax.set(xlabel='h value', ylabel='s value', ylim = (0,1), title=f"c={c}")
    plt.legend()
    plt.show()


for c in np.linspace(0.1,0.9,9) : 
    c = np.round(c,2)
    homing = 'germline'      
    precision = 500  
    res = np.zeros((precision,precision))
    values = np.linspace(0.01,0.99,precision)
    for s_index in range(precision) :
        for h_index in range(precision) :
            s = values[s_index]
            h = values[h_index]
            s_1 =  c/(1-h*(1-c))   
            if homing == "zygote" : s_2 = c/(2*c + h*(1-c))
            if homing == "germline" : s_2 = c/(2*c*h + h*(1-c))
            if s < min(s_1,s_2) : 
                res[h_index, s_index] = 1 # Drive monostable
            elif s > max(s_1,s_2) : 
                res[h_index, s_index] = 4 # Wild monostable
            elif s_1 < s_2 : 
                res[h_index, s_index] = 2 # Coexistence
            elif s_2 < s_1 : 
                res[h_index, s_index] = 3  # Bistability
            else : 
                print('problème...')
            
    fig, ax = plt.subplots()
    im = ax.imshow(res)
    ax.figure.colorbar(im, ax=ax)
    ax.set_title(f"c = {c}")
    plt.gca().invert_yaxis()  
    fig.tight_layout()
    plt.show()
    
    
    
# trouver les bons couples
    
for c in np.linspace(0.1,0.9,17) :
    for h in np.linspace(0.1,0.9,17) : 
        c = np.round(c,2)
        h = np.round(h,2)
        if h < 0.5 : 
                s1 = c/(1-h*(1-c))
                s3 = c/(2*c + h*(1-c))
                s2 = c/(2*c*h + h*(1-c))               
                if s1 > 0.2 and s1 < 0.8 and  s2 > 0.2 and s2 <0.8 and abs(s1-s2) > 0.2 :
                    print('Coexistance c :',c,'and h :',h)
                                    
        if h > 0.5 :
                s1 = c/(1-h*(1-c))
                s3 = c/(2*c + h*(1-c))
                s2 = c/(2*c*h + h*(1-c))               
                if s1 > 0.2 and s1 < 0.8 and  s2 > 0.2 and s2 <0.8 and s3 > 0.2 and abs(s1-s2) > 0.2 :
                    print('Bistabilite c :',c,'and h :',h)
            
    

for c in np.linspace(0.1,0.9,9) :
    h=np.linspace(0.1,1,9)
    fig, ax = plt.subplots()
    ax.plot(h,c/(1-h*(1-c)), label='s1')
    #ax.plot(h, c/(2*c + h*(1-c)), label='s2 zyg')
    ax.plot(h, c/(2*c*h + h*(1-c)), label='s2 ger')
    ax.vlines(0.4,0,1)
    ax.vlines(0.6,0,1)
    #ax.plot(h,2*np.sqrt(c*(1-2*s*h)- s*h*(1-c)))
    ax.set(xlabel='h value', ylabel='s value', ylim = (0,1), title=f"c={c}")
    plt.legend()
    plt.show()