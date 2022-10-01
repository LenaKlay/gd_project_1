#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 4 15:32:47 2021

@author: lena
"""


import numpy as np
import matplotlib.pyplot as plt

    
# trouver les bons couples
    
homing = "germline"
for h in np.linspace(0.1,0.9,17) : 
    for c in np.linspace(0.1,0.9,17) :
        c = np.round(c,2)
        h = np.round(h,2)
        s1 = c/(1-h*(1-c))
        s2 = c/(2*c + h*(1-c))    # zygote
        s3 = c/(2*c*h + h*(1-c))  # germline
        
        if homing == "zygote" :
            if (1-h)*(1-c) > 0.5 :
                if s1 > 0.2 and s1 < 0.8 and  s2 > 0.2 and s2 <0.8 and abs(s1-s3) > 0.2 :
                    print('Coexistance c :',c,'and h :',h)                                    
            if (1-h)*(1-c) < 0.5 :                          
                if s1 > 0.2 and s1 < 0.8 and  s2 > 0.2 and s2 <0.8 and abs(s1-s2) > 0.2 :
                    print('Bistabilite c :',c,'and h :',h)
        
        if homing == "germline" :
            if h < 0.5 :
                if s1 > 0.2 and s1 < 0.8 and  s3 > 0.2 and s3 <0.8 and abs(s1-s3) > 0.2 :
                    print('Coexistance c :',c,'and h :',h, abs(s1-s3))           
            if h > 0.5 :                           
                if s1 > 0.2 and s1 < 0.8 and  s3 > 0.2 and s3 <0.8 and abs(s1-s3) > 0.2 :
                    print('Bistabilite c :',c,'and h :',h, abs(s1-s3))

# Zygote 
#c = 0.25; h = 0.1           
#c = 0.75; h = 0.1    OU   c = 0.25; h = 0.75

# Germline  
#c = 0.25; h = 0.3            
#c = 0.25; h = 0.75            
       
    
    
    
    

# Divers

for i in range(3) : 
    c = [0.25, 0.75, 0.25][i]
    homing = ['zygote', 'zygote', 'germline'][i]      
    precision = 1000  
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
    ax.set_title(f"{homing} c = {c}")
    ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
    ax.set_ylabel("h (dominance)", fontsize=12)
    plt.gca().invert_yaxis()  
    fig.tight_layout()
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
