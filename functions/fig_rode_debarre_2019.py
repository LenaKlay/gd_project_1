#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 4 15:32:47 2021

@author: lena
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

    
# trouver les bons couples
    
conversion_timing = "germline"
for h in np.linspace(0.1,0.9,17) : 
    for c in np.linspace(0.1,0.9,17) :
        c = np.round(c,2)
        h = np.round(h,2)
        s1 = c/(1-h*(1-c))
        s2 = c/(2*c + h*(1-c))    # zygote
        s3 = c/(2*c*h + h*(1-c))  # germline
        
        if conversion_timing == "zygote" :
            if (1-h)*(1-c) > 0.5 :
                if s1 > 0.2 and s1 < 0.8 and  s2 > 0.2 and s2 <0.8 and abs(s1-s3) > 0.2 :
                    print('Coexistance c :',c,'and h :',h)                                    
            if (1-h)*(1-c) < 0.5 :                          
                if s1 > 0.2 and s1 < 0.8 and  s2 > 0.2 and s2 <0.8 and abs(s1-s2) > 0.2 :
                    print('Bistabilite c :',c,'and h :',h)
        
        if conversion_timing == "germline" :
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
    conversion_timing = ['zygote', 'zygote', 'germline'][i]      
    precision = 1000  
    res = np.zeros((precision,precision))
    values = np.linspace(0.01,0.99,precision)
    for s_index in range(precision) :
        for h_index in range(precision) :
            s = values[s_index]
            h = values[h_index]
            s_1 =  c/(1-h*(1-c))   
            if conversion_timing == "zygote" : s_2 = c/(2*c + h*(1-c))
            if conversion_timing == "germline" : s_2 = c/(2*c*h + h*(1-c))
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
    cmap = ListedColormap(["#ff7c7ce2", "#ea77eba1",  "#7ab552a2", "#bbcfffff"])  
    im = ax.imshow(res, cmap=cmap)
    ax.set_xticks(np.linspace(0,precision,5)); ax.set_yticks(np.linspace(0,precision,5))  
    ax.set_xticklabels(np.linspace(0,1,5)); ax.set_yticklabels(np.linspace(0,1,5))  
    ax.set_title(f"{conversion_timing} c = {c}")
    ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
    ax.set_ylabel("h (dominance)", fontsize=12)
    plt.gca().invert_yaxis()  
    fig.tight_layout()
    fig.savefig(f"../figures/rode_debarre_2019/{conversion_timing}_c_{c}.png", format='png')
    fig.savefig(f"../figures/rode_debarre_2019/{conversion_timing}_c_{c}.svg", format='svg')
    plt.show()
    

# Legend for color
x = np.linspace(0,1,10)
y1 = np.ones(10) 
y2 = np.ones(10)*2
y3 = np.ones(10)*3
y4 = np.ones(10)*4 
plt.plot(x,y4, color="#bbcfffff", label='wt')
plt.plot(x,y3, color="#7ab552a2", label='bist')
plt.plot(x,y2, color="#ea77eba1", label='coex')
plt.plot(x,y1, color="#ff7c7ce2", label='drive')
plt.legend()
plt.show()
    
  
    
    
    
    
    
