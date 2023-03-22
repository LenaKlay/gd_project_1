#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 16:34:21 2023

@author: lena
"""

#import matplotlib
#import colorsys
#colorNames = list(matplotlib.colors.cnames.keys())

import matplotlib.pyplot as plt
import numpy as np

# Font
plt.rcParams.update({'font.family':'serif'})
# Graph parameters 
label_size = 12 
line_size = 2


what_to_do = "density_nd" # "A_fct_of_a"  "ext_bist_line" "density_nd"

## Print the extinction and bistability line for all possible a values
if what_to_do == "ext_bist_line" : 
    cas = "d"; bist = False
    fig, ax = plt.subplots() 
    s = np.linspace(0, 1, 2000)
    acolor=['navy','darkblue','royalblue','cornflowerblue','lightskyblue','greenyellow','gold', 'orange','red','crimson','purple']
    if cas in ["a", "c"] :      
        ax.semilogy(s, s/(1-s), color='black', linewidth=line_size)
    if cas in ["b","d"] :
        for i in range(11):
            a = np.round(np.linspace(-1,1,11)[i],2)
            A = ((1-a)/2)**2  
            if cas == "b" :  
                if bist : 
                    ax.semilogy(s, -s/(a*(1-s)), label=f"a={a}", color=acolor[i], linewidth=line_size)        
                else :
                    ax.semilogy(s, s/(A*(1-s)), label=f"a={a}", color=acolor[i], linewidth=line_size)
            if cas == "d" :  
                if bist : 
                    s_under = s[np.where(s>-a)]; s_above = s[np.where(s<-a)]
                    ax.semilogy(s_under, -s_under/(a+s_under), label=f"a={a}", color=acolor[i], linewidth=line_size) 
                    ax.semilogy(s_above, -s_above/(a+s_above), color=acolor[i], linewidth=line_size)
                else :
                    s_under = s[np.where(s>A)]; s_above = s[np.where(s<A)]
                    ax.semilogy(s_under, s_under/(A-s_under), label=f"a={a}", color=acolor[i], linewidth=line_size)
                    ax.semilogy(s_above, s_above/(A-s_above), color=acolor[i], linewidth=line_size)            
    plt.legend()
    ax.set(xlabel='s (fitness disadvantage for drive)', ylabel="r (intrinsic growth rate)", xlim = (0,1), ylim = (0.01,10)) 
    ax.xaxis.label.set_size(label_size)
    ax.yaxis.label.set_size(label_size)     
    plt.grid()
    #fig.savefig(f"../outputs/extinction_cas_{cas}_bist_{bist}.png", format='png')
    plt.show()

# Value A = ((1-a)/2)**2 fct of a
if what_to_do == "ext_bist_line" :
    a_vect = np.linspace(-1, 1, 200)
    fig, ax = plt.subplots()
    ax.plot(a_vect, ((1-a_vect)/2)**2, label='A')
    ax.set(xlabel='s (fitness disadvantage for drive)', ylabel="r (intrinsic growth rate)" , ylim = (0,1)) 
    plt.legend()
    plt.show()


# Heatmap density n_D^*
if what_to_do == "density_nd" :
    a=-0.2
    for cas in ["a","b","c","d"] : 
        precision = 1000
        res = np.zeros((precision,precision))
        s_values = np.linspace(0.01,0.99,precision)
        r_values = np.logspace(-2, 1, num=precision)
        for s_index in range(precision) :
            for r_index in range(precision) :
                s = s_values[s_index]
                r = r_values[r_index]
                if cas == "a" :
                    res[r_index, s_index] = 1 - s/(r*(1-s))
                if cas == "b" :
                    res[r_index, s_index] = 0.5*(1+a+np.sqrt((1+a)**2-4*(a+(s/(r*(1-s))))))
                if cas == "c" : 
                    res[r_index, s_index] = 1 - s*(r+1)/r
                if cas == "d" : 
                    res[r_index, s_index] = 0.5*(1+a+np.sqrt((1+a)**2-4*(a+(s*(r+1)/r))))    
                if res[r_index, s_index] < 0 :
                    res[r_index, s_index] = np.nan
        fig, ax = plt.subplots()
        im = ax.imshow(res, cmap='YlGnBu', vmin=0, vmax=1)
        ax.figure.colorbar(im, ax=ax)  # Add the colorbar
        ax.set_xticks(np.linspace(0,precision,5)); ax.set_yticks(np.linspace(0,1,4)*(precision-1))    
        ax.set_xticklabels(np.linspace(0,1,5)); ax.set_yticklabels(np.around(np.logspace(-2, 1, num=4),2))  
        #ax.set_title(f"cas = {cas}, a = {a}, A = {A}")
        ax.set_xlabel("s (fitness disadvantage for drive)", fontsize=12) 
        ax.set_ylabel("r (intrinsic growth rate)", fontsize=12)
        plt.gca().invert_yaxis()  
        fig.tight_layout()
        #fig.savefig(f"../outputs/density_a_{a}_cas_{cas}.png", format='png')
        plt.show()
    
