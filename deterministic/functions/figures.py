#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 01:49:29 2022

@author: lena
"""

# Load libraries
import numpy as np
import matplotlib.pyplot as plt


s = 0.25  # when selection only acts on survival with a fitness cost s 
h = 0.1    # and sh for heterozygous individuals
c = 0.25              # homing rate
homing = "zygote"   # "zygote" or "germline"
label_size = 17


s_1 = c/(1-h*(1-c))   
if homing == "zygote" :
    s_2 = c/(2*c + h*(1-c))
    A_str = s*(2*(1-c)*(1-h)-1)
    B_str = s*((1-(1-c)*(1-h)))
    p_star = (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h)))       
    if s_1 < s_2 : 
        print("\nCoexistence", ":", np.round(s_1,3), "< s <", np.round(s_2,3))
    else :  
        print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        if c*(1-2*s)-(1-c)*s*h > 0 : 
            print("Linear speed :", 2*np.sqrt(c*(1-2*s)-(1-c)*s*h))                  
if homing == "germline" :         
    s_2 = c/(2*c*h + h*(1-c))
    A_str = s*(1-2*h)
    B_str = s*h
    p_star = ((1-s*h)*(1+c)-1)/(s*(1-2*h))
    if s_1 < s_2 :
        print("\nCoexistence", ":", np.round(s_1,3) , "< s <", np.round(s_2,3))
    else :  
        print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
    if c*(1-2*s*h)-(1-c)*s*h > 0 : 
        print("Linear speed :", 2*np.sqrt(c*(1-2*s*h)-(1-c)*s*h))

print("p_star :", np.round(p_star,2))   
print("s_1 :", np.round(s_1,2))   
print("s_2 :", np.round(s_2,2))   
    
p = np.linspace(-1,2,300)
f = (-A_str*(p-p_star)*p*(1-p))/(-A_str*p**2 - 2*B_str*p + 1)

fig, ax = plt.subplots()
ax.plot(p, f)
plt.hlines(y=0, color='dimgray', xmin=-1, xmax=2)   
plt.vlines(x=0, color='dimgray', ymin=-0.1, ymax=0.1)   
plt.vlines(x=1, color='dimgray', ymin=-0.1, ymax=0.1, linestyles ="dashed")     
ax.set(xlabel='p', ylabel='f(p)', xlim=[-1,2], ylim=[-0.02,0.02]) 
ax.xaxis.label.set_size(label_size); ax.yaxis.label.set_size(label_size)       
#fig.savefig(f"../outputs/{new_dir}/poly_{homing}_s_{s}_h_{h}_c_{c}.svg", format='svg')
plt.show()

h = np.linspace(0.01,0.99,100)
fig, ax = plt.subplots()
ax.plot(h, (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h))), label="zygote")
ax.plot(h, ((1-s*h)*(1+c)-1)/(s*(1-2*h)), label="germline" )
plt.hlines(y=0, color='dimgray', xmin=0, xmax=1)   
plt.hlines(y=-1, color='dimgray', xmin=0, xmax=1)  
plt.hlines(y=1, color='dimgray', xmin=0, xmax=1)    
ax.set(xlabel='h', ylabel='p_star', xlim=[0,1], ylim=[-3,3]) 
ax.xaxis.label.set_size(label_size); ax.yaxis.label.set_size(label_size) 
ax.legend()         
#fig.savefig(f"../outputs/{new_dir}/poly_{homing}_s_{s}_h_{h}_c_{c}.svg", format='svg')
plt.show()

