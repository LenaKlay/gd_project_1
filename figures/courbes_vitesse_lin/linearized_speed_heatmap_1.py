#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 22:37:57 2022

@author: lena
"""
import numpy as np
import matplotlib.pyplot as plt
label_size = 14
# Change font to serif
plt.rcParams.update({'font.family':'serif'})

c=0.75
h=0.1
s = np.linspace(0.1,0.9,100001)
s1 = c/(1-h*(1-c))
s2 = c/(2*c + h*(1-c))

fig, ax = plt.subplots()
ax.plot(s,2*np.sqrt(c*(1-2*s)-(1-c)*s*h),  '#ff7f0e', label=r'$v^{*}_{+}(s,1/4,1/10)$')
ax.plot(s,-2*np.sqrt((1-s*h)*(1-c)/(1-s)-1), '#1f77b4', label=r'$v^{*}_{-}(s,1/4,1/10)$')
plt.axhline(y=0,linewidth=1, color = "gray")   
plt.axvline(x=s1,linewidth=1, linestyle='--', color = "gray")   
plt.axvline(x=s2, linewidth=1, linestyle='--', color = "gray")  
ax.set(xlabel='s (fitness disadvantage for drive)', ylabel='Speed of the wave') 
ax.xaxis.label.set_size(label_size); ax.yaxis.label.set_size(label_size)      
plt.legend()
fig.savefig(f"speed.png", format='png')
plt.show()


print("zygote")
print("s1", c/(1-h*(1-c)))   
print("s2", c/(2*c + h*(1-c)))