#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 23:16:54 2022

@author: lena
"""


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors 


# Change font to serif
plt.rcParams.update({'font.family':'serif'})


fig = plt.figure()
ax = fig.gca(projection='3d')


# Make data.
X = np.arange(-3, 3, 0.01)
Y = np.arange(-3, 3, 0.01)
X, Y = np.meshgrid(X, Y)
Z = np.exp(-X**2)*np.exp(-Y**2)


colors = plt.cm.Spectral(np.flip(np.linspace(0, 0.3, 128*2)))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z,cmap=mymap,linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(0, 1.1)

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)
fig.savefig(f"../outputs/2D.svg", format='svg')
plt.show()





title_size = 15
label_size = 17
legend_size = 12
number_y_size = 20

x = np.arange(-3, 3, 0.1)
y = np.exp(-x**2)
fig, ax = plt.subplots()
ax.plot(x, y, color = 'crimson', label='Drive', linewidth = 4)
plt.rc('legend', fontsize=legend_size)  
ax.legend(bbox_to_anchor=(0.553,1.13), ncol=2)
defaultylim = (-0.03,1.03)
ax.set(xlabel='Space', ylabel="Allele densities" , ylim = defaultylim)
ax.set_title(f"t = 0", fontsize = title_size, loc='right')
ax.xaxis.label.set_size(label_size)
ax.yaxis.label.set_size(label_size)  
ax.yaxis.set_tick_params(labelsize=number_y_size)
ax.yaxis.set_ticks(np.arange(0, 2, 1))   
ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
plt.grid()
fig.savefig(f"../outputs/vague.svg", format='svg')
plt.show()



        