#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:29:46 2021

@author: lena
"""

# Libraries
import matplotlib.pyplot as plt
import os
import numpy as np


# Graph parameters 
title_size = 15
label_size = 17
legend_size = 12
line_size = 10
number_on_x_axe = False
number_x_size = 10
number_y_size = 20

# Saving figures
save_column = True


def graph(X,W,H,D,t,graph_para,bio_para,num_para, model_para):
       
        r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
        T,L,M,N,theta = num_para[0:5]
        wild, heterozygous, drive, mod, grid, semilogy, xlim, graph_type, show_graph_ini, show_graph_end, save_fig = graph_para

        fig, ax = plt.subplots()
        
        # Plot evolution for wild, heterozygous, drive (nb of individuals or proportions)
        if graph_type == "Genotype densities" : 
            nb = 3; Y = [W, D, H]
        if graph_type == "Genotype proportions" : 
            nb = 3; Y = [W/(W+H+D), D/(W+H+D), H/(W+H+D)]            
        if graph_type == "Allele densities" : 
            nb = 2
            if conversion_timing == "zygote" : 
                Y = [W+0.5*H, D+0.5*H]
            if conversion_timing == "germline" : 
                Y = [W+0.5*(1-c)*H, D+0.5*(1+c)*H]
        if graph_type == "Allele proportions" : 
            nb = 2
            if conversion_timing == "zygote" : 
                Y = [(W+0.5*H)/(W+H+D), (D+0.5*H)/(W+H+D) ]
            if conversion_timing == "germline" : 
                Y = [(W+0.5*(1-c)*H)/(W+H+D), (D+0.5*(1+c)*H)/(W+H+D) ]
                
        # what to plot
        plot = [wild, drive, heterozygous]
        # color for each
        col = ['cornflowerblue','crimson','orange']
        # label for each
        lab = ['Wild-type','Drive','Heterozygous']    
        # plot considering a log y-scale or not
        for i in range(nb) :
            if plot[i] :
                if semilogy : ax.semilogy(X, Y[i], color = col[i], label=lab[i], linewidth = line_size)
                else : ax.plot(X, Y[i], color = col[i], label=lab[i], linewidth = line_size)
        #ax.plot(X, np.ones(len(X))*threshold)        
         
        # Graphic size, title and labels
        if semilogy : defaultylim = (0.00001,1.1)
        else : defaultylim = (-0.03,1.03)  
        if xlim == None : 
            ax.set(xlabel='Space', ylabel=graph_type, ylim = defaultylim)
            ax.set_title(f"{graph_type} t = {t}", fontsize = title_size, loc='right')
        else : 
            ax.set(xlabel='Space', xlim = xlim, ylabel=graph_type, ylim = defaultylim)
            ax.set_title(f"{graph_type} t = {t}", fontsize = title_size, loc='right')
            
        # Grid
        if grid == True : 
            ax.grid()
            
        # Labels and legend sizes
        ax.xaxis.label.set_size(label_size)
        ax.yaxis.label.set_size(label_size)     
        if number_on_x_axe :
            ax.xaxis.set_tick_params(labelsize=number_x_size)
        else : 
            ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
            ax.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        ax.yaxis.set_tick_params(labelsize=number_y_size)
        ax.yaxis.set_ticks(np.arange(0, 2, 1))
        plt.rc('legend', fontsize=legend_size)  
        #ax.legend(bbox_to_anchor=(1.02,1.15), ncol=2)
        ax.legend(bbox_to_anchor=(0.553,1.13), ncol=2)
        
        # Show the graph      
        plt.show()
        
        # Saving figures and datas
        if save_fig : 
            directory = f"evolution/{conversion_timing}_r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"
            save_fig_or_data(directory, fig, [], f"t_{t}", bio_para, num_para, model_para)
            #num = str(int(t)//mod)
            #if len(num)==1: num = '0'+'0'+num
            #if len(num)==2: num = '0'+num
            #save_fig_or_data(directory, fig, [], f"{num}", bio_para, num_para, model_para)
            #columns = [X,W,D]; np.savetxt(f"../outputs/{directory}/t_{t}.txt", np.column_stack(columns), fmt='%.3e', delimiter="  ") 
        
    
    
  
def graph_2D(t, W, H, D, N, graph_para, bio_para, num_para, model_para):
    # Parameters
    wild, heterozygous, drive = graph_para[0:3]; save_fig = graph_para[-1]
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    # Figure
    for i in range(3):
        if [wild, heterozygous, drive][i] :
            genotype = ["WW","DW","DD"][i]
            heatmap_values = np.resize([W,H,D][i],(N+1,N+1)).transpose()
            fig, ax = plt.subplots() 
            im = ax.imshow(heatmap_values,cmap='Blues', aspect='auto', vmin=0, vmax=1)  
            ax.figure.colorbar(im, ax=ax)   
            fig.suptitle(f"Genotype {genotype} at time {np.round(t,2)}", fontsize=14)
            if save_fig :
                directory = f"evolution_2D/{conversion_timing}_r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"                       
                save_fig_or_data(directory, fig, [], f"{genotype}_t_{t}", bio_para, num_para, model_para)
            plt.show() 
    
def graph_2D_contour(t, W, H, D, N, Z_list, nb_graph, graph_para, bio_para, num_para, model_para):
    # Parameters
    wild, heterozygous, drive, mod = graph_para[0:4]; save_fig = graph_para[-1]
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    # Figure
    contour_threshold = 0.2
    for i in range(3):
        if [wild, heterozygous, drive][i] :
            genotype = ["WW","DW","DD"][i]
            heatmap_values = np.resize([W,H,D][i],(N+1,N+1)).transpose()
            fig, ax = plt.subplots() 
            g1 = lambda x,y: heatmap_values[int(y),int(x)] 
            g2 = np.vectorize(g1)
            x = np.linspace(0,heatmap_values.shape[1], 1001)[:-1]
            y = np.linspace(0,heatmap_values.shape[0], 1001)[:-1]
            X, Y= np.meshgrid(x,y)
            Z = g2(X,Y)  
            Z_list[nb_graph-1] = Z
            #x = np.linspace(0,heatmap_values.shape[1], heatmap_values.shape[1]*100)
            #y = np.linspace(0,heatmap_values.shape[0], heatmap_values.shape[0]*100)
            #X, Y= np.meshgrid(x[:-1],y[:-1])
            #Z = g2(X[:-1],Y[:-1])  
            #Z_list[nb_graph-1] = Z[:,1:]
            ax.set_aspect('equal', adjustable='box')
            #im = ax.imshow(heatmap_values,cmap='Blues', aspect='auto', vmin=0, vmax=1)  
            #ax.figure.colorbar(im, ax=ax)   
            #ax.contour(np.arange(N+1), np.arange(N+1), heatmap_values, levels=[0.4]) 
            for i in range(nb_graph) : 
                label = f'{int(mod*i)}'
                if i == nb_graph - 1 : label = f'{int(t)}'
                if np.any(Z_list[i] > contour_threshold) and np.any(Z_list[i] < contour_threshold): contour = ax.contour(Z_list[i], [contour_threshold], linewidths=3, extent=[0-0.5, x[:-1].max()-0.5,0-0.5, y[:-1].max()-0.5])
                fmt = {}; fmt[contour_threshold] = label
                ax.clabel(contour, np.ones(1)*contour_threshold, inline=True, fmt=fmt, fontsize=10) 
            if save_fig :
                directory = f"evolution_2D/{conversion_timing}_r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"       
                save_fig_or_data(directory, fig, [], f"contour_{genotype}_t_{t}", bio_para, num_para, model_para)
            plt.show()   
    return(Z_list)
             
    
    
    
    
    


# Create a new directory at the localisation f"../outputs{path}"
def create_directory(path) :
    # Directory to create
    new_dir = f"../outputs{path}"
    # ... if it doesn't already exist
    if not os.path.exists(new_dir): 
        try:
            os.mkdir(new_dir)
        except OSError:
            print ("Fail : %s " % new_dir)

   
    
def create_para_txt(path, bio_para, num_para, model_para):
    T,L,M,N,theta = num_para[0:5]
    file = open(f"../outputs/{path}/0_parameters.txt", "w") 
    file.write(f"Parameters : \nT = {T} \nL = {L} \nM = {M} \nN = {N} \ntheta = {theta}")  
    if bio_para != None :
        r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
        file.write(f"\nr = {r} \ns = {s} \nh = {h} \nc = {c} \nconversion_timing = {conversion_timing} \na = {a} \ndifWW = {difWW} \ndifDW = {difDW} \ndifDD = {difDD}")                     
    if model_para != None : 
        CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
        file.write(f"\ngrowth_dynamic = {growth_dynamic} \ndeath_dynamic = {death_dynamic}")  
    file.close()
    

# Create the tree of directories at the localisation f"../outputs{path}" with parameters.txt
def create_path(directory):
    #actual_dir = os.getcwd()
    #print ("\nWorking directory : %s" % actual_dir) 
    lst_slash_index = [-1]
    path = ""          
    for pos,char in enumerate(directory):
        if(char == "/"):
            lst_slash_index.append(pos)
            path = path+"/"+directory[lst_slash_index[-2]+1:lst_slash_index[-1]]
            create_directory(path)
    path = path+"/"+directory[lst_slash_index[-1]+1:len(directory)]    
    create_directory(path)
                 
    
# Save figures and datas regarding what we are doing. 
# fig is the figure to save (title.png and title.pdf)
# data is the data to save (title.txt)
def save_fig_or_data(directory, fig, data, title, bio_para, num_para, model_para):
    # Create the tree of directories in outputs (if it does not exist) with parameters.txt
    create_path(directory)
    if not os.path.exists(f"../outputs/{directory}/0_parameters.txt") :
        create_para_txt(directory, bio_para, num_para, model_para)
    # Save figure
    if fig != [] :
        fig.savefig(f"../outputs/{directory}/{title}.png", format='png')
        fig.savefig(f"../outputs/{directory}/{title}.svg", format='svg')
        #fig.savefig(f"../outputs/{directories}/{title}.pdf", format='pdf')       
    # Save datas
    if data != [] :
        np.savetxt(f"../outputs/{directory}/{title}.txt", data)   
              