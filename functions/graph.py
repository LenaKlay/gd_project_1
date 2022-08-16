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


def graph(X,W,H,D,t,graph_para,bio_para,num_para):
       
        r,s,h,a,difW,difH,difD,c,homing = bio_para
        T,L,M,N,mod,theta = num_para
        graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion, show_graph_ini, show_graph_fin = graph_para

        fig, ax = plt.subplots()
        
        # Plot evolution for wild, heterozygous, drive (nb of individuals or proportions)
        if graph_type == "Population densities" : 
            nb = 3; Y = [W, D, H]
        if graph_type == "Population proportions" : 
            nb = 3; Y = [W/(W+H+D), D/(W+H+D), H/(W+H+D)]            
        if graph_type == "Allele densities" : 
            nb = 2
            if homing == "zygote" : 
                Y = [W+0.5*H, D+0.5*H]
            if homing == "germline" : 
                Y = [W+0.5*(1-c)*H, D+0.5*(1+c)*H]
        if graph_type == "Allele proportions" : 
            nb = 2
            if homing == "zygote" : 
                Y = [(W+0.5*H)/(W+H+D), (D+0.5*H)/(W+H+D) ]
            if homing == "germline" : 
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
            #ax.set_title(f"t = {t}", fontsize = title_size, loc='right')
        else : 
            ax.set(xlabel='Space', xlim = xlim, ylabel=graph_type, ylim = defaultylim)
            ax.set_title(f"t = {t}", fontsize = title_size, loc='right')
            
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
            directory = f"evolution/{homing}/r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"
            #save_fig_or_data(directory, fig, [], f"t_{t}", bio_para, num_para)
            num = str(int(t)//mod)
            if len(num)==1: num = '0'+'0'+num
            if len(num)==2: num = '0'+num
            save_fig_or_data(directory, fig, [], f"{num}", bio_para, num_para)
            #columns = [X,W,D]; np.savetxt(f"../outputs/{directory}/t_{t}.txt", np.column_stack(columns), fmt='%.3e', delimiter="  ") 
        
    


# Create a new directory at the localisation f"../outputs{path}"
def create_directory(path, bio_para, num_para, parameters_txt) :
    # Directory to create
    new_dir = f"../outputs{path}"
    # ... if it doesn't already exist
    if not os.path.exists(new_dir): 
        try:
            os.mkdir(new_dir)
        except OSError:
            print ("Fail : %s " % new_dir)
        #else:
            #print ("Success : %s " % new_dir)
                        
            # Write parameters in the new directory
            if parameters_txt : 
                T,L,M,N,mod,theta = num_para
                file = open(f"../outputs{path}/0_parameters.txt", "w") 
                
                if bio_para == None :  #  no bio_para for tanaka
                    file.write(f"Parameters : \nT = {T} \nL = {L} \nM = {M} \nN = {N} \ntheta = {theta} \nmod = {mod}")  
                else : 
                    r,s,h,a,difW,difH,difD,c,homing = bio_para
                    file.write(f"Parameters : \nr = {r} \ns = {s} \nh = {h} \nc = {c} \nhoming = {homing} \nT = {T} \nL = {L} \nM = {M} \nN = {N} \ntheta = {theta} \nmod = {mod}")                     
                file.close() 
    #else : 
    #    print("Already exists : %s" % new_dir)
    
    

# Create the tree of directories at the localisation f"../outputs{path}" with parameters.txt
def create_path(directories, bio_para, num_para):
    #actual_dir = os.getcwd()
    #print ("\nWorking directory : %s" % actual_dir) 
    lst_slash_index = [-1]
    path = ""          
    for pos,char in enumerate(directories):
        if(char == "/"):
            lst_slash_index.append(pos)
            path = path+"/"+directories[lst_slash_index[-2]+1:lst_slash_index[-1]]
            create_directory(path, bio_para, num_para, False)
    path = path+"/"+directories[lst_slash_index[-1]+1:len(directories)]    
    create_directory(path, bio_para, num_para, True)
                 
    
# Save figures and datas regarding what we are doing. 
# fig is the figure to save (title.png and title.pdf)
# data is the data to save (title.txt)
def save_fig_or_data(directories, fig, data, title, bio_para, num_para):
    # Create the tree of directories in outputs (if it does not exist) with parameters.txt
    create_path(directories, bio_para, num_para)
    # Save figure
    if fig != [] :
        fig.savefig(f"../outputs/{directories}/{title}.png", format='png')
        fig.savefig(f"../outputs/{directories}/{title}.pdf", format='pdf') #; fig.savefig(f"../outputs/{directories}/{title}.png") 
        fig.savefig(f"../outputs/{directories}/{title}.svg", format='svg')
    # Save datas
    if data != [] :
        np.savetxt(f"../outputs/{directories}/{title}.txt", data)   
              