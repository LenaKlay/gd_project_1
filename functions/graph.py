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
legend_size = 15
label_size = 20
number_on_x_axe = False
number_x_size = 10
number_y_size = 20
line_size = 4

# Saving figures
save_column = True


def graph(X,W,H,D,t,graph_para,bio_para,num_para,directory,file,title_fig):
    
        graph_type, wild, heterozygous, drive, grid, semilogy, xlim, save_fig, speed_proportion, show_graph_ini = graph_para

        fig, ax = plt.subplots()
        
        # Plot evolution for wild, heterozygous, drive (nb of individuals or proportions)
        if graph_type == "Population density" : Y = [W,H,D]
        if graph_type == "Proportions" : Y = [W/(W+H+D),H/(W+H+D),D/(W+H+D)]
        # what to plot
        plot = [wild, heterozygous, drive]
        # color for each
        col = ['green','orange','deeppink']
        # label for each
        lab = ['Wild-type','Heterozygous','Drive']    
        # plot considering a log y-scale or not
        for i in range(3) :
            if plot[i] :
                if semilogy : ax.semilogy(X, Y[i], color = col[i], label=lab[i], linewidth = line_size)
                else : ax.plot(X, Y[i], color = col[i], label=lab[i], linewidth = line_size)

        # Graphic size, title and labels
        if semilogy : defaultylim = (0.00001,1.1)
        else : defaultylim = (-0.03,1.03)  
        if xlim == None : 
            ax.set(xlabel='Space', ylabel=graph_type, ylim = defaultylim, title=title_fig)   
        else : 
            ax.set(xlabel='Space', xlim = xlim, ylabel=graph_type, ylim = defaultylim, title=title_fig)
            
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
        ax.legend()
        
        # Show the graph      
        plt.show()
        
        # Saving figures and data
        if save_fig : 
            save_figure(t, fig, directory, file, bio_para, num_para) 
            if save_column : 
                columns = [X,W,D]; np.savetxt(f"../outputs/{directory}/columns_t_{t}.txt", np.column_stack(columns), fmt='%.3e', delimiter="  ") 
        
        
 





# Create a new directory "path" in outputs         
def create_directory(path, bio_para, num_para, parameters_txt) :
    # Directory to create
    new_dir = f"../outputs{path}"
    # ... if it doesn't already exist
    if not os.path.exists(new_dir): 
        try:
            os.mkdir(new_dir)
        except OSError:
            print ("Fail : %s " % new_dir)
        else:
            print ("Success : %s " % new_dir)
                        
            # Write parameters in the new directory
            if parameters_txt : 
                T,L,M,N,mod,theta = num_para
                file = open(f"../outputs{path}/0_parameters.txt", "w") 
                if bio_para != None :   
                    r,s,h,a,difW,difH,difD,c,homing = bio_para
                    file.write(f"Parameters : \nr = {r} \ns = {s} \nh = {h} \nc = {c} \nhoming = {homing} \nT = {T} \nL = {L} \nM = {M} \nN = {N} \ntheta = {theta} \nmod = {mod}") 
                else :   #  no bio_para for tanaka
                    file.write(f"Parameters : \nT = {T} \nL = {L} \nM = {M} \nN = {N} \ntheta = {theta} \nmod = {mod}")  
                file.close()
    #else : 
    #    print("Already exists : %s" % new_dir)

                 
    
# Create the tree of directories in outputs and save pdf figures
def save_figure(t, fig, directories, pdf_title, bio_para, num_para) : 
   
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
      
    # Save figures .pdf in the new directory
    if t != None : 
        fig.savefig(f"../outputs{path}/t_{t}_{pdf_title}.pdf")  
        fig.savefig(f"../outputs{path}/t_{t}_{pdf_title}.png") 
    else : 
        fig.savefig(f"../outputs{path}/{pdf_title}.pdf") 
        fig.savefig(f"../outputs{path}/{pdf_title}.png")
        
    
              
                