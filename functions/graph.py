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



def graph(X,W,H,D,t,speed,graph_para,r,s,L):
    
        wild, heterozygous, drive, theorical_wave, origin, grid, semilogy, ylim, xlim, mod, graph_type, save_figure, speed_proportion = graph_para

        fig, ax = plt.subplots()
        
        # Plot evolution for wild, heterozygous, drive (nb of individuals or proportions)
        if graph_type == "Individuals" : Y = [W,H,D]
        if graph_type == "Proportions" : Y = [W/(W+H+D),H/(W+H+D),D/(W+H+D)]
        # what to plot
        plot = [wild, heterozygous, drive]
        # color for each
        col = ['green','gold','deeppink']
        # label for each
        lab = ['Wild-type','Heterozygous','Drive']    
        # plot considering a log y-scale or not
        for i in range(3) :
            if plot[i] :
                if semilogy : ax.semilogy(X, Y[i], color = col[i], label=lab[i], linewidth = 2)
                else : ax.plot(X, Y[i], color = col[i], label=lab[i], linewidth = 2)
        
        
        # Plot the analytical solutions
        if speed != None and theorical_wave :   
           
            A_sty = (-speed+np.sqrt(speed**2 + 4))/2
            B_sty = (-speed+np.sqrt(speed**2 - 4*((1-s)*(r+1)-1)))/2                       
            b_1 = 1/3; b_2 = 1/3 
            
            if r == 0 :
                C_sty_0 = -speed/2
                D_sty_0 = (-speed-np.sqrt(speed**2 - 4*(1-2*s)))/2  
            
            if r != 0 :
                alpha = (-(1-2*s)*(1-r) + np.sqrt((1-2*s)**2 * (1-r)**2 + 8 * r**2 * (1-s)))/(4 * r * (1-s))
                C_sty = (-speed-np.sqrt(speed**2 + 4*r*(2*(1-s)*alpha+1)))/2
                D_sty = (-speed-np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2  
                D_sty_plus = (-speed+np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2   
            
                A_scr = (2*(1-s)*r)/(2*(1-s)*(1-r+2*alpha*r)-(1-r)) 
                B_scr = 1 + (-2*alpha*(1-s)*r)/(2*(1-s)*(1-r+2*alpha*r)-(1-r)) 
                
                a_3 =  (alpha-2)/3; a_4 =  (1/3)*(1 + A_scr*(2-alpha))    
                
                
            abscisse_1 = X[0:origin*len(X)//L]     
            abscisse_2 = X[origin*len(X)//L:len(X)] 
            
            if wild == True :
                if semilogy == False : ax.plot(abscisse_1, b_1*np.exp((abscisse_1-origin)*A_sty), color = 'blue', label='NO2', linestyle = '--')
                else : ax.semilogy(abscisse_1, b_1*np.exp((abscisse_1-origin)*A_sty), color = 'blue', label='NO2', linestyle = '--')
                if r != 0 : 
                    if semilogy == False : ax.plot(abscisse_2, 1-alpha*a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*B_scr*np.exp((abscisse_2-origin)*C_sty), color = 'yellowgreen', label='NO7 racine -', linestyle = '--')
                    else : ax.semilogy(abscisse_2, 1-alpha*a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*B_scr*np.exp((abscisse_2-origin)*C_sty), color = 'yellowgreen', label='NO7 racine -', linestyle = '--')
                    if semilogy == False : ax.plot(abscisse_2, 1-alpha*a_4*np.exp((abscisse_2-origin)*D_sty_plus) + a_3*B_scr*np.exp((abscisse_2-origin)*C_sty), color = 'deepskyblue', label='NO7 racine +', linestyle = '--')
                    else : ax.semilogy(abscisse_2, 1-alpha*a_4*np.exp((abscisse_2-origin)*D_sty_plus) + a_3*B_scr*np.exp((abscisse_2-origin)*C_sty), color = 'deepskyblue', label='NO7 racine +', linestyle = '--')
                    #if semilogy == False : ax.plot(abscisse_2, 1+a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*D_scr*np.exp((abscisse_2-origin)*C_sty), color = 'deepskyblue', label='NO7', linestyle = '--')
                    #else : ax.semilogy(abscisse_2, 1+a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*D_scr*np.exp((abscisse_2-origin)*C_sty), color = 'deepskyblue', label='NO7', linestyle = '--')
                if r == 0 : 
                    if semilogy == False : ax.plot(abscisse_2, 1-(2/3)*np.exp((abscisse_2-origin)*C_sty_0), color = 'deepskyblue', label='NO7', linestyle = '--')
                    else : ax.semilogy(abscisse_2, 1-(2/3)*np.exp((abscisse_2-origin)*C_sty_0), color = 'deepskyblue', label='NO7', linestyle = '--')
            
            if drive == True :
                if semilogy == False : ax.plot(abscisse_1, b_2*np.exp((abscisse_1-origin)*B_sty), color = 'orange', label='ND2', linestyle = '--')
                else : ax.semilogy(abscisse_1, b_2*np.exp((abscisse_1-origin)*B_sty), color = 'orange', label='ND2', linestyle = '--')
                if r != 0 : 
                    if semilogy == False : ax.plot(abscisse_2, a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*A_scr*np.exp((abscisse_2-origin)*C_sty), color = 'gold', label='ND7 racine -', linestyle = '--')
                    else :ax.semilogy(abscisse_2, a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*A_scr*np.exp((abscisse_2-origin)*C_sty), color = 'gold', label='ND7 racine -', linestyle = '--')
                    if semilogy == False : ax.plot(abscisse_2, a_4*np.exp((abscisse_2-origin)*D_sty_plus) + a_3*A_scr*np.exp((abscisse_2-origin)*C_sty), color = 'pink', label='ND7 racine +', linestyle = '--')
                    else :ax.semilogy(abscisse_2, a_4*np.exp((abscisse_2-origin)*D_sty_plus) + a_3*A_scr*np.exp((abscisse_2-origin)*C_sty), color = 'pink', label='ND7 racine +', linestyle = '--')
                    #if semilogy == False : ax.plot(abscisse_2, a_3*A_scr*np.exp((abscisse_2-origin)*C_sty), color = 'blue', label='ND7 sol part', linestyle = '--')
                    #else :ax.semilogy(abscisse_2, a_3*A_scr*np.exp((abscisse_2-origin)*C_sty), color = 'blue', label='ND7 sol part', linestyle = '--')
                    #if semilogy == False : ax.plot(abscisse_2, (1/alpha)*(-a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*(1-D_scr)*np.exp((abscisse_2-origin)*C_sty)), color = 'deepskyblue', label='ND7', linestyle = '--')
                    #else :ax.semilogy(abscisse_2, (1/alpha)*(-a_4*np.exp((abscisse_2-origin)*D_sty) + a_3*(1-D_scr)*np.exp((abscisse_2-origin)*C_sty)), color = 'deepskyblue', label='ND7', linestyle = '--')
                    
                if r == 0 : 
                    if semilogy == False : ax.plot(abscisse_2, (1/3)*np.exp((abscisse_2-origin)*D_sty_0), color = 'deepskyblue', label='ND7', linestyle = '--')
                    else :ax.semilogy(abscisse_2, (1/3)*np.exp((abscisse_2-origin)*D_sty_0), color = 'deepskyblue', label='ND7', linestyle = '--')
         
        if semilogy == True and ylim == False : defaultylim = (0.00001,1.1)
        if semilogy == False and ylim == False : defaultylim = (-0.03,1.03)  
                
        if ylim == False and xlim == False: 
            ax.set(xlabel='Space', ylabel=graph_type, ylim=defaultylim)   
        if ylim == False and xlim != False : 
            ax.set(xlabel='Space', ylabel=graph_type, xlim = xlim, ylim=defaultylim)
        if ylim != False and xlim == False: 
            ax.set(xlabel='Space', ylabel=graph_type, ylim = ylim)
        if ylim != False and xlim != False: 
            ax.set(xlabel='Space', ylabel=graph_type, xlim = xlim, ylim = ylim, title=f'r={np.round(r,4)}, s={np.round(s,4)}, t={np.round(t,2)}')
        if grid == True : 
            ax.grid()
        ax.xaxis.label.set_size(20)
        ax.yaxis.label.set_size(20)
        ax.yaxis.set_ticks(np.arange(0, 2, 1))
        ax.yaxis.set_tick_params(labelsize=20)
        # Pour enlever les chiffres des axes :
        #ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        #ax.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
        
        ax.legend()
        
        
        # Saving figures
        if save_figure :            
            if t == 0 : 
                actual_dir = os.getcwd()
                print ("The current working directory is %s" % actual_dir)
                new_dir = f"../outputs/evo_r_{r}_s_{s}_semilogy_{semilogy}"
                try:
                    os.mkdir(new_dir)
                except OSError:
                    print ("Creation of the directory %s failed" % new_dir)
                else:
                    print ("Successfully created the directory %s " % new_dir)
                    
            fig.savefig(f"../outputs/evo_r_{r}_s_{s}_semilogy_{semilogy}/t_{t}.pdf")   
            #columns = [X,W,D]; np.savetxt(f"t_{t}.txt", np.column_stack(columns), fmt='%.3e', delimiter="  ")  
            
        plt.show() 
        
        if speed != None :
            if r == 0 : 
                if speed**2 - 4*(1-2*s) < 0 :
                    print("\n Attention ! La condition 2 n'est pas vérifiée (ni 2a, ni 2b).")
                if speed**2 - 4*(1-2*s) < 0.00001 and  speed**2 - 4*(1-2*s) > -0.00001 :
                    print("\n Attention ! On est (à peu près) dans le cas de la condition 2b.")
            if r != 0 : 
                if speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1) < 0 :
                    print("\n Attention ! La condition 2 n'est pas vérifiée (ni 2a, ni 2b).")
                if speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1) < 0.00001 and  speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1) > -0.00001 :
                    print("\n Attention ! On est (à peu près) dans le cas de la condition 2b.")
                    
     
        
