#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:29:46 2021

@author: lena
"""

def graph(X,W,H,D,t,graph_type,speed,origin):
        fig, ax = plt.subplots()
        if graph_type == "Individuals" :  
            if wild == True :
                if semilogy == False : ax.plot(X, W, color = 'green', label='Wild-type', linewidth = 4)
                else : ax.semilogy(X, W, color = 'green', label='Wild-type')
            if heterozygous == True :
                if semilogy == False : ax.plot(X, H, color = 'gold', label='Heterozygous', linewidth = 4)
                else : ax.semilogy(X, H, color = 'gold', label='Heterozygous')         
            if drive == True :
                if semilogy == False : ax.plot(X, D, color = 'orange', label='Drive', linewidth = 4)
                else : ax.semilogy(X, D, color = 'orange', label='Drive') 
            if N_alpha == True and r!=0 :
                alpha = (-(1-2*s)*(1-r) + np.sqrt((1-2*s)**2 * (1-r)**2 + 8 * r**2 * (1-s)))/(4 * r * (1-s))
                if semilogy == False : ax.plot(X, alpha*D+W, color = 'pink', label='N_alpha')
                else : ax.semilogy(X, alpha*D+W, color = 'pink', label='N_alpha') 
        if graph_type == "Proportions" :  
            if wild == True :
                ax.plot(X, W/(W+H+D), color = 'green', label='Wild-type', linewidth = 4)
            if heterozygous == True :
                ax.plot(X, H/(W+H+D), color = 'gold', label='Heterozygous', linewidth = 4)
            if drive == True :
                ax.plot(X, D/(W+H+D), color = 'orange', label='Drive', linewidth = 4) 
        
        if speed != None :   
            
            A_sty = (-speed+np.sqrt(speed**2 + 4))/2
            B_sty = (-speed+np.sqrt(speed**2 - 4*((1-s)*(r+1)-1)))/2            
            
            b_1 = 1/3
            b_2 = 1/3 
            
            if r == 0 :
                C_sty_0 = -speed/2
                D_sty_0 = (-speed-np.sqrt(speed**2 - 4*(1-2*s)))/2  
            
            if r != 0 :
                alpha = (-(1-2*s)*(1-r) + np.sqrt((1-2*s)**2 * (1-r)**2 + 8 * r**2 * (1-s)))/(4 * r * (1-s))
                C_sty = (-speed-np.sqrt(speed**2 + 4*r*(2*(1-s)*alpha+1)))/2
                D_sty = (-speed-np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2     #  -0.03  # -0.38  
                D_sty_plus = (-speed+np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2   
                
                #if semilogy == False : 
                #    print('s = ',s)
                #    print('r = ',r, '\n')
                #    print('C_sty =', (-speed-np.sqrt(speed**2 + 4*r*(2*(1-s)*alpha+1)))/2)
                #    print('D_sty+ =', (-speed+np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2)
                #    print('D_sty- =', (-speed-np.sqrt(speed**2 - 4*(2*(1-s)*(1-r+alpha*r)-1)))/2)
            
                A_scr = (2*(1-s)*r)/(2*(1-s)*(1-r+2*alpha*r)-(1-r)) 
                B_scr = 1 + (-2*alpha*(1-s)*r)/(2*(1-s)*(1-r+2*alpha*r)-(1-r)) 
                #D_scr = r/(1-2*s*(1-r)*alpha)
                
                a_3 =  (alpha-2)/3 
                a_4 =  (1/3)*(1 + A_scr*(2-alpha))      #   (1/3)*(-2+D_scr*(2-alpha)) 
                
                if semilogy == "True" :
                    print("speed = ", speed)
                    print("A_sty =", A_sty)
                    print("B_sty =", B_sty)
                    print("alpha =", alpha)
                    print("C_sty =", C_sty)
                    print("D_sty =", D_sty)
                    print("D_sty_plus =", D_sty_plus)
                
            abscisse_1 = X[0:origin*len(X)//L]      # X[0:861*len(X)//L+1]   
            abscisse_2 = X[origin*len(X)//L:len(X)]   #  X[867*len(X)//L:len(X)]  
           
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
             
            if N_alpha == True and r != 0 :
                if semilogy == False : ax.plot(abscisse_2, a_3*np.exp((abscisse_2-origin)*C_sty)+1, color = 'blue', label='Nalpha7')
                else :ax.semilogy(abscisse_2, a_3*np.exp((abscisse_2-origin)*C_sty)+1, color = 'blue', label='Nalpha7')
                if semilogy == False : ax.plot(abscisse_1, (b_1*np.exp((abscisse_1-origin)*A_sty))+alpha*(b_2*np.exp((abscisse_1-origin)*B_sty)), color = 'orange', label='Nalpha2')
                else : ax.semilogy(abscisse_1, (b_1*np.exp((abscisse_1-origin)*A_sty))+alpha*(b_2*np.exp((abscisse_1-origin)*B_sty)), color = 'orange', label='Nalpha2')
                    
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
        #fig.savefig(f"r_{r}_s_{s}_T_{T}_L_{L}_M_Tfois{int(M/T)}_N_Lfois{int(N/L)}_semilogy_{semilogy}.pdf")
        fig.savefig(f"r_{r}_s_{s}_t_{t}_semilogy_{semilogy}.pdf")   
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
                    
     
        
