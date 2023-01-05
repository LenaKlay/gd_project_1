#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:03:42 2021

@author: lena
"""

# Libraries
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg  as la
import matplotlib.pyplot as plt

# External functions 
from graph import graph, graph_2D, graph_2D_contour, save_fig_or_data


# equ is the allele concerned by the equation.
# oth1 and oth2 are the two other alleles, indifferently.

# Growth function
def growth(equ, oth1, oth2, bio_para, model_para):
    # Parameters
    r,s,h,a,difW,difH,difD,c,conversion_timing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    # Function
    n = equ+oth1+oth2
    if growth_dynamic == "constant":
        return(r+1)
    if growth_dynamic == "allee_effect" :
        #return(r*(1-n)*(n-a)+1)
        return(np.maximum(r*(1-n)*(n-a)+1, np.zeros(len(n))))
    if growth_dynamic == "logistical" and linear_growth == False :
        return(r*(1-n)+1)
    if growth_dynamic == "logistical" and linear_growth == True :  
        return(r*np.minimum(1-(equ+oth1+oth2),equ) + equ)        
       


# Death function    
def death(equ, oth1, oth2, bio_para, model_para):
    # Parameters
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    # Function
    n = (equ+oth1+oth2)
    if death_dynamic == "constant" :
        return(1)
    if death_dynamic == "logistical" :
        return(1+r*n)
    if death_dynamic == "allee_effect" :
        return(r+1+r*(1-n)*(n-a))
    
    
# Mating function
def mating(W, H, D, bio_para, model_para):  
    # Parameters
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    
    # Function
    if not linear_mating : 
        
        WW = (W**2)/(W+H+D)
        HH = (H**2)/(W+H+D)
        DD = (D**2)/(W+H+D)    
        WH = (W*H)/(W+H+D)   # *2
        WD = (W*D)/(W+H+D)   # *2 
        HD = (D*H)/(W+H+D)   # *2
    
        if conversion_timing == "zygote" :       
            mat1 = WW + WH + 0.25*HH 
            mat2 = (1-c)*(WH + 2*WD + 0.5*HH + HD) 
            mat3 = c*WH + 2*c*WD + (0.5*c+0.25)*HH + (1+c)*HD + DD  
            
        if conversion_timing == "germline":
            mat1 = WW + (1-c)*WH + 0.25*(1-c)**2*HH 
            mat2 = (1+c)*WH + 2*WD + 0.5*(1-c**2)*HH + (1-c)*HD
            mat3 = 0.25*(c+1)**2*HH + (1+c)*HD + DD
            
    if linear_mating and conversion_timing == "zygote" : 
        mat1 = (W>D)
        mat2 = 0
        mat3 = ((W>D)+1)  

    return([mat1,mat2,mat3])   



def s1_s2(conversion_timing, c, h, s, graph_type) :
    
    threshold = 0.2 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    WT_allele_wave = False
    p_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    
    # Determine values of s1 and s2, and in which regime we are.
    s_1 = c/(1-h*(1-c))   
    if conversion_timing == "zygote" :
        s_2 = c/(2*c + h*(1-c))
        if s_1 < s_2 : 
            if graph_type != None: print("\nCoexistence", ":", np.round(s_1,3), "< s <", np.round(s_2,3))
            p_star = (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h)))
            if 0<p_star and p_star<1 : 
                WT_allele_wave = True; threshold = ((1-p_star)+1)/2 
        else :  
            if graph_type != None: print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        lin = c*(1-2*s)-(1-c)*s*h
        if lin > 0 : print("Linear speed :", 2*np.sqrt(lin))               
    if conversion_timing == "germline" : 
        s_2 = c/(2*c*h + h*(1-c))
        if s_1 < s_2 :
            if graph_type != None: print("\nCoexistence", ":", np.round(s_1,3) , "< s <", np.round(s_2,3))
            p_star = ((1-s*h)*(1+c)-1)/(s*(1-2*h))
            if 0<p_star and p_star<1 : 
                WT_allele_wave = True; threshold = ((1-p_star)+1)/2 
        else :  
            if graph_type != None: 
                print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        lin = c*(1-2*s*h)-(1-c)*s*h > 0 
        if lin > 0 : print("Linear speed :", 2*np.sqrt(lin))
        
    return(threshold, p_star, WT_allele_wave, s_1, s_2)
    
    

# Determine values of s1 and s2, and in which regime we are.   
def s1_s2_num(conversion_timing, c, h, s, graph_type) :
    s_1 = c/(1-h*(1-c))   
    if conversion_timing == "zygote" :
        s_2 = c/(2*c + h*(1-c))
        if s_1 < s_2 : 
            if graph_type != None: print("\nCoexistence", ":", np.round(s_1,3), "< s <", np.round(s_2,3))
            #p_star = (s*(1-(1-c)*(1-h)) - c*(1-s))/(s*(1-2*(1-c)*(1-h)))
        else :  
            if graph_type != None: print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        lin = c*(1-2*s)-(1-c)*s*h
        if lin > 0 : print("Linear speed :", 2*np.sqrt(lin))               
    if conversion_timing == "germline" : 
        s_2 = c/(2*c*h + h*(1-c))
        if s_1 < s_2 :
            if graph_type != None: print("\nCoexistence", ":", np.round(s_1,3) , "< s <", np.round(s_2,3))
            #p_star = ((1-s*h)*(1+c)-1)/(s*(1-2*h))
        else :  
            if graph_type != None: 
                print("\nBistability", ":", np.round(s_2,3), "< s <", np.round(s_1,3))
        lin = c*(1-2*s*h)-(1-c)*s*h > 0 
        if lin > 0 : print("Linear speed :", 2*np.sqrt(lin))
    return(s_1, s_2)
    

# Main evolution function (1D)     
def evolution(bio_para, model_para, num_para, graph_para) :  

    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    T,L,M,N,theta = num_para[0:5]
    wild, heterozygous, drive, mod, grid, semilogy, xlim, graph_type, show_graph_ini, show_graph_end, save_fig = graph_para
    
    # Parameters initialization
    position = np.array([])   # list containing the first position where the proportion of wild alleles is higher than the threshold value.
    speed_fct_of_time = np.array([])      # speed computed... 
    time = np.array([])       # ...for each value of time in this vector.
    threshold = 0.5 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    pw_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    coex_density = -1 
    
    # Regimes 
    if growth_dynamic == "logistical" and death_dynamic == "constant" : 
        s1, s2 = s1_s2_num(conversion_timing, c, h, s, graph_type)
    
    
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    X = np.linspace(0,N,N+1)*dx   
    
    # Initialization             
    if CI == "center" : 
        W = np.ones(N+1); W[0:N//2] = 0    # Wild individuals at t=0  
        H = np.zeros(N+1)                  # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//2] = 1   # Drive individuals at t=0
    if CI == "left" : 
        W = np.ones(N+1); W[0:N//10] = 0    # Wild individuals at t=0  
        H = np.zeros(N+1)                   # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//10] = 1   # Drive individuals at t=0
       
    if graph_type != None and show_graph_ini :
        graph(X,W,H,D,0,graph_para,bio_para,num_para, model_para)
    nb_graph = 1
        
    # Matrix
    C0 = -2*np.ones(N-1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N-1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N-1, N-1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    Bw = sp.identity(N-1)+((1-theta)*difWW*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    Bh = sp.identity(N-1)+((1-theta)*difDW*dt/dx**2)*A            # for W, H and D.
    Bd = sp.identity(N-1)+((1-theta)*difDD*dt/dx**2)*A     

    Bw_ = sp.identity(N-1)-(theta*difWW*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  
    Bh_ = sp.identity(N-1)-(theta*difDW*dt/dx**2)*A               # for W, H and D.
    Bd_ = sp.identity(N-1)-(theta*difDD*dt/dx**2)*A        

    # Evolution
    for t in np.linspace(dt,T,M) :
        
        t = round(t,2)
        
        f1 = growth(W,H,D,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[0] - death(W,H,D,bio_para,model_para)*W
        f2 = (1-s*h)*growth(H,W,D,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[1] - death(W,H,D,bio_para,model_para)*H
        f3 = (1-s)*growth(D,H,W,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[2] - death(W,H,D,bio_para,model_para)*D
      
        W[1:-1] = la.spsolve(Bw_, Bw.dot(W[1:-1]) + dt*f1[1:-1])
        H[1:-1] = la.spsolve(Bh_, Bh.dot(H[1:-1]) + dt*f2[1:-1])
        D[1:-1] = la.spsolve(Bd_, Bd.dot(D[1:-1]) + dt*f3[1:-1])
             
        W[0] = W[1]; H[0] = H[1]; D[0] = D[1]    # Neumann condition alpha=0
        W[-1] = W[-2]; H[-1] = H[-2]; D[-1] = D[-2]  # Neumann condition beta=0
        
        # Graph
        if t>mod*nb_graph and graph_type != None : 
            graph(X,W,H,D,t,graph_para,bio_para,num_para, model_para)
            nb_graph += 1
                                          

        # Compute the speed on WT allele density when coexistence density is close to zero
        # to avoid numerical errors (because then, proportions are very approximative), or when
        # there is no coexistence (coex_density = -1)         
        if coex_density < 0.1: 
            if conversion_timing == "zygote" : 
                WT = W+0.5*H
            if conversion_timing == "germline" : 
                WT = W+0.5*(1-c)*H
            p = WT/(W+H+D)
        # Else (coexistence state with final density > 0.1), compute the speed on WT allele density
        else : 
            if conversion_timing == "zygote" : 
                WT = (W+0.5*H)/(W+H+D) 
            if conversion_timing == "germline" : 
                WT = (W+0.5*(1-c)*H)/(W+H+D) 
            
        # Check if the WT wave is strictly increasing
        epsilon = -0.01 ; index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        if len(index_neg) > 4 : print("WT prop decreasing :", index_neg) 
               
        # Check if there is a coexistence state (if so adapt threshold and erase all speed values)
        # palliate centered
        center_index = np.arange(N//2-N//20,N//2+N//20); cond1 = np.all((p[center_index]>0.01) & (p[center_index]<0.99))
        # palliate left (if the wave move faster to the left)
        left_index = np.arange(N//2-N//10,N//2); cond2 = np.all((p[left_index]>0.01) & (p[left_index]<0.99))
        # palliate right (if the wave move faster to the left)
        right_index = np.arange(N//2,N//2+N//10); cond3 = np.all((p[right_index]>0.01) & (p[right_index]<0.99))
        if (cond1 or cond2 or cond3) and pw_star < 0 :    # if there exists a palliate in ]0,1[ around N//2 (and pw_star as not been changed yet)
            pw_star = np.mean(p[center_index])           
            coex_density = (W+H+D)[N//2]
            threshold = (pw_star+1)/2
            position = np.array([])   
            speed_fct_of_time = np.array([])   
            time = np.array([]) 
            print("coex ! p_w* =", pw_star)
            print("threshold =", threshold)          
        # we recorde the position only if the WT wave is still in the environment. We do not recorde the 0 position since the threshold value of the wave might be outside the window.
        if np.isin(True, WT>threshold) and np.isin(True, WT<0.99) and np.where(WT>threshold)[0][0] != 0:  
            # List containing the first position of WT where the proportion of wild alleles is higher than the threshold.  
            position = np.append(position, np.where(WT>threshold)[0][0])
            # Wait a minimum time and a minimum number of positions to compute the speed
            if len(position) > 20  : 
                # The speed is computed on the last fifth of the position vector.
                time = np.append(time,t)
                speed_fct_of_time = np.append(speed_fct_of_time, np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt)
                # if there is a big jump in the speed, stop the simulation. This sometimes happens in coexistance, when the palliate too close to 1.
                if len(speed_fct_of_time) > 100 : 
                    if abs(np.mean(speed_fct_of_time) - speed_fct_of_time[-1]) > 10: 
                        print("t =",t,"(coex : palliate very close to one)")
                        # erase the false value 
                        speed_fct_of_time = speed_fct_of_time[:-1]; time = time[:-1]
                        break
        # if the threshold value of the wave is outside the window, stop the simulation  
        if not(np.isin(False, WT>threshold) and np.isin(False, WT<threshold) ) :
            print("t =",t)
            break
        

    # last graph
    if show_graph_end :   
        graph(X,W,H,D,t,graph_para,bio_para,num_para, model_para)          
        
    # plot the speed function of time    
    if len(speed_fct_of_time) != 0 : 
        if graph_type!= None :
            fig, ax = plt.subplots()
            ax.plot(time, speed_fct_of_time) 
            ax.set(xlabel='Time', ylabel='Speed', title = f'Speed function of time')   
            if save_fig :
                directory = f"evolution/{conversion_timing}_r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"
                save_fig_or_data(directory, fig, speed_fct_of_time, "speed_fct_of_time", bio_para, num_para, model_para)
            plt.show() 
    else :
        print(f"No wave for r = {r} and s = {s}.") 

    return(W,H,D,time,speed_fct_of_time)



        
# Main evolution function      
def evolution_2D(bio_para, model_para, num_para, graph_para, CI_lenght) :  
    
    # Parameters loading
    r,s,h,a,difWW,difDW,difDD,c,conversion_timing = bio_para
    CI,growth_dynamic,death_dynamic,linear_growth,linear_mating = model_para
    T,L,M,N,theta = num_para[0:5]
    wild, heterozygous, drive, mod, grid, semilogy, xlim, graph_type, show_graph_ini, show_graph_end, save_fig = graph_para
    
    # Parameters initialization
    position = np.array([])   # list containing the first position where the proportion of wild alleles is higher than the threshold value.
    speed_fct_of_time = np.array([])      # speed computed... 
    time = np.array([])       # ...for each value of time in this vector.
    threshold = 0.5 # indicates which position of the wave we follow to compute the speed (first position where the WT wave come above the threshold)    
    pw_star = -1 # a value not admissible ; the value of the equilibrium will be change if there exists a coexistence stable state
    
    # Regimes 
    if growth_dynamic == "logistical" and death_dynamic == "constant" : 
        s1, s2 = s1_s2_num(conversion_timing, c, h, s, graph_type)
    
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Initialization  
    if CI_lenght >= N : print('CI_lenght bigger than the domain lenght')   
    else: edge = np.arange((N-CI_lenght)//2+1, (N+CI_lenght)//2+1); CI_where = np.tile(edge, len(edge))+np.repeat(edge, len(edge))*(N+1)
    W = np.ones((N+1)**2); W[CI_where] = 0    # Wild individuals at t=0  
    H = np.zeros((N+1)**2)                    # Heterozygous individuals at t=0
    D = np.zeros((N+1)**2); D[CI_where] = 1   # Drive individuals at t=0
       
    # First graph (initial condition)
    if graph_type != None and show_graph_ini :
        nb_graph = 1; Z_list = np.zeros((T//mod+2,1000,1000))
        graph_2D(0, W, H, D, N, graph_para, bio_para, num_para, model_para)
        Z_list = graph_2D_contour(0, W, H, D, N, Z_list, nb_graph, graph_para, bio_para, num_para, model_para) 
           
    # Laplacian matrix 
    index_diag_mat = np.arange(0,(N-1)**2+1,N-1)                  # (0, N-1, 2*(N-1), 3*(N-1), ....)
    C0 = -4*np.ones((N-1)**2); C0[np.array([0,N-2,(N-1)**2-(N-1),(N-1)**2-1])]=-2   # place -2
    C0[1:N-2] = -3; C0[(N-1)**2-(N-1)+1:(N-1)**2-1]=-3                              # place -3 in between -2      
    C0[index_diag_mat[1:-2]] = -3; C0[index_diag_mat[2:-1]-1] = -3                                  # place the others -3
    C1 = np.ones((N-1)**2+1); C1[index_diag_mat]=0
    C2 = np.ones((N-1)**2)
    A = sp.spdiags([C2,C1[1:],C0,C1[:-1],C2],[-N+1,-1,0,1,N-1], (N-1)**2, (N-1)**2) # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)
    # Matrix for the explicit side of the Crank Nicholson scheme for W, H and D. 
    Bw = sp.identity((N-1)**2)+((1-theta)*difWW*dt/dx**2)*A           
    Bh = sp.identity((N-1)**2)+((1-theta)*difDW*dt/dx**2)*A            
    Bd = sp.identity((N-1)**2)+((1-theta)*difDD*dt/dx**2)*A     
    # Matrix for the implicit side of the Crank Nicholson scheme for W, H and D.
    Bw_ = sp.identity((N-1)**2)-(theta*difWW*dt/dx**2)*A               
    Bh_ = sp.identity((N-1)**2)-(theta*difDW*dt/dx**2)*A               
    Bd_ = sp.identity((N-1)**2)-(theta*difDD*dt/dx**2)*A        
    # N,S,W,E the four cardinal points (Exterior index)
    index_N = np.arange(N+1,(N+1)*N,N+1); index_S = np.arange(N+1,(N+1)*N,N+1)+N; index_W = np.arange(1,N,1); index_E = np.arange(N*(N+1)+1,(N+1)**2-1,1) 
    index_NW = 0; index_NE = (N+1)**2-(N+1); index_SW = N; index_SE = ((N+1)**2-1)
    index_exterior = np.sort(np.concatenate((index_N, index_S, index_E, index_W, np.array([index_NW, index_NE, index_SW, index_SE]))))
    index_interior = np.array(list(set(np.arange(0,(N+1)**2,1)).difference(set(index_exterior))))
          
    
    # Evolution
    for t in np.linspace(dt,T,M) :
        
        t = round(t,2)
        
        f1 = growth(W,H,D,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[0] - death(W,H,D,bio_para,model_para)*W
        f2 = (1-s*h)*growth(H,W,D,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[1] - death(W,H,D,bio_para,model_para)*H
        f3 = (1-s)*growth(D,H,W,bio_para,model_para)*mating(W,H,D,bio_para,model_para)[2] - death(W,H,D,bio_para,model_para)*D
      
        # Interior
        W[index_interior] = la.spsolve(Bw_, Bw.dot(W[index_interior]) + dt*f1[index_interior])
        H[index_interior] = la.spsolve(Bh_, Bh.dot(H[index_interior]) + dt*f2[index_interior])
        D[index_interior] = la.spsolve(Bd_, Bd.dot(D[index_interior]) + dt*f3[index_interior])
        # Boundaries conditions : Neumann (null derivative)
        W[index_N] = W[index_N+1]; H[index_N] = H[index_N+1]; D[index_N] = D[index_N+1]
        W[index_S] = W[index_S-1]; H[index_S] = H[index_S-1]; D[index_S] = D[index_S-1];
        W[index_W] = W[index_W+N+1]; H[index_W] = H[index_W+N+1]; D[index_W] = D[index_W+N+1];
        W[index_E] = W[index_E-(N+1)]; H[index_E] = H[index_E-(N+1)]; D[index_E] = D[index_E-(N+1)]; 
        W[index_NW] = W[index_NW+1]; H[index_NW] = H[index_NW+1]; D[index_NW] = D[index_NW+1]; 
        W[index_NE] = W[index_NE+1]; H[index_NE] = H[index_NE+1]; D[index_NE] = D[index_NE+1];
        W[index_SW] = W[index_SW-1]; H[index_SW] = H[index_SW-1]; D[index_SW] = D[index_SW-1];
        W[index_SE] = W[index_SE-1]; H[index_SE] = H[index_SE-1]; D[index_SE] = D[index_SE-1];
        
        # Graph
        if t>=mod*nb_graph and graph_type != None : 
            nb_graph += 1
            graph_2D(t, W, H, D, N, graph_para, bio_para, num_para, model_para)
            Z_list = graph_2D_contour(t, W, H, D, N, Z_list, nb_graph, graph_para, bio_para, num_para, model_para) 
                                                   
        # Compute the speed on WT allele proportion or density
        if conversion_timing == "zygote" : 
            WT = (W+0.5*H)/(W+H+D) 
        if conversion_timing == "germline" : 
            WT = (W+0.5*(1-c)*H)/(W+H+D) 
        
        # Check if the WT wave is strictly increasing
        epsilon = -0.01 ; index_neg = np.where(WT[1:]-WT[:-1]<epsilon)[0]
        if len(index_neg) > 4 : print("WT prop decreasing :", index_neg) 
               
        # Check if there is a coexistence state (if so adapt threshold and erase all speed values)
        center_index = np.arange(N//2-N//10,N//2+N//10)
        if np.all((WT[center_index]>0.05) & (WT[center_index]<0.95)) and pw_star < 0 :    # if there exists a palliate in ]0,1[ around N//2 (and pw_star as not been changed yet)
            pw_star = np.mean(WT[center_index])                                            # if this palliate is above the threshold value
            threshold = (pw_star+1)/2
            position = np.array([])   
            speed_fct_of_time = np.array([])   
            time = np.array([]) 
            print("coex ! p* =", pw_star)
            print("threshold =", threshold)
        # we recorde the position only if the WT wave is still in the environment. We do not recorde the 0 position since the threshold value of the wave might be outside the window.
        if np.isin(True, WT>threshold) and np.isin(True, WT<0.99) and np.where(WT>threshold)[0][0] != 0:  
            # List containing the first position of WT where the proportion of wild alleles is higher than the threshold.  
            position = np.append(position, np.where(WT>threshold)[0][0])
            # Wait a minimum time and a minimum number of positions to compute the speed
            if len(position) > 20  : 
                # The speed is computed on the last fifth of the position vector.
                time = np.append(time,t)
                speed_fct_of_time = np.append(speed_fct_of_time, np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt)
        # if the threshold value of the wave is outside the window, stop the simulation  
        if not(np.isin(False, WT>threshold) and np.isin(False, WT<threshold) ) :
            print("t =",t)
            break
 
    
    # plot the speed function of time    
    if len(speed_fct_of_time) != 0 : 
        if graph_type!= None :
            fig, ax = plt.subplots()
            ax.plot(time, speed_fct_of_time) 
            ax.set(xlabel='Time', ylabel='Speed', title = f'Speed function of time')   
            if save_fig :
                directory = f"evolution/{conversion_timing}_r_{np.round(r,3)}_s_{np.round(s,3)}_h_{np.round(h,2)}_c_{np.round(c,2)}"
                save_fig_or_data(directory, fig, speed_fct_of_time, "speed_fct_of_time", bio_para, num_para, model_para)
            plt.show() 
    else :
        print(f"No wave for r = {r} and s = {s}.") 

   # last graph
    if show_graph_end :   
        nb_graph += 1
        graph_2D(t, W, H, D, N, graph_para, bio_para, num_para, model_para) 
        Z_list = graph_2D_contour(t, W, H, D, N, Z_list, nb_graph, graph_para, bio_para, num_para, model_para)                                  

    return(W,H,D,time,speed_fct_of_time)
