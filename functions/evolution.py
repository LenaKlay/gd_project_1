#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:03:42 2021

@author: lena
"""


################## Growth-Death-Mating Functions ##############################

# equ is the allele concerned by the equation.
# oth1 and oth2 are the two other alleles, indifferently.
    
def growth(equ, oth1, oth2, r, a, max_capacity, growth_dynamic, linear_growth):
    n = (equ+oth1+oth2)/max_capacity
    if growth_dynamic == "exponential":
        return(r+1)
    if growth_dynamic == "allee_effect" :
        return(np.max(r*(1-n)*(n-a)+1,0))
    if growth_dynamic == "logistical" and linear_growth == False :
        return(r*(1-n)+1)
    if growth_dynamic == "logistical" and linear_growth == True and max_capacity == 1 :  
        return(r*np.minimum(1-(equ+oth1+oth2),equ) + equ)         
    
def death(equ, oth1, oth2, r, a, max_capacity, death_dynamic):
    n = (equ+oth1+oth2)/max_capacity
    if death_dynamic == "exponential" :
        return(1)
    if death_dynamic == "logistical" :
        return(1+r*n)
    if death_dynamic == "allee_effect" :
        return(r+1+r*(1-n)*(n-a))
     
        
def mating(W, H, D, homing, c, max_capacity, linear_growth, linear_mating):
    
    if linear_mating == False : 
        
        WW = (W**2)/(W+H+D)
        HH = (H**2)/(W+H+D)
        DD = (D**2)/(W+H+D)    
        WH = (W*H)/(W+H+D)   # *2
        WD = (W*D)/(W+H+D)   # *2 
        HD = (D*H)/(W+H+D)   # *2
    
        if homing == "zygote" :       
            mat1 = WW + WH + 0.25*HH 
            mat2 = (1-c)*(WH + 2*WD + 0.5*HH + HD) 
            mat3 = c*WH + 2*c*WD + (0.5*c+0.25)*HH + (1+c)*HD + DD  
            
        if homing == "germline":
            mat1 = WW + (1-c)*WH + 0.25*(1-c)**2*HH 
            mat2 = (1+c)*WH + 2*WD + 0.5*(1-c**2)*HH + (1-c)*HD
            mat3 = 0.25*(c+1)**2*HH + (1+c)*HD + DD
            
    if linear_mating == True and homing == "zygote" and max_capacity == 1 : 
        mat1 = (W>D)
        mat2 = 0
        mat3 = ((W>D)+1)  

    return([mat1,mat2,mat3])          
        
    
    
########################## Evolution Function #################################
                     
        
def evolution(r,s,h,difW,difH,difD,c,T,L,M,N,theta,mod,homing): 
    
    # Steps
    dt = T/M    # time
    dx = L/N    # spatial
    
    # Spatial domain (1D)
    
    X = np.linspace(0,N,N+1)*dx   
    
    # Initialization 
            
    if CI == "center" : 
        W = np.ones(N+1)*max_capacity; W[0:N//2] = 0.05*max_capacity    # Wild individuals at t=0  
        H = np.zeros(N+1)                                               # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//2] = 0.95*max_capacity                # Drive individuals at t=0
    if CI == "left" : 
        W = np.ones(N+1)*max_capacity; W[0:N//10] = 0.05*max_capacity    # Wild individuals at t=0  
        H = np.zeros(N+1)                                                # Heterozygous individuals at t=0
        D = np.zeros(N+1); D[0:N//10] = 0.95*max_capacity                # Drive individuals at t=0
    
    if graph_type != None :
        graph(X,W,H,D,0,graph_type,None,None)
    nb_graph = 1
        
    position = np.array([])   # list containing the first position where the proportion of wild alleles is lower than 0.5.
    vitesse_en_fct_du_tps = np.zeros((2,T//mod))   # first line : time, second line : speed of the wave at the corresponding time
    
    # Matrix
    C0 = -2*np.ones(N+1); C0[0]=C0[0]+1; C0[-1]=C0[-1]+1               
    C1 = np.ones(N+1) 
    A = sp.spdiags([C1,C0,C1],[-1,0,1], N+1, N+1)                # 1D discrete Laplacian with Neumann boundary conditions (derivative=0)  
    
    Bw = sp.identity(N+1)+((1-theta)*difW*dt/dx**2)*A            # Matrix for the explicit side of the Crank Nicholson scheme  
    Bh = sp.identity(N+1)+((1-theta)*difH*dt/dx**2)*A            # for W, H and D.
    Bd = sp.identity(N+1)+((1-theta)*difD*dt/dx**2)*A     

    Bw_ = sp.identity(N+1)-(theta*difW*dt/dx**2)*A               # Matrix for the implicit side of the Crank Nicholson scheme  
    Bh_ = sp.identity(N+1)-(theta*difH*dt/dx**2)*A               # for W, H and D.
    Bd_ = sp.identity(N+1)-(theta*difD*dt/dx**2)*A        

    # Evolution
    for t in np.linspace(dt,T,M) : 
        
        t = round(t,2)
        
        f1 = growth(W,H,D,r,a)*mating(W,H,D)[0] - death(W,H,D,r,a)*W
        f2 = (1-s*h)*growth(H,W,D,r,a)*mating(W,H,D)[1] - death(W,H,D,r,a)*H
        f3 = (1-s)*growth(D,H,W,r,a)*mating(W,H,D)[2] - death(W,H,D,r,a)*D
      
        W = la.spsolve(Bw_, Bw.dot(W) + dt*f1)
        H = la.spsolve(Bh_, Bh.dot(H) + dt*f2)
        D = la.spsolve(Bd_, Bd.dot(D) + dt*f3)
        
        if t>=mod*nb_graph and graph_type != None :  
            graph(X,W,H,D,t,graph_type,None,None)
            nb_graph += 1
        
        if speed_proportion == False :
            if np.isin(True, W>0.5) and np.isin(True, W<0.99) :  # wild pop is still present somewhere in the environment
                position = np.append(position,np.where(W>0.5)[0][0])     # List containing the first position where the 
                                                                         # proportion of wild alleles is lower than 0.5.   
        if speed_proportion == True :                          # Idem on proportion of wild-type individuals
            if np.isin(True, W/(W+H+D)>0.5) and np.isin(True, W/(W+H+D)<0.99) : 
                position = np.append(position,np.where(W/(W+H+D)>0.5)[0][0])
   
        if t%mod == 0 and what_to_do == "vitesse en fonction du temps" : 
            #print("\nTemps :", t) ; print("Vitesse :", np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt, "\n")
            vitesse_en_fct_du_tps[:,int(t)//mod-1] = np.array([t,np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt])
            # first line : time, second line : speed of the wave at the corresponding time

    if np.shape(position)[0] != 0 :
        position = position[position!=0]   
        speed = np.mean(np.diff(position[int(4*len(position)/5):len(position)]))*dx/dt  # Speed of the wave   
    else :
        speed = None
        print(f"Can't determine the speed of the wave : For r = {r} and s = {s}, the environment is empty or not mixed enought (at least one place where W<0.5 and one where W>0.1), or T is too short.")                                                 

    if speed != None and graph_type != None and theorical_wave == True :  
            graph(X,W,H,D,t,graph_type,speed,origin)   

    if what_to_do == "vitesse en fonction du temps" : 
        return(speed,W,H,D,vitesse_en_fct_du_tps)
    if what_to_do == "system+evolution" : 
        return(speed,position,W,H,D)
    else :
        return(speed,W,H,D)


