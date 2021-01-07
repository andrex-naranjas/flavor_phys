#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# data preparation module

import numpy as np

def fetch_data(states):
    # selected states to do fit
    w_om=0.058892513
    w_ls=0.068087711
    wcas=0.062695053

    if(states == 'All'):
        param_v = np.array([0.00, 0.00, w_om, w_om, w_om, w_om, 0.00, 0.00, wcas, wcas, wcas, 0.00, 0.00, w_ls, 0.00, w_ls, w_ls, 0.00, wcas, wcas]) # coef infront kprim_c
        param_w = np.array([0.75, 3.75, 0.75, 3.75, 0.75, 3.75, 0.75, 3.75, 3.75, 0.75, 3.75, 0.75, 3.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75]) # coef infront A
        param_x = np.array([0.00, 0.00, -1.0, -2.5, 0.50, -1.0, 0.00, 0.00, -2.5, 0.50, -1.0, 0.00, 0.00, -1.0, 0.00, -1.0, 0.50, 0.00, -1.0, 0.50]) # coef infront B
        param_y = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.75, 0.75, 0.75, 0.75, 0.75, 2.00, 2.00, 2.00, 0.00, 0.00, 0.00, 0.75, 0.75, 0.75]) # coef infront E
        param_z = np.array([10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 10/3, 4/3 , 4/3 , 4/3 , 4/3 , 4/3 , 4/3 ]) # coef infront G
        mass_sum= np.array([2505, 2505, 2505, 2505, 2505, 2505, 2350, 2350, 2350, 2350, 2350, 2195, 2195, 2195, 2195, 2195, 2195, 2350, 2350, 2350]) # sum of masses
        return param_v,param_w,param_x,param_y,param_z,mass_sum
    elif(states == 'omega'):
        param_v = np.array([0.00, 0.00, w_om, w_om, w_om, w_om]) # coef infront kprim_c  
        param_w = np.array([0.75, 3.75, 0.75, 3.75, 0.75, 3.75]) # coef infront A        
        param_x = np.array([0.00, 0.00, -1.0, -2.5, 0.50, -1.0]) # coef infront B        
        param_y = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00]) # coef infront E        
        param_z = np.array([10/3, 10/3, 10/3, 10/3, 10/3, 10/3]) # coef infront G        
        mass_sum= np.array([2505, 2505, 2505, 2505, 2505, 2505]) # sum of masses         
        return param_v,param_w,param_x,param_y,param_z,mass_sum
    elif(states=='cascades'):
        param_v = np.array([0.00, 0.00, wcas, wcas, wcas, 0.00, wcas, wcas]) # coef infront kprim_c  
        param_w = np.array([0.75, 3.75, 3.75, 0.75, 3.75, 0.75, 0.75, 0.75]) # coef infront A        
        param_x = np.array([0.00, 0.00, -2.5, 0.50, -1.0, 0.00, -1.0, 0.50]) # coef infront B        
        param_y = np.array([0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75]) # coef infront E        
        param_z = np.array([10/3, 10/3, 10/3, 10/3, 10/3, 4/3 , 4/3 , 4/3 ]) # coef infront G        
        mass_sum= np.array([2350, 2350, 2350, 2350, 2350, 2350, 2350, 2350]) # sum of masses         
        return param_v,param_w,param_x,param_y,param_z,mass_sum
    elif(states=='sigmaLamb'):
        param_v = np.array([0.00, 0.00, w_ls, 0.00, w_ls, w_ls]) # coef infront kprim_c
        param_w = np.array([0.75, 3.75, 0.75, 0.75, 0.75, 0.75]) # coef infront A
        param_x = np.array([0.00, 0.00, -1.0, 0.00, -1.0, 0.50]) # coef infront B
        param_y = np.array([2.00, 2.00, 2.00, 0.00, 0.00, 0.00]) # coef infront E
        param_z = np.array([10/3, 10/3, 10/3, 4/3 , 4/3 , 4/3 ]) # coef infront G
        mass_sum= np.array([2195, 2195, 2195, 2195, 2195, 2195]) # sum of masses            
        return param_v,param_w,param_x,param_y,param_z,mass_sum
    

def hamiltonian_quantum_factors(state, sum_mass, J_tot, S_tot, L_tot, I_tot, SU_tot, HO_n, ModEx):
    # compute the (quantum) coeficients that multiplied the fitted/bootstrap parameters
    # this function needs to be fed directly with quantum numbers storeds in a numpy.array
    
    param_v,param_w,param_x,param_y,param_z,mass_sum = [],[],[],[],[],[]
    
    for i in range(len(sum_mass)):    
        omega_ho=0
        if(ModEx[i]=='lam'):
            if(state[i]=='omg'):   omega_ho=0.058892513
            elif(state[i]=='cas'): omega_ho=0.062695053
            elif(state[i]=='tri'): omega_ho=0.062695053            
            elif(state[i]=='sig'): omega_ho=0.068087711
            elif(state[i]=='Lam'): omega_ho=0.068087711
        elif(ModEx[i]=='rho'):
            if(state[i]=='omg'):   omega_ho=0.081649658
            elif(state[i]=='cas'): omega_ho=0.089742361
            elif(state[i]=='tri'): omega_ho=0.089742361
            elif(state[i]=='sig'): omega_ho=0.100843897
            elif(state[i]=='Lam'): omega_ho=0.100843897
                
        param_v.append(HO_n[i]  * omega_ho )    # coef infront kprim 
        param_w.append((S_tot[i] + 1)*S_tot[i]) # coef infront A
        param_x.append(0.5 * ( (J_tot[i] + 1)*J_tot[i] - (L_tot[i] + 1)*L_tot[i] - (S_tot[i] + 1)*S_tot[i] ) ) # coef infront B
        param_y.append((I_tot[i] + 1)*I_tot[i]) # coef infront E
        param_z.append(SU_tot[i])               # coef infront G
        mass_sum.append(sum_mass[i])            # sum of masses

    param_v = np.array(param_v)
    param_w = np.array(param_w)
    param_x = np.array(param_x)
    param_y = np.array(param_y)
    param_z = np.array(param_z)
    mass_sum= np.array(mass_sum)
    
    return param_v,param_w,param_x,param_y,param_z,mass_sum
