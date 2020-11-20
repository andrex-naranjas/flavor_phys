#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# data preparation module

import numpy as np

def fetch_data(states):
    
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
        param_v = np.array([0.00, 0.00, w_om, w_om, w_om, w_om])
        param_w = np.array([0.75, 3.75, 0.75, 3.75, 0.75, 3.75])
        param_x = np.array([0.00, 0.00, -1.0, -2.5, 0.50, -1.0])
        param_y = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
        param_z = np.array([10/3, 10/3, 10/3, 10/3, 10/3, 10/3])
        mass_sum= np.array([2505, 2505, 2505, 2505, 2505, 2505])
        return param_v,param_w,param_x,param_y,param_z,mass_sum
    elif(states=='cascades'):
        param_v = np.array([0.00, 0.00, wcas, wcas, wcas, 0.00, wcas, wcas])
        param_w = np.array([0.75, 3.75, 3.75, 0.75, 3.75, 0.75, 0.75, 0.75])
        param_x = np.array([0.00, 0.00, -2.5, 0.50, -1.0, 0.00, -1.0, 0.50])
        param_y = np.array([0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.75])
        param_z = np.array([10/3, 10/3, 10/3, 10/3, 10/3, 4/3 , 4/3 , 4/3 ])
        mass_sum= np.array([2350, 2350, 2350, 2350, 2350, 2350, 2350, 2350])
        return param_v,param_w,param_x,param_y,param_z,mass_sum
    elif(states=='sigma_lamb'):
        param_v = np.array([0.00, 0.00, w_ls, 0.00, w_ls, w_ls]) # coef infront kprim_c
        param_w = np.array([0.75, 3.75, 0.75, 0.75, 0.75, 0.75]) # coef infront A
        param_x = np.array([0.00, 0.00, -1.0, 0.00, -1.0, 0.50]) # coef infront B
        param_y = np.array([2.00, 2.00, 2.00, 0.00, 0.00, 0.00]) # coef infront E
        param_z = np.array([10/3, 10/3, 10/3, 4/3 , 4/3 , 4/3 ]) # coef infront G
        mass_sum= np.array([2195, 2195, 2195, 2195, 2195, 2195]) # sum of masses
        return param_v,param_w,param_x,param_y,param_z,mass_sum
