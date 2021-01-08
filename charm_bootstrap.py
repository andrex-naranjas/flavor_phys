#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# main module

# standard modules
import sys
from iminuit import Minuit
import numpy as np

# bootstrap
from sklearn.utils import resample

# framework modules
import data_visualization as dv
import data_preparation as dp
from data_results import CharmResults


if len(sys.argv) != 2:
    sys.exit('Provide charm states group name. Try again!')

#states = 'omega' # All, omega, cascades, sigmaLamb
states = sys.argv[1]

# input parameters
param_v,param_w,param_x,param_y,param_z,mass_sum = dp.fetch_data(states)

def model(mass_sum,v,w,x,y,z,k,a,b,e,g):
    return mass_sum + v*k + w*a + x*b + y*e + z*g

def least_squares(k, a, b, e, g):
    yvar = 1.
    pred_m = model(mass_sum, param_v, param_w, param_x, param_y, param_z, k, a, b, e, g)
    return np.sum((pred_m - exp_m)**2 / (yvar**2)) #**2

def fit(least_squares):
    m = Minuit(least_squares, k=0, a=0, b=0, e=0, g=0, error_k=1, error_a=1, error_b=1, error_e=1, error_g=1, errordef=1)
    m.migrad()
    return m

def sample_gauss(mu, sigma):
    return np.random.normal(mu, sigma, 10000)

def random(sample):
    return np.random.choice(sample, size=None)

# arrays to store the sampled parameters
sampled_k,sampled_a,sampled_b,sampled_e,sampled_g,sampled_masses = ([]),([]),([]),([]),([]),([])

# arrays to store sampled correlation coeficients
rho_ak,rho_bk,rho_ba,rho_ek,rho_ea,rho_eb,rho_gk,rho_ga,rho_gb,rho_ge=([]),([]),([]),([]),([]),([]),([]),([]),([]),([])

# gaussian pdf with the measured value and with experimental uncertainty
# Omega states
scale = 1
gauss_2695 = sample_gauss(2695.0, 2.0*scale)
gauss_2770 = sample_gauss(2766.0, 2.0*scale)
gauss_3000 = sample_gauss(3000.4, 0.3742*scale)
gauss_3050 = sample_gauss(3050.2, 0.3317*scale)
gauss_3066 = sample_gauss(3065.6, 0.4359*scale)
gauss_3090 = sample_gauss(3090.2, 0.6557*scale)
# Cascades six-plet
gauss_2578 = sample_gauss(2578.0, 2.9*scale)
gauss_2645 = sample_gauss(2645.9, 0.6*scale)
gauss_2923 = sample_gauss(2923.1, 0.4*scale)
gauss_2938 = sample_gauss(2938.6, 0.3*scale)
gauss_2964 = sample_gauss(2964.9, 0.3*scale)
# Sigma c
gauss_2453 = sample_gauss(2453.9, 0.14*scale)
gauss_2518 = sample_gauss(2518.0, 2.3*scale)
gauss_2801 = sample_gauss(2801.0, 6*scale)
# Lambda C
gauss_2286 = sample_gauss(2286.5, 0.14*scale)
gauss_2592 = sample_gauss(2592.3, 0.28*scale)
gauss_2628 = sample_gauss(2625.0, 0.19*scale)
# Cascade C anti-3-plet
gauss_2469 = sample_gauss(2469.0, 4*scale)
gauss_2792 = sample_gauss(2792.0, 3.3*scale)
gauss_2816 = sample_gauss(2815.0, 0.20*scale)

# quark-sum masses
gauss_2505 = sample_gauss(2505, 10.0)
gauss_2350 = sample_gauss(2350, 10.0)
gauss_2195 = sample_gauss(2195, 10.0)

# construct the simulated sampling distribution (bootstrap technique)
for _ in range(100):
    # measured and quark-sum sampled masses
    if(states=='All'):
        exp_m = np.array([random(gauss_2695), random(gauss_2770), random(gauss_3000),
                          random(gauss_3050), random(gauss_3066), random(gauss_3090),
                          random(gauss_2578), random(gauss_2645), random(gauss_2923),
                          random(gauss_2938), random(gauss_2964), random(gauss_2453),
                          random(gauss_2518), random(gauss_2801), random(gauss_2286),
                          random(gauss_2592), random(gauss_2628), random(gauss_2469),
                          random(gauss_2792), random(gauss_2816)])
        # mass_sum = np.array([random(gauss_2505),random(gauss_2505),random(gauss_2505),
        #                      random(gauss_2505),random(gauss_2505),random(gauss_2505),
        #                      random(gauss_2350),random(gauss_2350),random(gauss_2350),
        #                      random(gauss_2350),random(gauss_2350),
        #                      random(gauss_2195),random(gauss_2195),random(gauss_2195),
        #                      random(gauss_2195),random(gauss_2195),random(gauss_2195),
        #                      random(gauss_2350),random(gauss_2350),random(gauss_2350)])
        
    elif(states=='omega'):
        exp_m = np.array([random(gauss_2695), random(gauss_2770), random(gauss_3000),
                          random(gauss_3050), random(gauss_3066), random(gauss_3090)])
        # mass_sum = np.array([random(gauss_2505),random(gauss_2505),random(gauss_2505),
        #                      random(gauss_2505),random(gauss_2505),random(gauss_2505)])

    elif(states=='cascades'):
        exp_m = np.array([random(gauss_2578), random(gauss_2645), random(gauss_2923),
                          random(gauss_2938), random(gauss_2964), random(gauss_2469),
                          random(gauss_2792), random(gauss_2816)])
        # mass_sum = np.array([random(gauss_2350),random(gauss_2350),random(gauss_2350),
        #                      random(gauss_2350),random(gauss_2350),random(gauss_2350),
        #                      random(gauss_2350),random(gauss_2350)])
        
    elif(states=='sigmaLamb'):
        exp_m = np.array([random(gauss_2453), random(gauss_2518), random(gauss_2801),
                          random(gauss_2286), random(gauss_2592), random(gauss_2628)])
        # mass_sum = np.array([random(gauss_2195),random(gauss_2195),random(gauss_2195),
        #                      random(gauss_2195),random(gauss_2195),random(gauss_2195)])

               
    # perform the parameter fitting (via minimizing squared distance)
    m = fit(least_squares)

    sampled_k = np.append(sampled_k, m.values['k'])
    sampled_a = np.append(sampled_a, m.values['a'])
    sampled_b = np.append(sampled_b, m.values['b'])
    sampled_e = np.append(sampled_e, m.values['e'])
    sampled_g = np.append(sampled_g, m.values['g'])

    if states != 'omega':
        # correlation matrix
        corr = m.np_matrix(correlation=True)
        
        rho_ak = np.append(rho_ak, corr[1,0])
        rho_bk = np.append(rho_bk, corr[2,0])
        rho_ba = np.append(rho_ba, corr[2,1])
        rho_ek = np.append(rho_ek, corr[3,0])
        rho_ea = np.append(rho_ea, corr[3,1])
        rho_eb = np.append(rho_eb, corr[3,2])
        rho_gk = np.append(rho_gk, corr[4,0])
        rho_ga = np.append(rho_ga, corr[4,1])
        rho_gb = np.append(rho_gb, corr[4,2])
        rho_ge = np.append(rho_ge, corr[4,3])

    # store the sampled experimental masses
    sampled_masses = np.append(sampled_masses, exp_m)


# create dictionaries
param   = {'mass':mass_sum, 'V':param_v, 'W':param_w, 'X':param_x, 'Y':param_y, 'Z':param_z}
sampled = {'sampled_k':sampled_k, 'sampled_a':sampled_a, 'sampled_b':sampled_b, 'sampled_e':sampled_e, 'sampled_g':sampled_g}
corr_mat= {'rho_ak':rho_ak,'rho_bk':rho_bk,'rho_ba':rho_ba,'rho_ek':rho_ek,'rho_ea':rho_ea,'rho_eb':rho_eb,'rho_gk':rho_gk,'rho_ga':rho_ga,'rho_gb':rho_gb,'rho_ge':rho_ge}

# calculate the results using bootstrap simulation above
results = CharmResults(param, sampled, corr_mat, asymmetric=True, name=states)
#results.mass_prediction_compare()
results.correlation_matrix()
results.param_comparison()
results.plot()
# results.execute_decay_width()

# omegas,cascades,sigmas,lambdas,cascades_anti3
results.paper_results_predictions(baryons='omegas', bootstrap=True, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='cascades', bootstrap=True, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='sigmas', bootstrap=True, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='lambdas', bootstrap=True, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='cascades_anti3', bootstrap=True, prev_params=True, decay_width=True)
