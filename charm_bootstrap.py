#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# authors: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)
#          H. Garcia-Tecocoatzi

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
    yvar = 0.1
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
gauss_2695 = sample_gauss(2695.00, 1.7)  # OK
gauss_2770 = sample_gauss(2765.90, 2.0)  # OK
gauss_3000 = sample_gauss(3000.41, 0.22) # OK
gauss_3050 = sample_gauss(3050.20, 0.13) # OK
gauss_3066 = sample_gauss(3065.46, 0.28) # OK
gauss_3090 = sample_gauss(3090.00, 0.5)  # OK
gauss_3120 = sample_gauss(3119.10, 1.0)  # need to the fit (corresponds to predicted 3140)

# Cascades six-plet
gauss_2578 = sample_gauss(2578.4, 0.5) # OK, 2579.2 ± 0.5 MeV
gauss_2645 = sample_gauss(2645.56,0.27)# OK
gauss_2923 = sample_gauss(2923.1, 0.4) # OK (not in PDG)
gauss_2938 = sample_gauss(2938.6, 0.3) # OK (not in PDG) is it the same as 2942???
gauss_2970 = sample_gauss(2969.9, 0.5) # OK 

gauss_2942 = sample_gauss(2942.0, 5.0) # need add to fit (corresponds to predicted 2941 )  2** PDG,Sig(2930)
gauss_3055 = sample_gauss(3055.9, 0.4) # need add to fit (corresponds to predicted 3060)
gauss_3080 = sample_gauss(3079.9, 1.4) # need add to fit (corresponds to predicted 3096)


# Sigma c
gauss_2453 = sample_gauss(2453.07, 0.4) # OK (corresponds to predicted 2453)
gauss_2518 = sample_gauss(2517.5, 2.3)  # OK (corresponds to predicted 2517)
gauss_2801 = sample_gauss(2801.0, 5)    # OK (corresponds to predicted 2819)
# Lambda C
gauss_2286 = sample_gauss(2286.46, 0.14) # OK (corresponds to 2283)
gauss_2592 = sample_gauss(2592.25, 0.28) # OK (corresponds to 2649)
gauss_2628 = sample_gauss(2628.11, 0.19) # OK (corresponds to 2685)
# Cascade C anti-3-plet
gauss_2469 = sample_gauss(2470.9,  0.29) # OK (corresponds to predicted 2461)
gauss_2792 = sample_gauss(2792.4,  0.5)  # OK (corresponds to predicted 2796)
gauss_2816 = sample_gauss(2816.74, 0.23) # OK (corresponds to predicted 2832)

# construct the simulated sampling distribution (bootstrap technique)
for _ in range(100):
    # measured and quark-sum sampled masses
    if(states=='All'):
        exp_m = np.array([random(gauss_2695), random(gauss_2770), random(gauss_3000),
                          random(gauss_3050), random(gauss_3066), random(gauss_3090),
                          random(gauss_2578), random(gauss_2645), random(gauss_2923),
                          random(gauss_2938), random(gauss_2970), random(gauss_2453),
                          random(gauss_2518), random(gauss_2801), random(gauss_2286),
                          random(gauss_2592), random(gauss_2628), random(gauss_2469),
                          random(gauss_2792), random(gauss_2816)])        
    elif(states=='omega'):
        exp_m = np.array([random(gauss_2695), random(gauss_2770), random(gauss_3000),
                          random(gauss_3050), random(gauss_3066), random(gauss_3090)])
    elif(states=='cascades'):
        exp_m = np.array([random(gauss_2578), random(gauss_2645), random(gauss_2923),
                          random(gauss_2938), random(gauss_2970), random(gauss_2469),
                          random(gauss_2792), random(gauss_2816)])
    elif(states=='sigmaLamb'):
        exp_m = np.array([random(gauss_2453), random(gauss_2518), random(gauss_2801),
                          random(gauss_2286), random(gauss_2592), random(gauss_2628)])
               
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
results.mass_prediction_compare(bootstrap=True)
results.correlation_matrix()
results.param_comparison()
results.plot()

# omegas,cascades,sigmas,lambdas,cascades_anti3
results.paper_results_predictions(baryons='omegas',         bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True)
results.paper_results_predictions(baryons='cascades',       bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True)
results.paper_results_predictions(baryons='sigmas',         bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True)
results.paper_results_predictions(baryons='lambdas',        bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True)
results.paper_results_predictions(baryons='cascades_anti3', bootstrap=True, bootstrap_width=True, prev_params=False, decay_width=True)


# nominal results, expected to be the same as the previous paper, very important check
# results.paper_results_predictions(baryons='omegas',         bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='cascades',       bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='sigmas',         bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='lambdas',        bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)
# results.paper_results_predictions(baryons='cascades_anti3', bootstrap=False, bootstrap_width=False, prev_params=True, decay_width=True)


# quark-sum masses
# gauss_2505 = sample_gauss(2505, 10.0)
# gauss_2350 = sample_gauss(2350, 10.0)
# gauss_2195 = sample_gauss(2195, 10.0)
# mass_sum = np.array([random(gauss_2505),random(gauss_2505),random(gauss_2505),
#                      random(gauss_2505),random(gauss_2505),random(gauss_2505),
#                      random(gauss_2350),random(gauss_2350),random(gauss_2350),
#                      random(gauss_2350),random(gauss_2350),
#                      random(gauss_2195),random(gauss_2195),random(gauss_2195),
#                      random(gauss_2195),random(gauss_2195),random(gauss_2195),
#                      random(gauss_2350),random(gauss_2350),random(gauss_2350)])
# mass_sum = np.array([random(gauss_2505),random(gauss_2505),random(gauss_2505),
#                      random(gauss_2505),random(gauss_2505),random(gauss_2505)])
# mass_sum = np.array([random(gauss_2350),random(gauss_2350),random(gauss_2350),
#                      random(gauss_2350),random(gauss_2350),random(gauss_2350),
#                      random(gauss_2350),random(gauss_2350)])
# mass_sum = np.array([random(gauss_2195),random(gauss_2195),random(gauss_2195),
#                      random(gauss_2195),random(gauss_2195),random(gauss_2195)])

        
