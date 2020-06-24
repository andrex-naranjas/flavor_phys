#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Code to obtain mass uncertainties of quarkonium mass spectrum
#author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# main module
from iminuit import Minuit
import numpy as np
from sklearn.utils import resample

# call module for plots
import data_visualization as dv

# input parameters 
param_o = np.array([0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00]) # coef infront omega_lambda, omega_rho
param_x = np.array([0.75, 3.75, 0.75, 3.75, 0.75, 3.75, 3.75]) # coef infront A
param_y = np.array([0.00, 0.00, -1.0, -2.5, 0.50, -1.0, 1.50]) # coef infront B
param_z = np.array([3.33, 3.33, 3.33, 3.33, 3.33, 3.33, 3.33]) # coef infront G
#exp_m   = np.array([2695, 2766, 3000, 3050, 3066, 3090, 3188]) # experimental masses

def model(w,x,y,z,o,a,b,g):
    return 2505 + w*o + x*a + y*b + z*g

def least_squares(o, a, b, g):
    yvar = 0.01
    comp_m = model(param_o, param_x, param_y, param_z, o, a, b, g)
    return np.sum((comp_m - exp_m)**2 / yvar)

def fit(least_squares):
    m = Minuit(least_squares, o=0, a=0, b=0, g=0, error_o=1, error_a=1, error_b=1, error_g = 1, errordef=1)
    m.migrad()
    return m

def sample_gauss(mu, sigma):
    return np.random.normal(mu, sigma, 10000)

def random(sample):
    return np.random.choice(sample, size=1).mean()


# arrays to store the sampled parameters
sampled_o,sampled_a,sampled_b,sampled_g = ([]),([]),([]),([])

# gaussian pdf with the measured value and with experimental uncertainty
gauss_2695 = sample_gauss(2695.0, 2.0000)
gauss_2770 = sample_gauss(2766.0, 2.0000)
gauss_3000 = sample_gauss(3000.4, 0.3742)
gauss_3050 = sample_gauss(3050.2, 0.3317)
gauss_3066 = sample_gauss(3065.6, 0.4359)
gauss_3090 = sample_gauss(3090.2, 0.6557)
gauss_3180 = sample_gauss(3188.0, 13.923)


# construct the simulated sampling distribution (boostrap technique)
for _ in range(100000):
    # experimental sampled masses
    exp_m = np.array([random(gauss_2695), random(gauss_2770), random(gauss_3000),
                      random(gauss_3050), random(gauss_3066), random(gauss_3090),
                      random(gauss_3180)])

    # perform the parameter fitting (via minimizing squared distance)
    m = fit(least_squares)

    sampled_o = np.append(sampled_o, m.values['o'])
    sampled_a = np.append(sampled_a, m.values['a'])
    sampled_b = np.append(sampled_b, m.values['b'])
    sampled_g = np.append(sampled_g, m.values['g'])

print(np.mean(sampled_o),np.mean(sampled_a),np.mean(sampled_b),np.mean(sampled_g))
print(np.std(sampled_o),np.std(sampled_a),np.std(sampled_b),np.std(sampled_g))

# plot the simulated sampling distribution,
# under the Central Limit Theorem, it is expected normal
dv.plot(sampled_o,'o','Parameter omega','Charm')
dv.plot(sampled_a,'a','Parameter A','Charm')
dv.plot(sampled_b,'b','Parameter B','Charm')
dv.plot(sampled_g,'g','Parameter G','Charm')
