#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# data utils module

import numpy as np
from scipy.stats import norm, normaltest, shapiro, chisquare, kstest
from statsmodels.stats.diagnostic import lilliefors
from oct2py import octave


def linear_algebra_check(name):
    A = octave.load('./octave/matrix_'+ name + '.mat')
    b = octave.load('./octave/vector_'+ name + '.mat')
    return octave.linsolve(A,b)


def normal_test(sample,alpha,verbose):
    # hypothesis test: null hypothesis, the data is gaussian distributed

    # Shapiro-Wilk
    stat, p = shapiro(sample)
    if verbose:
        if p > alpha: print('Shapiro this is Gaussian', p)
        else:         print('Shapiro this is NOT Gaussian', p)

    # chisquare
    stat, p = chisquare(sample)
    if verbose:
        if p > alpha: print('Chisquare this is Gaussian', p)
        else:         print('Chisquare this is NOT Gaussian', p)

    # lilliefors
    stat, p = lilliefors(sample)
    if verbose:
        if p > alpha: print('Lilliefors this is Gaussian', p)
        else:         print('Lilliefors this is NOT Gaussian', p)

    # kolmogorov
    stat, p = kstest(sample, 'norm')
    if verbose:
        if p > alpha: print('Kolmogorov this is Gaussian', p)
        else:         print('Kolmogorov this is NOT Gaussian', p)

    # Angostino
    k2, p = normaltest(sample)
    if verbose:
        if p > alpha: print('Angostino this is Gaussian', p)
        else:         print('Angostino this is NOT Gaussian', p)    

    return p,alpha


def names_values(name):
    if name=='All':
        quantum = ['$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,3/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,1,0,10/3 \\rangle$ &',
                   '$\\vert ssc,1/2,3/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',
                   '$\\vert suc,1/2,1/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,1/2,3/2,1,1/2,10/3 \\rangle$ &',
                   '$\\vert suc,3/2,1/2,1,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,1,1/2,10/3 \\rangle$ &',
                   '$\\vert uuc,1/2,1/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,3/2,3/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,1/2,1/2,1,1,10/3 \\rangle$ &',
                   '$\\vert udc,1/2,1/2,0,0,4/3  \\rangle$ &',   '$\\vert udc,1/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert udc,3/2,1/2,1,0,4/3  \\rangle$ &',
                   '$\\vert ssc,1/2,1/2,0,1/2,4/3 \\rangle$ &',  '$\\vert ssc,1/2,1/2,1,1/2,4/3 \\rangle$ &',  '$\\vert ssc,3/2,1/2,1,1/2,10/3 \\rangle$ &' ]

        old = np.array([2702.4, 2767.0, 3015.8, 3044.5, 3051.6, 3080.4, 2570.1, 2634.8, 2934.1, 2941.2,
                        2969.9, 2453.1, 2517.7, 2819.0, 2283.7, 2649.7, 2685.6, 2461.2, 2796.5, 2832.4])
        exp = np.array([2695.0, 2766.0, 3000.4, 3050.2, 3065.6, 3090.2, 2578.0, 2645.9, 2923.0, 2938.5,
                        2964.9, 2453.9, 2518.0, 2801.0, 2286.5, 2592.3, 2625.0, 2469.0, 2792.0, 2815.0])
        delta_exp = np.array([2.0,2.0,0.4,0.3,0.4,0.7,2.9,0.6,0.4,0.3,0.3,0.1,2.3,6.0,0.1,0.4,0.2,4.0,3.3,0.2])

    if name=='omega':        
        quantum = ['$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,3/2,0,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,1,0,10/3 \\rangle$ &',
                   '$\\vert ssc,1/2,3/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,3/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert ssc,1/2,1/2,0,0,10/3 \\rangle$ &']
        old = np.array([2702.4, 2767.0, 3015.8, 3044.5, 3051.6, 3080.4])
        exp = np.array([2695.0, 2766.0, 3000.4, 3050.2, 3065.6, 3090.2])
        delta_exp = np.array([2.0,2.0,0.4,0.3,0.4,0.7])

    if name=='cascades':
        quantum = ['$\\vert suc,1/2,1/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,0,1/2,10/3 \\rangle$ &', '$\\vert suc,1/2,3/2,1,1/2,10/3 \\rangle$ &',
                   '$\\vert suc,3/2,1/2,1,1/2,10/3 \\rangle$ &', '$\\vert suc,3/2,3/2,1,1/2,10/3 \\rangle$ &',                   
                   '$\\vert ssc,1/2,1/2,0,1/2,4/3 \\rangle$ &',  '$\\vert ssc,1/2,1/2,1,1/2,4/3 \\rangle$ &',  '$\\vert ssc,3/2,1/2,1,1/2,10/3 \\rangle$ &' ]
        old = np.array([2570.1, 2634.8, 2934.1, 2941.2, 2969.9, 2461.2, 2796.5, 2832.4])
        exp = np.array([2578.0, 2645.9, 2923.0, 2938.5, 2964.9, 2469.0, 2792.0, 2815.0])
        delta_exp = np.array([2.9,0.6,0.4,0.3,0.3,4.0,3.3,0.2])

    if name=='sigmaLamb':
        quantum = ['$\\vert uuc,1/2,1/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,3/2,3/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,1/2,1/2,1,1,10/3 \\rangle$ &',
                   '$\\vert udc,1/2,1/2,0,0,4/3  \\rangle$ &',   '$\\vert udc,1/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert udc,3/2,1/2,1,0,4/3  \\rangle$ &']
        old = np.array([2453.1, 2517.7, 2819.0, 2283.7, 2649.7, 2685.6])
        exp = np.array([2453.9, 2518.0, 2801.0, 2286.5, 2592.3, 2625.0])
        delta_exp = np.array([0.1,2.3,6.0,0.1,0.4,0.2])
        
    return quantum, old, exp, delta_exp
