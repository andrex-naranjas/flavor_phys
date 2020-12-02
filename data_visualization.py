#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# data visalization module

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

import data_utils as du

# sample plots
def plot(sample,name,xlabel,quark, states):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(sample, 100, density=True, label = 'Sampling')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    y = norm.pdf(x, np.mean(sample), np.std(sample))
    p, alpha = du.normal_test(sample, alpha=0.05, verbose=False)        
    plt.plot(x,y,label='Gaussian fit')
    plt.text(0.15, 0.9,'$\\mu$={}, $\\sigma$={}'.format(round(np.mean(sample),1), round(np.std(sample),1)),
             ha='center', va='center', transform=ax.transAxes)
    plt.text(0.15, 0.8,'$p_{{val}}$={}, $\\alpha$={}'.format(round(p,3), alpha),
             ha='center', va='center', transform=ax.transAxes)
    plt.legend(frameon=False)
    plt.xlabel(xlabel)
    plt.ylabel('Arbitrary Units')
    plt.title(quark+' mesons')
    plt.savefig('./plots/'+quark+'_bootstrap_'+name+'_'+states+'.pdf')
    plt.close()
