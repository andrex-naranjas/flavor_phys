#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Code to obtain mass uncertainties of quarkonium mass spectrum
#author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

#data visalization module

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm

# sample plots
def plot(sample,name,xlabel,quark):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(sample, 100, density=True, label = 'Sampling')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    y = norm.pdf(x, np.mean(sample), np.std(sample))
    plt.plot(x,y,label='Gaussian fit')
    plt.text(0.15, 0.9,'$\mu$={}, $\sigma$={}'.format(round(np.mean(sample),1), round(np.std(sample),1)),
             ha='center', va='center', transform=ax.transAxes)
    plt.legend(frameon=False)
    plt.xlabel(xlabel)
    plt.ylabel('Arbitrary Units')
    plt.title(quark+' mesons')
    plt.savefig('./plots/'+quark+'_bootstrap_'+name+'.png')
    plt.close()
