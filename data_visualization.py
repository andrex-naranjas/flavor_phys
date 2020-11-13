#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# data visalization module

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, normaltest, shapiro, chisquare, kstest
from statsmodels.stats.diagnostic import lilliefors
from oct2py import octave

# sample plots
def plot(sample,name,xlabel,quark, states):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(sample, 100, density=True, label = 'Sampling')
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    y = norm.pdf(x, np.mean(sample), np.std(sample))
    p, alpha = normal_test(sample, alpha=0.05, verbose=True)        
    plt.plot(x,y,label='Gaussian fit')
    plt.text(0.15, 0.9,'$\\mu$={}, $\\sigma$={}'.format(round(np.mean(sample),1), round(np.std(sample),1)),
             ha='center', va='center', transform=ax.transAxes)
    plt.text(0.15, 0.8,'$p_{{val}}$={}, $\\alpha$={}'.format(round(p,3), alpha),
             ha='center', va='center', transform=ax.transAxes)
    plt.legend(frameon=False)
    plt.xlabel(xlabel)
    plt.ylabel('Arbitrary Units')
    plt.title(quark+' mesons')
    plt.savefig('./plots/'+quark+'_bootstrap_TEST_'+name+'_'+states+'.pdf')
    plt.close()


# mass prediction
def mass_prediction(sum_mass, v_param, w_param, x_param, y_param, z_param,
                    k_sampled, a_sampled, b_sampled, e_sampled, g_sampled,
                    rho_ak,rho_bk,rho_ba,rho_ek,rho_ea,rho_eb,rho_gk,rho_ga,rho_gb,rho_ge,
                    bootstrap, asymmetric, name):

    Kp, delta_Kp = np.mean(k_sampled), np.std(k_sampled)
    A, delta_A = np.mean(a_sampled), np.std(a_sampled)
    B, delta_B = np.mean(b_sampled), np.std(b_sampled)
    E, delta_E = np.mean(e_sampled), np.std(e_sampled)
    G, delta_G = np.mean(g_sampled), np.std(g_sampled)

    if bootstrap:
        bootstrap_masses, bootstrap_errors, delta_up, delta_dn = sampled_prediction(sum_mass, v_param, w_param, x_param, y_param, z_param,
                                                                                      k_sampled, a_sampled, b_sampled, e_sampled, g_sampled,
                                                                                      rho_ak,rho_bk,rho_ba,rho_ek,rho_ea,rho_eb,rho_gk,rho_ga,rho_gb,rho_ge, name)       
    quantum, old, exp, delta_exp = values(name)
        
    LA = linear_algebra(name)
    f = open('./tables/parameters_'+name+'.tex', "w")
    
    K = np.sqrt(Kp)/1000000000
    print("\\begin{tabular}{c | c  c  c  c  c }\hline \hline", file=f)
    print("          & $K$        & $A$             & $B$     & $E$        & $G$            \\\ \hline", file=f)
    print("Paper     & 5727.12$\pm x.xx$ & 21.54$\pm 0.37$ & 23.91$\pm 0.31$ & 30.34 $\pm 0.23$ & 54.37$\pm 0.58$ \\\ ", file=f)
    print("Sampled   &",round(Kp,1),' $\\pm',round(delta_Kp,2),'$ &', round(A,2),' $\\pm',round(delta_A,2),'$ &', round(B,2),' $\\pm',round(delta_B,2),'$ &', round(E,2), ' $\\pm',round(delta_E,2),'$ &', round(G,2), ' $\\pm',round(delta_G,2),'$ \\\ ', file=f)
    print("L.Algebra &",round(LA[0][0],1),' &', round(LA[1][0],2),' &', round(LA[2][0],2),' &', round(LA[3][0],2), '&', round(LA[4][0],2), ' \\\ ', file=f)
    print("\hline\hline", file=f)
    print("\end{tabular}", file=f)
    print("\caption{Model pararameters in MeV, for states: $", name, "$}",file=f)
    print("    ")
    f.close()
    #print(Kp/5727.12, A/21.54, B/23.91, E/30.34,G/54.37)
    
    f = open('./tables/masses_'+name+'.tex', "w")
    print("\\begin{tabular}{c | c  c  c  c  c}\hline \hline", file=f)
    print("Mass State & Experiment  &   Predicted mass  &    Predicted mass & diff pred & diff sampl\\\ ", file=f)
    print("           & (MeV)       &   old (MeV)       &    sampled (MeV)  &     (\\%) &     (\\%)  \\\ \hline", file=f)
    # print("Mass State & Experiment  &   Predicted mass  &    Predicted mass \\\  ")
    # print("           & (MeV)       &   old (MeV)       &    sampled (MeV) \\\ \hline  ")

    tot_diff_pred, tot_diff_sample = 0,0
    for i in range(len(v_param)):
        if bootstrap:
            mass = bootstrap_masses[i]
            if not asymmetric:
                error= bootstrap_errors[i]
            else:
                error_up = delta_up[i]
                error_dn = delta_dn[i]
        else:
            mass = sum_mass[i] + Kp*v_param[i] + A*w_param[i] + B*x_param[i] + E*y_param[i] + G*z_param[i]
            error = np.sqrt( (delta_Kp*v_param[i])**2 + (delta_A*w_param[i])**2 + (delta_B*x_param[i])**2 + (delta_E*y_param[i])**2 + (delta_G*z_param[i])**2 )
            
        color = "red"
        old_exp = exp[i]-old[i]
        new_exp = exp[i]-mass
        tot_diff_pred+=abs(old_exp)
        tot_diff_sample+=abs(new_exp)
        old_exp_abs = np.abs(1 - old[i] / exp[i])
        new_exp_abs = np.abs(1 - mass / exp[i])        
        if(old_exp_abs > new_exp_abs):
            color = "blue"

        if not asymmetric:
            print(quantum[i], exp[i], '$\\pm', delta_exp[i],'$  &', old[i], '$\\pm xx$  &', '\\textcolor{'+color+'}{', round(mass,1),' $\\pm',round(error,1),'$}  &  ',
                  round(old_exp,1), '(', round(old_exp_abs*100,1),')  &  ', round(new_exp,1), '(',round(new_exp_abs*100,1),') \\\  ', file=f)
        else:
            print(quantum[i], exp[i], '$\\pm', delta_exp[i],'$  &', old[i], '$\\pm xx$  &', '\\textcolor{'+color+'}{', round(mass,1),' $^{+',round(error_up,1),'}_{',round(error_dn,1),'}$}  &  ',
                  round(old_exp,1), '(', round(old_exp_abs*100,1),')  &  ', round(new_exp,1), '(',round(new_exp_abs*100,1),') \\\  ', file=f)

        
    print('\hline', file=f)
    print("  &  &  & Total diff & ", round(tot_diff_pred), " &" , round(tot_diff_sample,1), "\\\ ", file=f)
    print('\hline \hline', file=f)
    print('\end{tabular}', file=f)
    print("\caption{Every quantity is in MeV, except for percentage differences. States:", name, "}",file=f)


def sampled_prediction(sum_mass, v_param, w_param, x_param, y_param, z_param,
                       sampled_k, sampled_a, sampled_b, sampled_e, sampled_g,
                       rho_ak,rho_bk,rho_ba,rho_ek,rho_ea,rho_eb,rho_gk,rho_ga,rho_gb,rho_ge, name):

    bootstrap_masses,sorted_masses,symm_errors = [],[],[]

    for i in range(len(sum_mass)):
        dummy = ([])
        for j in range(len(sampled_k)):
            mass = sum_mass[i] + sampled_k[j]*v_param[i] + sampled_a[j]*w_param[i] + sampled_b[j]*x_param[i] + sampled_e[j]*y_param[i] + sampled_g[j]*z_param[i]
            dummy = np.append(dummy,mass)
        bootstrap_masses.append(dummy.mean())
        symm_errors.append(dummy.std(ddof=1))
        sorted_masses.append(np.sort(dummy))
                
    bootstrap_masses = np.array(bootstrap_masses)
    symm_errors      = np.array(symm_errors)
    sorted_masses    = np.array(sorted_masses)

    # asymmetric error calculation via 68% quantile method
    N = len(sorted_masses[0])
    quantile_dn = int(np.floor(N*0.1587))
    quantile_up = int(np.floor(N*0.8413))
    asymmetric_up, asymmetric_dn = ([]), ([])
    for i in range(len(sum_mass)):
        asymmetric_up = np.append(asymmetric_up, sorted_masses[i][quantile_up-1] - np.mean(sorted_masses[i]))
        asymmetric_dn = np.append(asymmetric_dn, sorted_masses[i][quantile_dn-1] - np.mean(sorted_masses[i]))

        
    dummy = ([])
    for i in range(len(sampled_k)):
        mass = sum_mass[0] + sampled_k[i]*v_param[0] + sampled_a[i]*w_param[0] + sampled_b[i]*x_param[0] + sampled_e[i]*y_param[0] + sampled_g[i]*z_param[0]
        dummy = np.append(dummy,mass)
        
    plot(dummy,'mass','mass','charm', name)
    print(dummy.mean(), dummy.std(ddof=1), len(dummy))
    dummy = np.sort(dummy)

    fK = v_param[0]
    fA = w_param[0]
    fB = x_param[0]
    fE = y_param[0]
    fG = z_param[0]

    rho_ak = np.mean(rho_ak)
    rho_bk = np.mean(rho_bk)
    rho_ba = np.mean(rho_ba)
    rho_ek = np.mean(rho_ek)
    rho_ea = np.mean(rho_ea)
    rho_eb = np.mean(rho_eb)
    rho_gk = np.mean(rho_gk)
    rho_ga = np.mean(rho_ga)
    rho_gb = np.mean(rho_gb)
    rho_ge = np.mean(rho_ge)

    delta_K = np.std(sampled_k,ddof=1)
    delta_A = np.std(sampled_a,ddof=1)
    delta_B = np.std(sampled_b,ddof=1)
    delta_E = np.std(sampled_e,ddof=1)
    delta_G = np.std(sampled_g,ddof=1)
        
    error_diag  =  (delta_K*fK)**2 + (delta_A*fA)**2 + (delta_B*fB)**2 + (delta_E*fE)**2 + (delta_G*fG)**2
    error_off1  =  fA*fK*rho_ak*delta_A*delta_K + fB*fK*rho_bk*delta_B*delta_K + fB*fA*rho_ba*delta_B*delta_A + fE*fK*rho_ek*delta_E*delta_K
    error_off2  =  fE*fA*rho_ea*delta_E*delta_A + fE*fB*rho_eb*delta_E*delta_B + fG*fK*rho_gk*delta_G*delta_K + fG*fA*rho_ga*delta_G*delta_A
    error_off3  =  fG*fB*rho_gb*delta_G*delta_B + fG*fE*rho_ge*delta_G*delta_E 

    error_total =  np.sqrt( error_diag + 2*(error_off1+error_off2+error_off3))

    print(error_total)
    return bootstrap_masses, symm_errors, asymmetric_up, asymmetric_dn

def linear_algebra(name):
    A = octave.load('./octave/matrix_'+ name + '.mat')
    b = octave.load('./octave/vector_'+ name + '.mat')
    return octave.linsolve(A,b)


def normal_test(sample,alpha,verbose):
    # # hypothesis test: null hypothesis, the data is gaussian distributed

    # # Shapiro-Wilk
    # stat, p = shapiro(sample)
    # if verbose:
    #     if p > alpha: print('Shapiro this is Gaussian', p)
    #     else:         print('Shapiro this is NOT Gaussian', p)

    # # chisquare
    # stat, p = chisquare(sample)
    # if verbose:
    #     if p > alpha: print('Chisquare this is Gaussian', p)
    #     else:         print('Chisquare this is NOT Gaussian', p)

    # # lilliefors
    # stat, p = lilliefors(sample)
    # if verbose:
    #     if p > alpha: print('Lilliefors this is Gaussian', p)
    #     else:         print('Lilliefors this is NOT Gaussian', p)

    # # kolmogorov
    # stat, p = kstest(sample, 'norm')
    # if verbose:
    #     if p > alpha: print('Kolmogorov this is Gaussian', p)
    #     else:         print('Kolmogorov this is NOT Gaussian', p)

    # Angostino
    k2, p = normaltest(sample)
    if verbose:
        if p > alpha: print('Angostino this is Gaussian', p)
        else:         print('Angostino this is NOT Gaussian', p)    

    return p,alpha


def values(name):
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

    if name=='sigma_lamb':        
        quantum = ['$\\vert uuc,1/2,1/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,3/2,3/2,0,1,10/3 \\rangle$ &',   '$\\vert uuc,1/2,1/2,1,1,10/3 \\rangle$ &',
                   '$\\vert udc,1/2,1/2,0,0,4/3  \\rangle$ &',   '$\\vert udc,1/2,1/2,1,0,10/3 \\rangle$ &',   '$\\vert udc,3/2,1/2,1,0,4/3  \\rangle$ &']
        old = np.array([2453.1, 2517.7, 2819.0, 2283.7, 2649.7, 2685.6])
        exp = np.array([2453.9, 2518.0, 2801.0, 2286.5, 2592.3, 2625.0])
        delta_exp = np.array([0.1,2.3,6.0,0.1,0.4,0.2])
        
    return quantum, old, exp, delta_exp


def correlation_matrix(rho_ak,rho_bk,rho_ba,rho_ek,rho_ea,rho_eb,rho_gk,rho_ga,rho_gb,rho_ge,name):

    rho_ak=round(np.mean(rho_ak),2)
    rho_bk=round(np.mean(rho_bk),2)
    rho_ba=round(np.mean(rho_ba),2)
    rho_ek=round(np.mean(rho_ek),2)
    rho_ea=round(np.mean(rho_ea),2)
    rho_eb=round(np.mean(rho_eb),2)
    rho_gk=round(np.mean(rho_gk),2)
    rho_ga=round(np.mean(rho_ga),2)
    rho_gb=round(np.mean(rho_gb),2)
    rho_ge=round(np.mean(rho_ge),2)

    
    f = open('./tables/correlation_'+name+'.tex', "w")
    print("\\begin{tabular}{c  c  c  c  c  c}\hline \hline", file=f)
    print("     &  $K$      &     $A$   &      $B$   &      $E$  & $G$ \\\ \hline", file=f)
    print(" $K$ &     1     &           &            &           &   \\\ ", file=f)
    print(" $A$ &",rho_ak, "&      1    &            &           &   \\\ ", file=f)
    print(" $B$ &",rho_bk, "&", rho_ba,"&      1     &           &   \\\ ", file=f)
    print(" $E$ &",rho_ek, "&", rho_ea,"&", rho_eb, "&      1    &   \\\ ", file=f)
    print(" $G$ &",rho_gk, "&", rho_ga,"&", rho_gb, "&", rho_ge,"& 1 \\\ \hline \hline", file=f)
    print('\end{tabular}', file=f)
