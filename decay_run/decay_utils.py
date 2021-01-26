'''
---------------------------------------------------------------
 Code to calcualte charm-baryon decay widths
 Authors: A. Ramirez-Morales and H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''

import numpy as np


def reduced_masses(baryons):
    m_lam_omega = 864.97
    m_lam_casca = 763.228
    m_lam_sigma = 647.118
    m_rho_omega = 450.
    m_rho_casca = 372.5
    m_rho_sigma = 295.    
    # fetch reduced masses of the harmonic oscillator
    if(baryons==1):
        return m_lam_omega,   m_rho_omega        
    elif(baryons==2 or baryons ==5):
        return m_lam_casca,   m_rho_casca            
    elif(baryons==3 or baryons==4):
        return m_lam_sigma,   m_rho_sigma
            

def alphas(k_prim, m_lam_rho):
    value1 = (np.sqrt(3./m_lam_rho)) * k_prim
    value2 = value1*m_lam_rho
    return np.sqrt(value2)/1000. # transform from GeV -> MeV


def append_dic(baryons, state1, state2,state3,
               state4, state5, state6, state7):
    baryons.append(state1)
    baryons.append(state2)
    baryons.append(state3)
    baryons.append(state4)
    baryons.append(state5)
    baryons.append(state6)
    baryons.append(state7)
    return baryons

def state_labels(baryon, ModEx, decPr):
    
    if(baryon==1):
        baryon_name = "omega"
        if(decPr==1):
            decPr_name = "Xi+K"
        elif(decPr==2):
            decPr_name = "Xi'+K"
        elif(decPr==3):
            decPr_name = "Xi*+K"
    elif(baryon==2 or baryon==5):
        if(baryon==2): baryon_name = 'cas_6'
        if(baryon==5): baryon_name = 'cas_3'        
        if(decPr==1):
            decPr_name = "lamb+K"
        elif(decPr==2):
            decPr_name = "Xi+Pi"
        elif(decPr==3):
            decPr_name = "Xi'+Pi"
        elif(decPr==4):
            decPr_name = "Xi*+Pi"
        elif(decPr==5):
            decPr_name = "Sig+K"
        elif(decPr==6):
            decPr_name = "Sig*+K"
        elif(decPr==7):
            decPr_name = "Xi+eta"
    elif(baryon==3):
        baryon_name = 'sigma'
        if(decPr==1):
            decPr_name = "Sig+Pi"
        elif(decPr==2):
            decPr_name = "Sig*+Pi"
        elif(decPr==3):
            decPr_name = "lamb+Pi"
        elif(decPr==4):
            decPr_name = "Sig+Eta"
        elif(decPr==5):
            decPr_name = "Xi+K"    
    elif(baryon==4):
        baryon_name = 'lamda'
        if(decPr==1):
            decPr_name = "Sig+Pi"
        elif(decPr==2):
            decPr_name = "Sig*+Pi"
        elif(decPr==3):
            decPr_name = "lamb+eta"


    if(ModEx==1):  ModEx_name = 'lam'
    elif(ModEx==2): ModEx_name = 'rho'

    return baryon_name,ModEx_name,decPr_name

def latex_decay_label(baryon,decPr):

    if(baryon==1):
        baryon_name = "omega"
        if(decPr==1):
            decPr_name = "$\Gamma(\Xi_{c} K)$"
        elif(decPr==2):
            decPr_name = "$\Gamma(\Xi'_{c} K)$"
        elif(decPr==3):
            decPr_name = "$\Gamma(\Xi^{*}_{c} K)$"
    elif(baryon==2 or baryon==5):
        if(baryon==2): baryon_name = 'cas_6'
        if(baryon==5): baryon_name = 'cas_3'        
        if(decPr==1):
            decPr_name = "$\Gamma(\Lambda_{c} K)$"
        elif(decPr==2):
            decPr_name = "$\Gamma(\Xi_{c} \pi)$"
        elif(decPr==3):
            decPr_name = "$\Gamma(\Xi'_{c} \pi)$"
        elif(decPr==4):
            decPr_name = "$\Gamma(\Xi^{*}_{c} \pi)$"
        elif(decPr==5):
            decPr_name = "$\Gamma(\Sigma_{c} K)$"
        elif(decPr==6):
            decPr_name = "$\Gamma(\Sigma^{*}_{c} K)$"
        elif(decPr==7):
            decPr_name = "$\Gamma(\Xi_{c} \eta)$"
    elif(baryon==3):
        baryon_name = 'sigma'
        if(decPr==1):
            decPr_name = "$\Gamma(\Sigma_{c} \pi)$"
        elif(decPr==2):
            decPr_name = "$\Gamma(\Sigma^{*}_{c} \pi)$"
        elif(decPr==3):
            decPr_name = "$\Gamma(\Lambda_{c} \pi)$"
        elif(decPr==4):
            decPr_name = "$\Gamma(\Sigma_{c} \eta)$"
        elif(decPr==5):
            decPr_name = "$\Gamma(\Xi_{c} K)$"    
    elif(baryon==4):
        baryon_name = 'lamda'
        if(decPr==1):
            decPr_name = "$\Gamma(\Sigma_{c} \pi)$"
        elif(decPr==2):
            decPr_name = "$\Gamma(\Sigma^{*}_{c} \pi)$"
        elif(decPr==3):
            decPr_name = "$\Gamma(\Lambda_{c} \eta)$"

    return decPr_name

def print_row_latex(state_name,state_decays, f_out):
    nstate=len(state_decays)

    print("state",state_name, "  &  ", end='',file=f_out)
    for i in range(nstate):
        value=0
        if(state_decays[i]==0.0): value = '$\\dagger$'
        else: value = round(state_decays[i],4)    
        print(value,"  &  ", end='', file=f_out)
        
    print(round(np.sum(state_decays),2), '\\\\', file=f_out)

def print_header_latex(name_states, f_out):
    nNames = len(name_states)    
    print("\\begin{tabular}{c |", end='',file=f_out)
    for i in range(nNames-1): print("  c", end='',file=f_out)
    print("} \hline \hline", file=f_out)
    for i in range(nNames-1): print(name_states[i]," & ", end='',file=f_out)
    print(name_states[nNames-1]," \\\\ \hline", file=f_out)

def print_bottom_latex(baryons,f_decay):
    print('\hline \hline', file=f_decay)
    print('\end{tabular}', file=f_decay)
    print("\caption{Decay widths in MeV, for states: ", baryons,"}",file=f_decay)
    print("\label{tab:gordo}", file=f_decay)
