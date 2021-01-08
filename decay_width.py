'''
---------------------------------------------------------------
 Code to calculate charm-baryon decay widths
 Authors: A. Ramirez-Morales and H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''
from charm_wrapper import decay
import decay_utils as du
import numpy as np

width  = decay()

# baryon FLAG: 1 -> omega, 2->cascade_6, 3->sigma,# 4 -> lambda, 5-> cascade_3 
# ModEx  FLAG: 1 -> lambda, 2->rho  Excitation
# decPr  FLAG: 3->Xi'+Pi, 5->Sigma+K  decayProduct Flag

def baryon_flag(baryons):
    #omegas,cascades,sigmas,lambdas,cascades_anti3
    if(baryons=='omegas'):           return 1    
    elif(baryons=='cascades'):       return 2
    elif(baryons=='sigmas'):         return 3
    elif(baryons=='lambdas'):        return 4
    elif(baryons=='cascades_anti3'): return 5

def ModEx_flag(ModEx_val):
    # lam =1 , rho=2
    if(ModEx_val=='lam'):   return 1
    elif(ModEx_val=='rho'): return 2

def decay_masses(baryons, decPr):
    # mass of the decay products
    pion_mass   = 0.140
    kaon_mass   = 0.493
    lambda_mass = 2.286
    xi_mass     = 2.469
    xi_p_mass   = 2.578
    xi_s_mass   = 2.645
    sigma_mass  = 2.455
    sigma_s_mass= 2.518 # 2.520
    eta_mass    = 0.548    

    if(baryons=='omegas'):
        if(decPr==1):   return xi_mass,   kaon_mass
        elif(decPr==2): return xi_p_mass, kaon_mass
        elif(decPr==3): return xi_s_mass, kaon_mass
    elif(baryons=='cascades'):
        if(decPr==1):   return lambda_mass, kaon_mass
        elif(decPr==2): return xi_mass,     pion_mass
        elif(decPr==3): return xi_p_mass,   pion_mass
        elif(decPr==4): return xi_s_mass,   pion_mass
        elif(decPr==5): return sigma_mass,  kaon_mass
        elif(decPr==6): return sigma_s_mass,kaon_mass
        elif(decPr==7): return xi_mass,     eta_mass        
    elif(baryons=='sigmas'):
        if(decPr==1):   return sigma_mass,   pion_mass
        elif(decPr==2): return sigma_s_mass, pion_mass
        elif(decPr==3): return lambda_mass,  pion_mass
        elif(decPr==4): return sigma_mass,   eta_mass
        elif(decPr==5): return xi_mass,      kaon_mass
    elif(baryons=='lambdas'):
        if(decPr==1):   return sigma_mass,   pion_mass
        elif(decPr==2): return sigma_s_mass, pion_mass
        elif(decPr==3): return lambda_mass,  eta_mass
    elif(baryons=='cascades_anti3'):
        if(decPr==1):   return lambda_mass, kaon_mass
        elif(decPr==2): return xi_mass,     pion_mass
        elif(decPr==3): return xi_p_mass,   pion_mass
        elif(decPr==4): return xi_s_mass,   pion_mass
        elif(decPr==5): return sigma_mass,  kaon_mass
        elif(decPr==6): return sigma_s_mass,kaon_mass
        elif(decPr==7): return xi_mass,     eta_mass

def n_channels(baryons):
    #omegas,cascades,sigmas,lambdas,cascades_anti3
    if(baryons=='omegas'):           return 3    
    elif(baryons=='cascades'):       return 7
    elif(baryons=='sigmas'):         return 5
    elif(baryons=='lambdas'):        return 3
    elif(baryons=='cascades_anti3'): return 7
    
        
def total_decay_width(baryons, mass, SA_val, L_val, JA_val, SL_val, ModEx_val):    
    MassA = mass/1000.0
    SA_qm = SA_val
    LA_qm = L_val
    JA_qm = JA_val
    SL_qm = SL_val
    baryon= baryon_flag(baryons)
    ModEx = ModEx_flag(ModEx_val)
    channel_widths = ([])
    nChannels = n_channels(baryons)

    for i in range(nChannels):        
        decPr = i+1
        MassB,MassC = decay_masses(baryons, decPr)
        single_decay_value = width.decay_width(MassA, MassB, MassC, SA_qm,
                                               LA_qm, JA_qm, SL_qm, baryon, ModEx, decPr)        
        channel_widths = np.append(channel_widths, single_decay_value)

    total_decay_width = np.sum(channel_widths)
    return total_decay_width


def single_decay_width(baryons, mass, SA_val, L_val, JA_val, SL_val, ModEx_val):
    #state1 = {'MA':2.797,'MB':2.286,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':5,'ModEx':1,'decPr':1}
    MassA = mass/1000.0
    MassB,MassC = decay_masses(baryons, decPr=3)
    SA_qm = SA_val
    LA_qm = L_val
    JA_qm = JA_val
    SL_qm = SL_val
    baryon= baryon_flag(baryons)
    ModEx = ModEx_flag(ModEx_val)
    decPr = 3    
    decay_value = width.decay_width(MassA, MassB, MassC, SA_qm,
                                    LA_qm, JA_qm, SL_qm, baryon, ModEx, decPr)
    return single_decay_value


    # baryon_name, ModEx_name, decPr_name = du.state_labels(baryon,ModEx,decPr)
    # print('%6s |  %4s | %7s |  %5.3f |  %5.3f | %5.3f |  %5.1f |  %5.1f |  %5.1f |  %5.1f | %5.6f '
    #       %(baryon_name, ModEx_name, decPr_name, MassA, MassB, MassC, JA_qm, LA_qm, SA_qm, SL_qm,  decay_value))
    # if(i%7==6 and i !=0):
    #     print('--------------------------------------------------------------------------------------------------')

