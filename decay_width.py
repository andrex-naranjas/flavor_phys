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

class DecayWidths:

    def __init__(self, bootstrap=False, baryons=''):
        self.fetch_decay_masses(bootstrap)
        self.set_m_lam_rho(bootstrap)

        
    def total_decay_width(self, bootstrap, baryons, k_prim, massA, SA_val, L_val, JA_val, SL_val, ModEx_val):
        MassA = massA/1000.0
        SA_qm = SA_val
        LA_qm = L_val
        JA_qm = JA_val
        SL_qm = SL_val
        baryon= self.baryon_flag(baryons)
        ModEx = self.ModEx_flag(ModEx_val)
        channel_widths = ([])
        nChannels = self.n_channels(baryons)
        m_lam, m_rho = self.reduced_masses(bootstrap, baryons)
        
        alpha_lam = self.alphas(k_prim, m_lam)
        alpha_rho = self.alphas(k_prim, m_rho)

        for i in range(nChannels):        
            decPr = i+1
            MassB,MassC = self.decay_masses(bootstrap, baryons, decPr)            
            single_decay_value = width.decay_width(MassA, MassB, MassC, SA_qm,
                                                   LA_qm, JA_qm, SL_qm, alpha_lam, alpha_rho,
                                                   baryon, ModEx, decPr)
            channel_widths = np.append(channel_widths, single_decay_value)
            
        # sum the individual width to obtain total width
        total_decay_width = np.sum(channel_widths)
        #print(total_decay_width, 'Total width', baryons, SA_qm, JA_qm)
        
        return total_decay_width
    

    def baryon_flag(self, baryons):
        #omegas,cascades,sigmas,lambdas,cascades_anti3
        if(baryons=='omegas'):           return 1    
        elif(baryons=='cascades'):       return 2
        elif(baryons=='sigmas'):         return 3
        elif(baryons=='lambdas'):        return 4
        elif(baryons=='cascades_anti3'): return 5
        

    def ModEx_flag(self, ModEx_val):
        # lam =1 , rho=2
        if(ModEx_val=='lam'):   return 1
        elif(ModEx_val=='rho'): return 2


    def n_channels(self, baryons):
        # set number of decay channels has each baryon
        # omegas,cascades,sigmas,lambdas,cascades_anti3
        if(baryons=='omegas'):           return 3
        elif(baryons=='cascades'):       return 7
        elif(baryons=='sigmas'):         return 5
        elif(baryons=='lambdas'):        return 3
        elif(baryons=='cascades_anti3'): return 7
        

    def alphas(self, k_prim, m_lam_rho):
        value1 = (np.sqrt(3./m_lam_rho)) * k_prim
        value2 = value1*m_lam_rho
        return np.sqrt(value2)/1000. # transform from GeV -> MeV

            
    def decay_masses(self, bootstrap, baryons, decPr):
        # fetch mass of the decay products        
        if(baryons=='omegas'):
            if(decPr==1):
                if not bootstrap: return self.xi_mass,   self.kaon_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==2):                
                if not bootstrap: return self.xi_p_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_s_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_kaon, size=None)                
        elif(baryons=='cascades'):
            if(decPr==1):
                if not bootstrap: return self.lambda_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==2):
                if not bootstrap: return self.xi_mass,     self.pion_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_p_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_s_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==5):
                if not bootstrap: return self.sigma_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==6):
                if not bootstrap: return self.sigma_s_mass,self.kaon_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==7):
                if not bootstrap: return self.xi_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_eta, size=None)
        elif(baryons=='sigmas'):
            if(decPr==1):
                if not bootstrap: return self.sigma_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==2):
                if not bootstrap: return self.sigma_s_mass, self.pion_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.lambda_mass,  self.pion_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==4):
                if not bootstrap: return self.sigma_mass,   self.eta_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_eta, size=None)
            elif(decPr==5):
                if not bootstrap: return self.xi_mass,      self.kaon_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_kaon, size=None)
        elif(baryons=='lambdas'):
            if(decPr==1):
                if not bootstrap: return self.sigma_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==2):
                if not bootstrap: return self.sigma_s_mass, self.pion_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.lambda_mass,  self.eta_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_eta, size=None)
        elif(baryons=='cascades_anti3'):
            if(decPr==1):
                if not bootstrap: return self.lambda_mass, self.kaon_mass
                else: return np.random.choice(self.gauss_lambda, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==2):
                if not bootstrap: return self.xi_mass,     self.pion_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==3):
                if not bootstrap: return self.xi_p_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_p, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==4):
                if not bootstrap: return self.xi_s_mass,   self.pion_mass
                else: return np.random.choice(self.gauss_xi_s, size=None), np.random.choice(self.gauss_pion, size=None)
            elif(decPr==5):
                if not bootstrap: return self.sigma_mass,  self.kaon_mass
                else: return np.random.choice(self.gauss_sigma, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==6):
                if not bootstrap: return self.sigma_s_mass,self.kaon_mass
                else: return np.random.choice(self.gauss_sigma_s, size=None), np.random.choice(self.gauss_kaon, size=None)
            elif(decPr==7):
                if not bootstrap: return self.xi_mass,     self.eta_mass
                else: return np.random.choice(self.gauss_xi, size=None), np.random.choice(self.gauss_eta, size=None)
                

    def reduced_masses(self, bootstrap, baryons):
        # fetch reduced masses of the harmonic oscillator
        if(baryons=='omegas'):
            if not bootstrap: return self.m_lam_omega, self.m_rho_omega            
            else: return np.random.choice(self.gauss_m_lam_omega, size=None), np.random.choice(self.gauss_m_rho_omega, size=None)
        elif(baryons=='cascades' or baryons =='cascades_anti3'):
            if not bootstrap: return self.m_lam_casca, self.m_rho_casca
            else: return np.random.choice(self.gauss_m_lam_casca, size=None), np.random.choice(self.gauss_m_rho_casca, size=None)
        elif(baryons=='sigmas' or baryons=='lambdas'):
            if not bootstrap: return self.m_lam_sigma, self.m_rho_sigma
            else: return np.random.choice(self.gauss_m_lam_sigma, size=None), np.random.choice(self.gauss_m_rho_sigma, size=None)


    def fetch_decay_masses(self, bootstrap):
        self.pion_mass   = 0.140
        self.kaon_mass   = 0.493
        self.lambda_mass = 2.286
        self.xi_mass     = 2.469
        self.xi_p_mass   = 2.578
        self.xi_s_mass   = 2.645
        self.sigma_mass  = 2.455
        self.sigma_s_mass= 2.518 # 2.520
        self.eta_mass    = 0.548

        if(bootstrap):
            self.gauss_pion   = np.random.normal(0.140, 0.01, 10000)
            self.gauss_kaon   = np.random.normal(0.493, 0.01, 10000)
            self.gauss_lambda = np.random.normal(2.286, 0.01, 10000)
            self.gauss_xi     = np.random.normal(2.469, 0.01, 10000)
            self.gauss_xi_p   = np.random.normal(2.578, 0.01, 10000)
            self.gauss_xi_s   = np.random.normal(2.645, 0.01, 10000)
            self.gauss_sigma  = np.random.normal(2.455, 0.01, 10000)
            self.gauss_sigma_s= np.random.normal(2.518, 0.01, 10000)
            self.gauss_eta    = np.random.normal(0.548, 0.01, 10000)
            

    def set_m_lam_rho(self, bootstrap):
        self.m_lam_omega = 864.97
        self.m_lam_casca = 763.228
        self.m_lam_sigma = 647.118
        self.m_rho_omega = 450.
        self.m_rho_casca = 372.5
        self.m_rho_sigma = 295.

        if(bootstrap):
            self.gauss_m_lam_omega = np.random.normal(864.97 , 10, 10000)
            self.gauss_m_lam_casca = np.random.normal(763.228, 10, 10000)
            self.gauss_m_lam_sigma = np.random.normal(647.118, 10, 10000)
            self.gauss_m_rho_omega = np.random.normal(450.   , 10, 10000)
            self.gauss_m_rho_casca = np.random.normal(372.5  , 10, 10000)
            self.gauss_m_rho_sigma = np.random.normal(295.   , 10, 10000)

            
        
    # def single_decay_width(self, baryons, mass, SA_val, L_val, JA_val, SL_val, ModEx_val):
    #     #state1 = {'MA':2.797,'MB':2.286,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':5,'ModEx':1,'decPr':1}
    #     MassA = mass/1000.0
    #     MassB,MassC = self.decay_masses(baryons, decPr=3)
    #     SA_qm = SA_val
    #     LA_qm = L_val
    #     JA_qm = JA_val
    #     SL_qm = SL_val
    #     baryon= self.baryon_flag(baryons)
    #     ModEx = self.ModEx_flag(ModEx_val)
    #     decPr = 3    
    #     decay_value = width.decay_width(MassA, MassB, MassC, SA_qm,
    #                                 LA_qm, JA_qm, SL_qm, baryon, ModEx, decPr)
    #     return single_decay_value
    # # baryon_name, ModEx_name, decPr_name = du.state_labels(baryon,ModEx,decPr)
    # # print('%6s |  %4s | %7s |  %5.3f |  %5.3f | %5.3f |  %5.1f |  %5.1f |  %5.1f |  %5.1f | %5.6f '
    # #       %(baryon_name, ModEx_name, decPr_name, MassA, MassB, MassC, JA_qm, LA_qm, SA_qm, SL_qm,  decay_value))
    # # if(i%7==6 and i !=0):
    # #     print('--------------------------------------------------------------------------------------------------')
