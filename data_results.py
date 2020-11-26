#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# code to obtain uncertainties of quarkonium mass spectrum
# author: A. Ramirez-Morales (andres.ramirez.morales@cern.ch)

# data visalization module

import data_visualization as dv
import data_utils as du
import numpy as np

# Results class

class CharmResults:
    
    def __init__(self, params, sampled, corr_mat, bootstrap, asymmetric, name):
        self.params = params
        self.sampled = sampled
        self.corr_mat = corr_mat
        self.bootstrap = bootstrap
        self.asymmetric = asymmetric
        self.name = name
        self.fetch_values()

    
    def model_mass(self, i, j, sampled=False):
        # compute baryon mass according to model
        
        if not sampled: # here the mass is calculated using the average value of parameters
            return self.sum_mass[i] + self.Kp*self.v_param[i] + self.A*self.w_param[i] \
                + self.B*self.x_param[i] + self.E*self.y_param[i] + self.G*self.z_param[i]
        else: # here the mass is calculated using every single simulated parameter
            return self.sum_mass[i] + self.sampled_k[j]*self.v_param[i] + self.sampled_a[j]*self.w_param[i] \
                + self.sampled_b[j]*self.x_param[i] + self.sampled_e[j]*self.y_param[i] + self.sampled_g[j]*self.z_param[i]

                
    def mass_prediction(self):
        # mass prediction (compare and make tables)
        
        if self.bootstrap:
            bootstrap_masses, bootstrap_errors, delta_up, delta_dn = self.sampled_prediction()
            
        quantum, old, exp, delta_exp = du.names_values(self.name)
            
        f_note  = open('./tables/masses_'+self.name+'_note.tex', "w")
        f_paper = open('./tables/masses_'+self.name+'_paper.tex', "w")
        self.latex_header(table_file=f_note,  paper=False) # write the header of the latex table
        self.latex_header(table_file=f_paper, paper=True)  # write the header of the latex table
        
        tot_diff_pred, tot_diff_sample = 0,0
        for i in range(len(self.sum_mass)): # run over exp mass states
            if self.bootstrap:
                mass = bootstrap_masses[i]
                if not self.asymmetric:
                    error = bootstrap_errors[i]
                else:
                    error_up = delta_up[i]
                    error_dn = delta_dn[i]
            else:                
                mass  = self.model_mass(i, 0, sampled=False)
                error = self.analytical_error(i)
            
            color = "red"
            old_exp = exp[i]-old[i]
            new_exp = exp[i]-mass
            tot_diff_pred+=abs(old_exp)
            tot_diff_sample+=abs(new_exp)
            old_exp_abs = np.abs(1 - old[i] / exp[i])
            new_exp_abs = np.abs(1 - mass / exp[i])        
            if(old_exp_abs > new_exp_abs):
                color = "blue"

            if not self.asymmetric:
                print(quantum[i], exp[i], '$\\pm', delta_exp[i],'$  &', old[i], '$\\pm xx$  &', '\\textcolor{'+color+'}{', round(mass,1),' $\\pm',round(error,1),'$}  &  ',
                      round(old_exp,1), '(', round(old_exp_abs*100,1),')  &  ', round(new_exp,1), '(',round(new_exp_abs*100,1),') \\\  ', file=f_note)
                print(quantum[i], exp[i], '$\\pm', delta_exp[i],'$  &',  round(mass,1),' $\\pm',round(error,1),'\\\ ', file=f_paper)
            else:
                print(quantum[i], exp[i], '$\\pm', delta_exp[i],'$  &', old[i], '$\\pm xx$  &', '\\textcolor{'+color+'}{', round(mass,1),' $^{+',round(error_up,1),'}_{',round(error_dn,1),'}$}  &  ',
                      round(old_exp,1), '(', round(old_exp_abs*100,1),')  &  ', round(new_exp,1), '(',round(new_exp_abs*100,1),') \\\  ', file=f_note)
                print(quantum[i], exp[i], '$\\pm', delta_exp[i],'$  &', round(mass,1),' $^{+',round(error_up,1),'}_{',round(error_dn,1),'}$ \\\ ', file=f_paper)

        self.latex_bottom(f_note,tot_diff_pred,tot_diff_sample,paper=False) # write bottom's table
        self.latex_bottom(f_paper,0,0,paper=True) # write bottom's table

                
    def sampled_prediction(self):

        bootstrap_masses,sorted_masses,symm_errors = [],[],[]

        for i in range(len(self.sum_mass)):
            dummy = ([])
            for j in range(len(self.sampled_k)):
                mass = self.model_mass(i, j, sampled=True)                
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
        for i in range(len(self.sum_mass)):
            asymmetric_up = np.append(asymmetric_up, sorted_masses[i][quantile_up-1] - np.mean(sorted_masses[i]))
            asymmetric_dn = np.append(asymmetric_dn, sorted_masses[i][quantile_dn-1] - np.mean(sorted_masses[i]))
        
        dummy = ([])
        for j in range(len(self.sampled_k)):
            mass = self.model_mass(0, j, sampled=True)            
            dummy = np.append(dummy,mass)
        
        dv.plot(dummy,'mass','mass','charm', self.name)
        print(dummy.mean(), dummy.std(ddof=1), len(dummy))
        dummy = np.sort(dummy)
        error_total = self.analytical_error(0) # error for the first mass

        print(error_total)
        return bootstrap_masses, symm_errors, asymmetric_up, asymmetric_dn
        

    def analytical_error(self, iter_a):
        # propagate the error using the analytical method,
        # and the covariance matrix    
        fK, delta_K = self.v_param[iter_a], self.delta_Kp
        fA, delta_A = self.w_param[iter_a], self.delta_A 
        fB, delta_B = self.x_param[iter_a], self.delta_B 
        fE, delta_E = self.y_param[iter_a], self.delta_E 
        fG, delta_G = self.z_param[iter_a], self.delta_G 
                
        error_diag  =  (delta_K*fK)**2 + (delta_A*fA)**2 + (delta_B*fB)**2 + (delta_E*fE)**2 + (delta_G*fG)**2
        error_off1  =  fA*fK*self.rho_ak*delta_A*delta_K + fB*fK*self.rho_bk*delta_B*delta_K + fB*fA*self.rho_ba*delta_B*delta_A + fE*fK*self.rho_ek*delta_E*delta_K
        error_off2  =  fE*fA*self.rho_ea*delta_E*delta_A + fE*fB*self.rho_eb*delta_E*delta_B + fG*fK*self.rho_gk*delta_G*delta_K + fG*fA*self.rho_ga*delta_G*delta_A
        error_off3  =  fG*fB*self.rho_gb*delta_G*delta_B + fG*fE*self.rho_ge*delta_G*delta_E

        return np.sqrt(error_diag + 2*(error_off1 + error_off2 + error_off3))


    def correlation_matrix(self):
        # print correlation matrix        
        f = open('./tables/correlation_'+self.name+'.tex', "w")
        print("\\begin{tabular}{c  c  c  c  c  c}\hline \hline", file=f)
        print("     &  $K$           &     $A$        &      $B$        &      $E$       & $G$ \\\ \hline", file=f)
        print(" $K$ &     1          &                &                 &                &   \\\ ", file=f)
        print(" $A$ &",round(self.rho_ak,2), "&      1         &                 &                &   \\\ ", file=f)
        print(" $B$ &",round(self.rho_bk,2), "&", round(self.rho_ba,2),"&      1          &                &   \\\ ", file=f)
        print(" $E$ &",round(self.rho_ek,2), "&", round(self.rho_ea,2),"&", round(self.rho_eb,2), "&      1         &   \\\ ", file=f)
        print(" $G$ &",round(self.rho_gk,2), "&", round(self.rho_ga,2),"&", round(self.rho_gb,2), "&", round(self.rho_ge,2),"& 1 \\\ \hline \hline", file=f)
        print('\end{tabular}', file=f)
        print("\caption{Correlation parameters. States:", self.name, "}",file=f)
        print("\label{tab:"+self.name+"_corr}", file=f)


    def param_comparison(self):
        
        LA = du.linear_algebra_check(self.name)
        f = open('./tables/parameters_'+self.name+'.tex', "w")
        print("\\begin{tabular}{c | c  c  c  c  c }\hline \hline", file=f)
        print("          & $K$        & $A$             & $B$     & $E$        & $G$            \\\ \hline", file=f)
        print("Paper     & 5727.12$\pm x.xx$ & 21.54$\pm 0.37$ & 23.91$\pm 0.31$ & 30.34 $\pm 0.23$ & 54.37$\pm 0.58$ \\\ ", file=f)
        print("Sampled   &",round(self.Kp,1),' $\\pm',round(self.delta_Kp,2),'$ &', round(self.A,2),' $\\pm',round(self.delta_A,2),'$ &', round(self.B,2),' $\\pm',round(self.delta_B,2),'$ &', round(self.E,2), ' $\\pm',round(self.delta_E,2),'$ &', round(self.G,2), ' $\\pm',round(self.delta_G,2),'$ \\\ ', file=f)
        print("L.Algebra &",round(LA[0][0],1),' &', round(LA[1][0],2),' &', round(LA[2][0],2),' &', round(LA[3][0],2), '&', round(LA[4][0],2), ' \\\ ', file=f)
        print("\hline\hline", file=f)
        print("\end{tabular}", file=f)
        print("\caption{Model pararameters in MeV, for states: $", self.name, "$}",file=f)
        print("\label{tab:"+self.name+"_param}", file=f)
        f.close()
        #print(Kp/5727.12, A/21.54, B/23.91, E/30.34,G/54.37)
        

    def latex_header(self,table_file,paper):        
        if not paper:
            print("\\begin{tabular}{c | c  c  c  c  c}\hline \hline", file=table_file)
            print("Mass State & Experiment  &   Predicted mass  &    Predicted mass & diff pred & diff sampl\\\ ", file=table_file)
            print("           & (MeV)       &   old (MeV)       &    sampled (MeV)  &     (\\%) &     (\\%)  \\\ \hline", file=table_file)
        else:
            print("\\begin{tabular}{c | c  c  }\hline \hline", file=table_file)
            print("Mass State & Experiment  &   Predicted mass  \\\ ", file=table_file)
            print("           & (MeV)       &   sampled (MeV)   \\\ \hline", file=table_file)
                    

    def latex_bottom(self, table_file,diff_pred, diff_sample,paper):
        label = 'paper'
        if not paper:
            print('\hline', file=table_file)
            print("  &  &  & Total diff & ",round(diff_pred), " &" ,round(diff_sample,1),"\\\ ", file=table_file)
            label = 'note'
            
        print('\hline \hline', file=table_file)
        print('\end{tabular}', file=table_file)
        print("\caption{Every quantity is in MeV, except for percentage differences. States:", self.name, "}",file=table_file)
        print("\label{tab:"+self.name+"_mass_"+label+"}", file=table_file)


    def plot(self):
        # plot the simulated sampling distribution,
        # under the Central Limit Theorem, it is expected normal
        dv.plot(self.sampled_k,'k','Parameter omega','charm', self.name)
        dv.plot(self.sampled_a,'a','Parameter A','charm', self.name)
        dv.plot(self.sampled_b,'b','Parameter B','charm', self.name)
        dv.plot(self.sampled_e,'e','Parameter E','charm', self.name)
        dv.plot(self.sampled_g,'g','Parameter G','charm', self.name)

    
    def fetch_values(self):
        # fetch the values from the input dictionaries    
        # sampled
        self.sampled_k = self.sampled['sampled_k']
        self.sampled_a = self.sampled['sampled_a']
        self.sampled_b = self.sampled['sampled_b']
        self.sampled_e = self.sampled['sampled_e']
        self.sampled_g = self.sampled['sampled_g']
                
        self.Kp,self.delta_Kp= np.mean(self.sampled['sampled_k']), np.std(self.sampled['sampled_k'],ddof=1)
        self.A, self.delta_A = np.mean(self.sampled['sampled_a']), np.std(self.sampled['sampled_a'],ddof=1)
        self.B, self.delta_B = np.mean(self.sampled['sampled_b']), np.std(self.sampled['sampled_b'],ddof=1)
        self.E, self.delta_E = np.mean(self.sampled['sampled_e']), np.std(self.sampled['sampled_e'],ddof=1)
        self.G, self.delta_G = np.mean(self.sampled['sampled_g']), np.std(self.sampled['sampled_g'],ddof=1)

        # parameters
        self.sum_mass= self.params['mass']
        self.v_param = self.params['V']
        self.w_param = self.params['W']
        self.x_param = self.params['X']
        self.y_param = self.params['Y']
        self.z_param = self.params['Z']
        
        # correlation matrix
        self.rho_ak = np.mean(self.corr_mat['rho_ak'])
        self.rho_bk = np.mean(self.corr_mat['rho_bk'])
        self.rho_ba = np.mean(self.corr_mat['rho_ba'])
        self.rho_ek = np.mean(self.corr_mat['rho_ek'])
        self.rho_ea = np.mean(self.corr_mat['rho_ea'])
        self.rho_eb = np.mean(self.corr_mat['rho_eb'])
        self.rho_gk = np.mean(self.corr_mat['rho_gk'])
        self.rho_ga = np.mean(self.corr_mat['rho_ga'])
        self.rho_gb = np.mean(self.corr_mat['rho_gb'])
        self.rho_ge = np.mean(self.corr_mat['rho_ge'])


# print("Mass State & Experiment  &   Predicted mass  &    Predicted mass \\\  ")
# print("           & (MeV)       &   old (MeV)       &    sampled (MeV) \\\ \hline  ")
# mass = self.sum_mass[0] + self.sampled_k[i]*self.v_param[0] + self.sampled_a[i]*self.w_param[0] + self.sampled_b[i]*self.x_param[0] + self.sampled_e[i]*self.y_param[0] + self.sampled_g[i]*self.z_param[0]
# error = np.sqrt((self.delta_Kp*self.v_param[i])**2 + (self.delta_A*self.w_param[i])**2 + (self.delta_B*self.x_param[i])**2 + (self.delta_E*self.y_param[i])**2 + (self.delta_G*self.z_param[i])**2 )
