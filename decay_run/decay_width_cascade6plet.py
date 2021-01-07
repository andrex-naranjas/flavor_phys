'''
---------------------------------------------------------------
 Code to calcualte charm-baryon decay widths
 Authors: A. Ramirez-Morales and H. Garcia-Tecocoatzi
 ---------------------------------------------------------------
'''
from charm_wrapper import decay
import decay_utils as du
import numpy as np

# decay widths for A ->B+C
width  = decay()
baryons = []

# baryon FLAG: 1 -> omega, 2->cascade_6, 3->sigma,# 4 -> lambda, 5-> cascade_3 
# ModEx  FLAG: 1 -> lambda, 2->rho  Excitation
# decPr  FLAG: 3->Xi'+Pi, 5->Sigma+K  decayProduct Flag

# cascade sixplet states (->lambda+K)
state1 = {'MA':2.905,'MB':2.286,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':1,'state':1}
state2 = {'MA':2.934,'MB':2.286,'MC':0.493,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':1,'state':2}
state3 = {'MA':2.941,'MB':2.286,'MC':0.493,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':1,'state':3}
state4 = {'MA':2.969,'MB':2.286,'MC':0.493,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':1,'state':4}
state5 = {'MA':3.029,'MB':2.286,'MC':0.493,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':1,'state':5}
state6 = {'MA':3.060,'MB':2.286,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':1,'state':6}
state7 = {'MA':3.096,'MB':2.286,'MC':0.493,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':1,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7)# add states to dict array

# cascade sixplet states (->Xi+Pi)
state1 = {'MA':2.905,'MB':2.469,'MC':0.140,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':2,'state':1}
state2 = {'MA':2.934,'MB':2.469,'MC':0.140,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':2,'state':2}
state3 = {'MA':2.941,'MB':2.469,'MC':0.140,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':2,'state':3}
state4 = {'MA':2.969,'MB':2.469,'MC':0.140,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':2,'state':4}
state5 = {'MA':3.029,'MB':2.469,'MC':0.140,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':2,'state':5}
state6 = {'MA':3.060,'MB':2.469,'MC':0.140,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':2,'state':6}
state7 = {'MA':3.096,'MB':2.469,'MC':0.140,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':2,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7)# add states to dict array

# cascade sixplet states (->Xi'+Pi)
state1 = {'MA':2.905,'MB':2.578,'MC':0.140,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':3,'state':1}
state2 = {'MA':2.934,'MB':2.578,'MC':0.140,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':3,'state':2}
state3 = {'MA':2.941,'MB':2.578,'MC':0.140,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':3,'state':3}
state4 = {'MA':2.969,'MB':2.578,'MC':0.140,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':3,'state':4}
state5 = {'MA':3.029,'MB':2.578,'MC':0.140,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':3,'state':5}
state6 = {'MA':3.060,'MB':2.578,'MC':0.140,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':3,'state':6}
state7 = {'MA':3.096,'MB':2.578,'MC':0.140,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':3,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7)# add states to dict array

# cascade sixplet states (->Xi*+Pi)
state1 = {'MA':2.905,'MB':2.645,'MC':0.140,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':4,'state':1}
state2 = {'MA':2.934,'MB':2.645,'MC':0.140,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':4,'state':2}
state3 = {'MA':2.941,'MB':2.645,'MC':0.140,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':4,'state':3}
state4 = {'MA':2.969,'MB':2.645,'MC':0.140,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':4,'state':4}
state5 = {'MA':3.029,'MB':2.645,'MC':0.140,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':4,'state':5}
state6 = {'MA':3.060,'MB':2.645,'MC':0.140,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':4,'state':6}
state7 = {'MA':3.096,'MB':2.645,'MC':0.140,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':4,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7) # add states to array

# cascade sixplet states (->Sigma+K)
state1 = {'MA':2.905,'MB':2.453,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':5,'state':1}
state2 = {'MA':2.934,'MB':2.453,'MC':0.493,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':5,'state':2}
state3 = {'MA':2.941,'MB':2.453,'MC':0.493,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':5,'state':3}
state4 = {'MA':2.969,'MB':2.453,'MC':0.493,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':5,'state':4}
state5 = {'MA':3.029,'MB':2.453,'MC':0.493,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':5,'state':5}
state6 = {'MA':3.060,'MB':2.453,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':5,'state':6}
state7 = {'MA':3.096,'MB':2.453,'MC':0.493,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':5,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7) # add states to array

# cascade sixplet states (->Sigma+K)
state1 = {'MA':2.905,'MB':2.518,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':6,'state':1}
state2 = {'MA':2.934,'MB':2.518,'MC':0.493,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':6,'state':2}
state3 = {'MA':2.941,'MB':2.518,'MC':0.493,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':6,'state':3}
state4 = {'MA':2.969,'MB':2.518,'MC':0.493,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':6,'state':4}
state5 = {'MA':3.029,'MB':2.518,'MC':0.493,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':6,'state':5}
state6 = {'MA':3.060,'MB':2.518,'MC':0.493,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':6,'state':6}
state7 = {'MA':3.096,'MB':2.518,'MC':0.493,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':6,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7) # add states to array

# cascade sixplet states (->Xi+eta)
state1 = {'MA':2.905,'MB':2.469,'MC':0.548,'SA':1/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':7,'state':1}
state2 = {'MA':2.934,'MB':2.469,'MC':0.548,'SA':3/2,'LA':1,'JA':1/2,'SL':1,'baryon':2,'ModEx':1,'decPr':7,'state':2}
state3 = {'MA':2.941,'MB':2.469,'MC':0.548,'SA':1/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':7,'state':3}
state4 = {'MA':2.969,'MB':2.469,'MC':0.548,'SA':3/2,'LA':1,'JA':3/2,'SL':1,'baryon':2,'ModEx':1,'decPr':7,'state':4}
state5 = {'MA':3.029,'MB':2.469,'MC':0.548,'SA':3/2,'LA':1,'JA':5/2,'SL':1,'baryon':2,'ModEx':1,'decPr':7,'state':5}
state6 = {'MA':3.060,'MB':2.469,'MC':0.548,'SA':1/2,'LA':1,'JA':1/2,'SL':0,'baryon':2,'ModEx':2,'decPr':7,'state':6}
state7 = {'MA':3.096,'MB':2.469,'MC':0.548,'SA':1/2,'LA':1,'JA':3/2,'SL':0,'baryon':2,'ModEx':2,'decPr':7,'state':7}
baryons = du.append_dic(baryons,state1,state2,state3,state4,state5,state6,state7) # add states to array

name_decays=[]
name_decays.append('State')
state1_decays,state2_decays,state3_decays,state4_decays,\
    state5_decays,state6_decays,state7_decays = ([]), ([]),([]),([]),([]),([]),([])
print('--------------------------------------------------------------------------------------------------')
print('Baryon |  Exc. |  Mode   |  MassA |  MassB | MassC |     JA |     LA |     SA |     SL | DecayWidth')
print('--------------------------------------------------------------------------------------------------')
for i in range(len(baryons)):
    MassA = baryons[i]['MA']
    MassB = baryons[i]['MB']
    MassC = baryons[i]['MC']
    SA_qm = baryons[i]['SA']
    LA_qm = baryons[i]['LA']
    JA_qm = baryons[i]['JA']
    SL_qm = baryons[i]['SL']
    baryon= baryons[i]['baryon']
    ModEx = baryons[i]['ModEx']
    decPr = baryons[i]['decPr']
    state = baryons[i]['state']
    
    decay_value = width.decay_width(MassA, MassB, MassC, SA_qm,
                                    LA_qm, JA_qm, SL_qm, baryon, ModEx, decPr)
    baryon_name, ModEx_name, decPr_name = du.state_labels(baryon,ModEx,decPr)
    print('%6s |  %4s | %7s |  %5.3f |  %5.3f | %5.3f |  %5.1f |  %5.1f |  %5.1f |  %5.1f | %5.6f '
          %(baryon_name, ModEx_name, decPr_name, MassA, MassB, MassC, JA_qm, LA_qm, SA_qm, SL_qm,  decay_value))
    if(i%7==6 and i !=0):
        print('--------------------------------------------------------------------------------------------------')

    if(state==1):
        state1_decays = np.append(state1_decays, decay_value)
        state1_name = '1'
    elif(state==2):
        state2_decays = np.append(state2_decays, decay_value)
        state2_name = '2'
    elif(state==3):
        state3_decays = np.append(state3_decays, decay_value)
        state3_name = '3'
    elif(state==4):
        state4_decays = np.append(state4_decays, decay_value)
        state4_name = '4'
    elif(state==5):
        state5_decays = np.append(state5_decays, decay_value)
        state5_name = '5'
    elif(state==6):
        state6_decays = np.append(state6_decays, decay_value)
        state6_name = '6'
    elif(state==7):
        state7_decays = np.append(state7_decays, decay_value)
        state7_name = '7'
        name_decays.append(du.latex_decay_label(baryon,decPr))


print('TOTAL WIDTH')
print('state1 %5.6f |  state2 %5.6f |  state3 %5.6f |  state4 %5.6f |  state5 %5.6f |  state6 %5.6f |  state7 %5.6f'
      %(np.sum(state1_decays),np.sum(state2_decays),np.sum(state3_decays),np.sum(state4_decays),np.sum(state5_decays),np.sum(state6_decays),np.sum(state7_decays)))


name_decays.append('Tot $\\Gamma$')
baryons='casc6plet'
f_decay = open('./tables/'+baryons+'_decays.tex', "w")

du.print_header_latex(name_decays, f_decay)

du.print_row_latex(state1_name,state1_decays,f_decay)
du.print_row_latex(state2_name,state2_decays,f_decay)
du.print_row_latex(state3_name,state3_decays,f_decay)
du.print_row_latex(state4_name,state4_decays,f_decay)
du.print_row_latex(state5_name,state5_decays,f_decay)
du.print_row_latex(state6_name,state6_decays,f_decay)
du.print_row_latex(state7_name,state7_decays,f_decay)

du.print_bottom_latex(baryons,f_decay)    
