# program to fit parameters of a physics model

from iminuit import Minuit
import numpy as np

# input parameters 
#param_o = np.array([0.00, 0.00, 337.27, 337.27, 337.27, 337.27, 337.27])
param_o = np.array([0.00, 0.00, 1.00, 1.00, 1.00, 1.00, 1.00])
param_x = np.array([0.75, 3.75, 0.75, 3.75, 0.75, 3.75, 3.75])
param_y = np.array([0.00, 0.00, -1.0, -2.5, 0.50, -1.0, 1.50])
param_z = np.array([3.33, 3.33, 3.33, 3.33, 3.33, 3.33, 3.33])
param_m = np.array([2695, 2770, 3000, 3050, 3066, 3090, 3188])

def model(w,x,y,z,o,a,b,g):
    return 2505 + w*o + x*a + y*b + z*g

def least_squares(o, a, b, g):
    yvar = 0.01
    comp_m = model(param_o, param_x, param_y, param_z, o, a, b, g)
    return np.sum((param_m - comp_m)**2 / yvar)

m=Minuit(least_squares, o=0, a=0, b=0, g=0, error_o=1, error_a=1, error_b=1, error_g = 1, errordef=1)
m.migrad()

print(m.get_param_states())

print("W-> Paper value : {}, Fit value: {}, Difference: {}%".format(337.3,  round(m.values['o'],1),  int((1 - m.values['o']/337.3)*100)))
print("A-> Paper value : {}, Fit value: {}, Difference: {}%".format(21.540, round(m.values['a'],2),  int((1 - m.values['a']/21.54)*100)))
print("B-> Paper value : {}, Fit value: {}, Difference: {}%".format(23.910, round(m.values['b'],2),  int((1 - m.values['b']/23.91)*100)))
print("G-> Paper value : {}, Fit value: {}, Difference: {}%".format(54.370, round(m.values['g'],2),  int((1 - m.values['g']/54.37)*100)))
