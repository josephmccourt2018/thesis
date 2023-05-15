#############################################################################
# 20220506
# Joey McCourt
#
# Bilayer form factor simulation
#
#
#############################################################################

import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib.pyplot as plt
plt.style.use('veusz')
import sys
import os

rho_shell = 410 #e/nm^3
rho_core = 318
rho_solv = 334


q_values=np.arange(0.01,10,0.01)

def f(q,t,h):
    result= ((rho_shell-rho_solv)*2*np.sin(q*t/2)+(rho_core-rho_shell)*2*np.sin(q*(t-2*h)/2))/q
    return result


def form_factor_squared(q,t,h):
    return (f(q,t,h)/q)**2


###################################
#sample152.csv.fit
#t=3.45
#h=1.35
#t = float(sys.argv[1])
#h = float(sys.argv[2])
#scale_factor=0.42e-7
#background=2.3e-4
###################################

###################################
#sample181.csv.fit
#t=3.92
#h=1.68
#t = float(sys.argv[1])
#h = float(sys.argv[2])
#scale_factor=0.2e-7
#background=2.3e-4
###################################

###################################
#sample164.csv.fit
#t=4.2
#h=1.5
#t = float(sys.argv[1])
#h = float(sys.argv[2])
#scale_factor=.5e-7
#background=3.7e-4
###################################


####################################
#sample102.csv.fit
#[4.56573477e+00 1.46979083e+00 1.03226178e+00 2.40000000e-04]
#t= 4.56573477e+00
#h=1.46979083e+00
#scale_factor=1.03226178e+00
#background=2.40000000e-04


t= 4.565
h=1.469
scale_factor=1.03226178e+00
background=2.80000000e-04

intensity=[]
for q in q_values:
    intensity.append(form_factor_squared(q,t,h)*scale_factor+background)

plt.figure()
plt.loglog(q_values,[intensity for intensity in intensity],label='Total Thickness = {}, Headgroup Thickness = {}'.format(t, h))
#filename='/home/joey/Documents/bedzyk_research_group/reduced_data/C16K2/20220928/sample164.csv'



#df_sim=pd.DataFrame({"q":q_values, "Intensity":[scale_factor*i/np.max(intensity)+background for i in intensity]})
#df_sim=pd.DataFrame({"q":q_values, "Intensity":[i for i in intensity]})
#df_sim.to_csv(filename+".fit")




#df = pd.read_csv(filename,names=['q','I','error'],comment='#',header=1)
#plt.loglog(df[df.columns[0]]*10,df[df.columns[1]],label=filename,alpha=0.5)
plt.xlim(0.01,10)
plt.legend(fontsize=10)
plt.show()
