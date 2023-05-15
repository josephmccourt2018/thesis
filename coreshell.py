#############################################################################
# 20220506
# Joey McCourt
#
# Core shell (micelle) form factor simulation
#
#
#############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('veusz')
import sys
import os

rho_shell = 450 #e/nm^3
rho_core = 270
rho_solv = 334


q_values=np.arange(0.1,10,0.01)


#https://www.ncnr.nist.gov/resources/sansmodels/CoreShell.html
def sphere(q,R):
    return 3*(np.sin(q*R)-q*R*np.cos(q*R))/(q*R)**3

def vol(R):
    return 4*np.pi*R**3/3

def core(q,R):
    return sphere(q,R)*(rho_core-rho_shell)*vol(R)

def shell(q,R,R_shell):
    r = R+R_shell
    return sphere(q,r)*(rho_shell-rho_solv)*vol(r)

def form_factor_squared(q,R,R_shell):
    return (core(q,R)+shell(q,R,R_shell))**2/(vol(R+R_shell))

##############################
#sample146.csv.fit
#r_core=2.19
#r_shell=0.71
r_core = float(sys.argv[1])
r_shell = float(sys.argv[2])
scale_factor=5.2e-8
background=3e-4
#############################

##############################
#sample172.csv.fit
#r_core=2.14
#r_shell=0.59
#r_core = float(sys.argv[1])
#r_shell = float(sys.argv[2])
#scale_factor=1.2e-7
#background=1.7e-4
#############################


##############################
#sample172.csv.fit
#r_core=2.22
#r_shell=0.57
#r_core = float(sys.argv[1])
#r_shell = float(sys.argv[2])
#scale_factor=1.2e-7
#background=1.7e-4
#############################

##############################
#sample142.csv.fit
#r_core=1.88
#r_shell=0.78
#r_core = float(sys.argv[1])
#r_shell = float(sys.argv[2])
#scale_factor=1.25e-8
#background=3e-4
#############################


##############################
#sample140.csv.fit
#r_core=1.91
#r_shell=0.74
#r_core = float(sys.argv[1])
#r_shell = float(sys.argv[2])
#scale_factor=1.02e-8
#background=3e-4
#############################

##############################
#sample144.csv.fit
#r_core=2.09
#r_shell=0.75
#r_core = float(sys.argv[1])
#r_shell = float(sys.argv[2])
#scale_factor=2.52e-8
#background=3e-4
#############################


##############################
#sample160.csv.fit
#r_core=2.03
#r_shell=0.66
r_core = float(sys.argv[1])
r_shell = float(sys.argv[2])
scale_factor=9.72e-8
background=3.7e-4
#############################


intensity=[]
for q in q_values:
    intensity.append(form_factor_squared(q,r_core,r_shell)*scale_factor+background)

plt.figure()
plt.loglog(q_values,[intensity for intensity in intensity],label='Core Radius = {}, Shell Thickness = {}'.format(r_core,r_shell))
filename='/home/joey/Documents/bedzyk_research_group/reduced_data/C16K2/20220928/sample160.csv'



#df_sim=pd.DataFrame({"q":q_values, "Intensity":[scale_factor*i/np.max(intensity)+background for i in intensity]})
df_sim=pd.DataFrame({"q":q_values, "Intensity":[i for i in intensity]})
#df_sim.to_csv(filename+".fit")



df = pd.read_csv(filename,names=['q','I','error'],comment='#',header=1)
plt.loglog(df[df.columns[0]]*10,df[df.columns[1]],label=filename,alpha=0.5)
plt.xlim(0.1,5)
plt.legend(fontsize=10)
plt.show()
