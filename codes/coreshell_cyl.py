#############################################################################
# 20221209
# Joey McCourt
#
# Core shell cylindrical micelle form factor simulation
# https://www.ncnr.nist.gov/resources/sansmodels/CoreShellCylinder.html#:~:text=The%20shell%20thickness,%20t,%20is%20considered%20to%20be,average%20over%20all%20possible%20orientations%20of%20the%20cylinder.
#
#
#############################################################################

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('veusz')
import sys
import os
from scipy.special import jv
from scipy import integrate


q_values=np.arange(0.1,10,0.1)
##########################################
#sample977.csv.fit
#t=.90
#R=1.54
#H=200
#L=2*H
#rho_core=270
#rho_shell=450
#rho_solv=334
#scale=0.244
#bkg=0.0024
##########################################
##########################################
#sample176.csv.fit
#t=1.0
#R=1.64
#H=200
#L=2*H
#rho_core=270
#rho_shell=450
#rho_solv=334
#scale=0.224
#bkg=0.0014
##########################################
##########################################
#sample162.csv.fit
t=0.95
R=1.58
H=200
L=2*H
rho_core=270
rho_shell=450
rho_solv=334
scale=0.484
bkg=0.0026
##########################################

def j0(x):
    return np.sin(x)/x

def cyl(q,R,v,H,alpha):
    return v*j0(q*H*np.cos(alpha))*jv(1,q*R*np.sin(alpha))/(q*R*np.sin(alpha))

def vol_shell(R,t,Ltotal):
    return np.pi*Ltotal*(R+t)**2

def vol_core(R,L):
    return np.pi*L*R**2

def f1(q,alpha,t,R,H,L):
    return 2*(rho_core-rho_shell)*cyl(q,R,vol_core(R,L),H,alpha)+2*(rho_shell-rho_solv)*cyl(q,R+t,vol_shell(R,t,2*H+2*t),H+t,alpha)


def p(q,t,R,H,L):
    integrand = lambda x: np.sin(x)*f1(q,x,t,R,H,L)**2

    #for approximate integral instead of quad
    # Define bounds of integral
    a = 0.001
    b = np.pi/2
    # Generate function values
    x_range = np.arange(a,b+0.001,.001)
    fx = integrand(x_range)
    # Approximate integral
    approx = integral_approximation(fx,a,b)

    #print('\n')
    #print((1/vol_shell(R,t,2*H+2*t))*integrate.quad(integrand,0,np.pi/2)[0])
    #print((1/vol_shell(R,t,2*H+2*t))*approx)
    
    return (1/vol_shell(R,t,2*H+2*t))*approx


# Our integral approximation function
def integral_approximation(f, a, b):
    return (b-a)*np.mean(f)



#bkg=1e-3
intensity=[]
for q in q_values:
    intensity.append(p(q,t,R,H,L))

plt.figure()
plt.loglog(q_values,[scale*i/np.max(intensity)+bkg for i in intensity])
df_sim=pd.DataFrame({"q":q_values, "Intensity":[scale*i/np.max(intensity)+bkg for i in intensity]})
#plt.show()

#basedir='/media/joey/DATA/Office_OptiPlex_7050_Documents/bedzyk_research_group/reduced_data/C16K2/20220928/'
#f='sample148.csv'
#f='sample164.csv'
f='/home/joey/Documents/bedzyk_research_group/reduced_data/C16K2/20220928/sample162.csv'
#f='/home/joey/Documents/bedzyk_research_group/reduced_data/C16K2/20230208/sample977.csv'
#f='/media/joey/DATA/Office_OptiPlex_7050_Documents/bedzyk_research_group/reduced_data/C16K2/20220928/sample148.csv'
filename=f

df_sim.to_csv(filename+".fit")

df = pd.read_csv(filename,names=['q','I','error'],comment='#',header=1)
xlim=[0.005,1]
df=df[df[df.columns[0]]>xlim[0]]
df=df[df[df.columns[0]]<xlim[1]]
xdata=df[df.columns[0]].values*10
ydata=df[df.columns[1]].values/np.max(df[df.columns[1]])
plt.loglog(xdata,ydata,alpha=0.5)
plt.show()
