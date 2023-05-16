import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('veusz')
from scipy.optimize import curve_fit

df0=pd.read_excel('20220914_zeta/c16k2.xlsx',sheet_name='sample1-zeta')
df5=pd.read_excel('20220914_zeta/c16k2.xlsx',sheet_name='sample5-zeta')
df10=pd.read_excel('20220914_zeta/c16k2.xlsx',sheet_name='sample6-zeta')
df25=pd.read_excel('20220914_zeta/c16k2.xlsx',sheet_name='sample3-zeta') 
df50=pd.read_excel('20220914_zeta/c16k2.xlsx',sheet_name='sample4-zeta') 
df250=pd.read_excel('20220914_zeta/c16k2.xlsx',sheet_name='sample7-zeta')
df_list=[df0,df5,df10,df25,df50,df250]

def gauss2(x, x_max, Delta, del_eps_max):
    return del_eps_max*np.exp(-(x-x_max)**2/(Delta)**2)

def gauss2_super(x, x_max_1, Delta_1, del_eps_max_1, x_max_2, Delta_2, del_eps_max_2):
    return gauss2(x, x_max_1, Delta_1, del_eps_max_1) + gauss2(x, x_max_2, Delta_2, del_eps_max_2)

count = 0
ionic_strength=[8,13,18,33,58,258]
for df in df_list:
    #guess=[40,5,1e5,60,3,1e3]
    guess=[10,5,1e5,70,5,1e3] # for df250, need different initial guess
    xdata=df[df.columns[0]].values
    ydata=df[df.columns[1]].values
    popt, pcov = curve_fit(gauss2_super, xdata, ydata,p0=guess,bounds=[(-200,-10,-1e5,-200,-10,-1e5),(200,65,1e10,200,65,1e10)])
    plt.scatter(xdata,ydata)#,label=ionic_strength[count])
    plt.plot(xdata,ydata,label=ionic_strength[count])
    count+=1
    plt.plot(xdata,gauss2_super(xdata,*popt))
    print(popt)
    total_scale=popt[2]+popt[5]
    peak_area2=2.5609*(popt[5]*popt[4])
    peak_area1=2.5069*(popt[2]*popt[1])
    print(peak_area1/peak_area2)
    total_scale=peak_area1+peak_area2
    term1=popt[0]*(peak_area1/total_scale)
    term2=popt[3]*(peak_area2/total_scale)
    print("zeta={}".format((term1+term2)))
plt.xlabel('Zeta Potential [mV]')
plt.ylabel('Counts')
plt.legend(title='Ionic Strength')
plt.show()
    
    
