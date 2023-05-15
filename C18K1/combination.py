from scipy.interpolate import UnivariateSpline as us
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('veusz')
import pandas as pd

df_bil=pd.read_csv('sample102.csv.fit',comment='#',index_col=0)
df_bil_exp=pd.read_csv('sample102.csv',comment='#',index_col=0)
xdata=df_bil_exp[df_bil_exp.columns[0]]*10
ydata=df_bil_exp[df_bil_exp.columns[1]]
xsim=df_bil[df_bil.columns[0]]
ysim=df_bil[df_bil.columns[1]]/1.75
plt.loglog(xdata,ydata)
plt.loglog(xsim,ysim)
us_bil=us(xsim,ysim,s=0)
q_values=xdata[0:710]
plt.loglog(q_values,us_bil(q_values))
#pd.DataFrame({'q':q_values,'Intensity':us_bil(q_values)}).to_csv('sample102.csv.fit_2')
plt.show()


df_hel=pd.read_csv('sample108.csv.fit',comment='#',index_col=0)
df_hel_exp=pd.read_csv('sample108.csv',comment='#',index_col=0)
xdata=df_hel_exp[df_hel_exp.columns[0]]*10
ydata=df_hel_exp[df_hel_exp.columns[1]]
xsim2=df_hel[df_hel.columns[0]]
ysim2=df_hel[df_hel.columns[1]]*4.6
plt.loglog(xdata,ydata)
plt.loglog(xsim2,ysim2)
us_hel=us(xsim2,ysim2,s=0)
q_values=xdata[0:710]
plt.loglog(q_values,us_hel(q_values))
plt.show()

def combination(q,a):
    return a*us_hel(q)+(1-a)*us_bil(q)



df_mix_exp=pd.read_csv('sample105.csv',comment='#',index_col=0)
xdata=df_mix_exp[df_mix_exp.columns[0]]*10
ydata=df_mix_exp[df_mix_exp.columns[1]]
plt.loglog(xdata,ydata)
q_values=xdata[0:710]
plt.loglog(q_values,combination(q_values,0.35))
pd.DataFrame({'q':q_values,'Intensity':combination(q_values,0.35)}).to_csv('sample105.csv.fit_combination')
plt.show()
