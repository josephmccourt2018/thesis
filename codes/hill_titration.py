import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.set_option('mode.chained_assignment', None) # ***https://stackoverflow.com/questions/44028898/a-value-is-trying-to-be-set-on-a-copy-of-a-slice-from-a-dataframe-pandas***
plt.style.use('veusz')
from scipy.optimize import curve_fit
import sys
from symfit import Parameter, Variable
from symfit import Fit
from scipy.interpolate import UnivariateSpline as us


# Input parameters, initial volume of sample (acid, Va) and initial conc. of
# titrant (Cb)
Cb = 0.1
# Volume of PA (Va) and salt line up with files
# use command python c16k2/[5mM].dat c16k2/[50mM].dat c16k2/[500mM].dat ...
#
# python titration_fit.py c16k2/Chengrui.dat c16k2/L_C16KK_3.dat c16k2/L_C16KK_2.dat c16k2/L_C16KK_5.dat c16k2/L_C16KK_4.dat c16k2/L_C16KK.dat
#
# The volumes are a bit off due to the measurement of sample and actual concentration
Vas=[5000,5800,5000,6556,3500]
Vas=[4800]
Vas=[8430]
# These are artificial volumes to put the equivalence points of the different samples in the same place
salt=[5,50,100,250,500]
salt=[0]
salt=[0]
#where to stop fitting
vol_lim=1.5
pH_start=4.3
pH_finish=11

def main():
    #reading in files, fitting/plotting
    dfs=[read_file(sys.argv[i],"{} M NaOH".format(Cb)) for i in range(1,len(sys.argv))]
    fit_list=[]
    count = 0
    plt.figure()
    for df in dfs:
        #df = df[df[df.columns[1]]<pH_finish]
        #df[df.columns[0]] = df[df.columns[0]].map(lambda x: x/Vas[count])
        equiv_index=equivalence_point(df)
        equiv_vol=df[df.columns[0]].values[equiv_index]
        print(equiv_vol)
        df[df.columns[0]]= df[df.columns[0]]/equiv_vol
        df = df[df[df.columns[0]]<vol_lim]
        plot_df(df)
        #plt.scatter(df[df.columns[0]],df[df.columns[1]],label=salt[count])
        guess=[[9.2,1.0,0.04],[7.5,0.6,0.004]]
        guess=[[7.0,0.6,0.004]]
        fit_sym = fit_hill_symfit(df,guess[count],Va=Vas[count],equiv=equiv_vol)
        sim_mono_hill_equiv(pH_start,pH_finish,*fit_sym,Va=Vas[count],equiv=equiv_vol)
        count+=1
        fit_list.append(fit_sym)
    #plt.legend(title='NaCl [mM]')
    plt.xlim(-0.005,np.max(vol_lim))
    plt.xlabel('Vol. 0.1 M NaOH / Vol. Equivalence Point')
    plt.ylabel('pH')
    plt.show()
    plot_fit_params(salt,fit_list,index=0)
    plot_fit_params(salt,fit_list,index=1)

    ##ionic strengths
    #count = 0
    #for df in dfs:
    #    ionic_strengths(df,Vas[count],salt[count])
    #    count+=1

def plot_fit_params(salt,fit_list,index=0):
    plt.figure()
    plt.xlabel('NaCl [mM]')
    if index==0:
        plt.ylabel('pKa')
        for i in range(0,len(fit_list)):
            plt.scatter(salt[i],-np.log10(fit_list[i][index]),color='r',s=100)
    else:
        plt.ylabel('n')
        for i in range(0,len(fit_list)):
            plt.scatter(salt[i],fit_list[i][index],color='r',s=100)
    plt.show()
    
    
def equivalence_point(df,skip=18):
    # calculates the equivalence point by looking for maxium of the derivative (inflection point)
    # skips the first skip=18 points due to initial sharp increase for titration curves
    xdata=df[df.columns[0]].values
    ydata=df[df.columns[1]].values
    delx=[]
    dely=[]
    for i in range(1,len(xdata)-1):
        delx.append(xdata[i+1]-xdata[i])
        dely.append(ydata[i+1]-ydata[i])
    deriv = [dely[i]/delx[i] for i in range(0,len(delx))]
    delderiv=[]
    for i in range(0,len(delx)-1):
        delderiv.append(deriv[i+1]-deriv[i])
            
    delderiv = [delderiv[i]/delx[i] for i in range(0,len(delx)-1)] 
    f = us(xdata[skip:-1],deriv[skip-2:-1],s=0)
    f.set_smoothing_factor(0.00005)
    equivalence_index = skip+np.argmax(f(xdata[skip:-1]))
    return equivalence_index
    
    
def read_file(filename, titrant):
    # Col1: total added titrant (uL)
    # Col2: pH
    # Col3: time btwn addition of titrant and reading of pH meter (min)
    return pd.read_csv(filename,sep = '\s+',comment='#',\
                       names = ["Total Added {} (uL)".format(titrant),"pH",\
                                "time btwn (min)"])

def plot_df(*args):
    for df in args:
        plt.scatter(df[df.columns[0]],df[df.columns[1]])
        plt.xlabel(df.columns[0],fontsize=18)
        plt.ylabel(df.columns[1],fontsize=18)
        plot_default()

def plot_default():
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(fontsize=18)
    
# monoprotic hill model for titration: Vb([H]), volume of base as a function of
# H+ ion conc.
# params are Ka (pKa), n (hill), Ainit (initial conc. of acid), Va (initial vol
# of acid), Cb (conc. of titrant (base))
def mono_hill(H,Ka,n,Ainit,Va=Vas[0]):
    Kw = 10**(-14) #pK water
    Binit=Cb
    return -(((H**(2 + n) - Ainit*H*Ka**n + H**2*Ka**n - H**n*Kw - Ka**n*Kw)*Va)\
             /((H**n + Ka**n)*(Binit*H + H**2 - Kw)))

def mono_hill_red(H,Ka,n,Ainit,Va=Vas[0]):
    Kw = 10**(-14) #pK water
    Binit=Cb
    return (-(((H**(2 + n) - Ainit*H*Ka**n + H**2*Ka**n - H**n*Kw - Ka**n*Kw)*Va)\
             /((H**n + Ka**n)*(Binit*H + H**2 - Kw))))/Va

def mono_hill_equiv(H,Ka,n,Ainit,Va=Vas[0],equiv=1000):
    Kw = 10**(-14) #pK water
    Binit=Cb
    print(equiv)
    return -(((H**(2 + n) - Ainit*H*Ka**n + H**2*Ka**n - H**n*Kw - Ka**n*Kw)*Va)\
             /((H**n + Ka**n)*(Binit*H + H**2 - Kw)))/equiv


def sim_mono_hill(start_pH, end_pH, *args):
    pH_range = np.linspace(start_pH,end_pH,40)
    H = [10**(-pH) for pH in pH_range]
    params=[args[i] for i in range(0, len(args))]
    plt.xlabel('Added {} M  NaOH ($\mu$L)'.format(Cb),fontsize=18)
    plt.ylabel('pH',fontsize=18)
    pKa, n, Ainit=np.round(-np.log10(params[0]),2), np.round(params[1],2),\
                            np.round(params[2],4)
    plt.plot(mono_hill(10**(-pH_range),*params),pH_range, \
             label="pKa = {}, n = {}, conc (A) = {}".format(pKa,n, Ainit))
    plot_default()

def sim_mono_hill_red(start_pH, end_pH, *args):
    pH_range = np.linspace(start_pH,end_pH,40)
    H = [10**(-pH) for pH in pH_range]
    params=[args[i] for i in range(0, len(args))]
    plt.xlabel('Added {} M  NaOH ($\mu$L)'.format(Cb),fontsize=18)
    plt.ylabel('pH',fontsize=18)
    pKa, n, Ainit=np.round(-np.log10(params[0]),2), np.round(params[1],2),\
                            np.round(params[2],4)
    plt.plot(mono_hill_red(10**(-pH_range),*params),pH_range, \
             label="pKa = {}, n = {}, conc (A) = {}".format(pKa,n, Ainit))
    plot_default()


def sim_mono_hill_equiv(start_pH, end_pH, *args,**kwargs):
    pH_range = np.linspace(start_pH,end_pH,40)
    H = [10**(-pH) for pH in pH_range]
    params=[args[i] for i in range(0, len(args))]
    plt.xlabel('Added {} M  NaOH ($\mu$L)'.format(Cb),fontsize=18)
    plt.ylabel('pH',fontsize=18)
    pKa, n, Ainit=np.round(-np.log10(params[0]),2), np.round(params[1],2),\
                            np.round(params[2],4)
    plt.plot(mono_hill_equiv(10**(-pH_range),*params,**kwargs),pH_range, \
             label="pKa = {}, n = {}, conc (A) = {}".format(pKa,n, Ainit))
    plot_default()
    #df=pd.DataFrame({"vol_naoh":mono_hill_equiv(10**(-pH_range),*params),"pH":pH_range})
    #df.to_csv('sim_pKa_8.2_HH.csv')
    

def fit_hill(df,guess,**kwargs):
    fit,cov= curve_fit(mono_hill,10**(-df[df.columns[1]][1:]),df[df.columns[0]][1:],\
                       p0=guess,bounds=[[0,0,0.004],[1,4,0.01]],method='trf',\
                       sigma=(df[df.columns[0]][1:]))
    return fit

def fit_hill_red(df,guess,**kwargs):
    fit,cov= curve_fit(mono_hill_red,10**(-df[df.columns[1]][1:]),df[df.columns[0]][1:],\
                       p0=guess,bounds=[[10**(-10),0.5,0.005],[10**(-7),2.5,0.008]],method='trf',\
                       sigma=(df[df.columns[0]][1:]))
    return fit

def fit_hill_symfit(df,guess,**kwargs):
    Ka = Parameter('Ka',value=10**(-guess[0]),min=10**(-10),max=10**(-6))
    n = Parameter('n',value=guess[1],min=0.4,max=2.0)
    #n = Parameter('n',value=guess[1],min=0.95,max=1.05)
    Ainit = Parameter('Ainit',value=guess[2],min=0.004,max=0.1)
    H = Variable('H')
    model = mono_hill_equiv(H,Ka,n,Ainit,**kwargs)
    xdata=10**(-df[df.columns[1]][1:]).values
    ydata=df[df.columns[0]][1:].values
    fit = Fit(model,xdata,ydata)
    fit_result = fit.execute()
    fit_array = [fit_result.value(Ka),fit_result.value(n),fit_result.value(Ainit)]

    return fit_array


def fit_hill_symfit(df,guess,*args,**kwargs):
    Ka = Parameter('Ka',value=10**(-guess[0]),min=10**(-10),max=10**(-6))
    n = Parameter('n',value=guess[1],min=0.2,max=1.1)
    #n = Parameter('n',value=guess[1],min=0.95,max=1.05)
    Ainit = Parameter('Ainit',value=guess[2],min=0.001,max=0.01)
    H = Variable('H')
    model = mono_hill_equiv(H,Ka,n,Ainit,*args,**kwargs)
    xdata=10**(-df[df.columns[1]][1:]).values
    ydata=df[df.columns[0]][1:].values
    fit = Fit(model,xdata,ydata)
    fit_result = fit.execute()
    fit_array = [fit_result.value(Ka),fit_result.value(n),fit_result.value(Ainit)]

    return fit_array



def plot_vol_frac(label_list,*args):
    count=0 # counting number of arguments to match up to volumes
    for df in args:
         df = df[df[df.columns[1]]<pH_finish]
         Vb_frac=[v/Vas[count] for v in df[df.columns[0]].values]
         plt.scatter(Vb_frac, df[df.columns[1]],label="{}".format(label_list[count]))
         count+=1
    plt.xlabel('Vol Frac NaOH')
    plt.ylabel('pH')
    plt.legend(title='NaCl [mM]',fontsize=18)
    plt.show()

def ionic_strengths(df,Va,salt):
    # Calculating ionic strength
    Kw = 10**(-14)
    H = [10**(-pH) for pH in df[df.columns[1]]]
    OH = [Kw/h for h in H]
    X=[0.004 for h in H]
    #NaCl
    NaCl=salt/1000 # Mols/L
    Cl=NaCl
    Na=[Cb*x/1e6/Va + NaCl for x in df[df.columns[0]]]
    
    plt.figure()
    xdata=df[df.columns[1]]
    plt.plot(xdata,H,label="H")
    plt.plot(xdata,OH,label="OH")
    plt.plot(xdata,X,label="X")
    plt.plot(xdata,Na,label="Na")
    plt.xlabel('pH')
    plt.ylabel('conc [M]')
    plt.legend()
    plt.show()
    
    i1=[H[i]+OH[i]+X[i]+Na[i] +Cl for i in range(0,len(H))]
    Ds1 = [0.304/np.sqrt(i) for i in i1]
    i2=[H[i]+OH[i]+Na[i] + Cl for i in range(0,len(H))]
    Ds2 = [0.304/np.sqrt(i) for i in i2]
    plt.figure()
    plt.plot(df[df.columns[1]],i1,label="including X-")
    plt.plot(df[df.columns[1]],i2,label="not including X-")
    plt.xlabel('pH')
    plt.ylabel('ionic strength [M]')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(df[df.columns[1]],Ds1,label="including X-")
    plt.plot(df[df.columns[1]],Ds2,label="not including X-")
    plt.xlabel('pH')
    plt.ylabel('screening length (nm)')
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.plot(df[df.columns[1]],Ds1,label="including X-")
    plt.xlabel('pH')
    plt.ylabel('screening length (nm)')
    plt.legend()
    plt.show()

main()
