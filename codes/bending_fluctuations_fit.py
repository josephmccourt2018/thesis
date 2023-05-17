import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('veusz')
import sys
from scipy.optimize import curve_fit
import MDAnalysis as mda



#filename=sys.argv[1]
charge=["0"]
tail_length="5"

kT=4.142*10**(-21)
n=16



#for filename in sys.argv[1:]:
for c in charge:
    file1=c+"/"+tail_length+"/bilayer/4x_bilayer/md_2.tpr"
    file2=c+"/"+tail_length+"/bilayer/4x_bilayer/md_2_skip10.xtc"
    u = mda.Universe(file1,file2)
    # voxelizing the simulation box with rectangular prisms
    Lx=u.dimensions[0]#/10 #nm
    Ly=u.dimensions[1]#/10 #nm
    Lz=u.dimensions[2]#/10 #nm
    def power(q,kappa):
        #n**4 to take into account finite fourier transform of dimension n
        #/4 to account for the fact that I used np.pi/L instead of 2*np.pi/L in
        #the FFT2
        #return (kT/(kappa*(q)**4))/(Lx*Ly*10**(-20))
        return (n**4/4)*(kT/(kappa*(q)**4))/(Lx*Ly*10**(-20))


    
    #df_list=[pd.read_csv(c+'/'+tail_length+'/bilayer/4x_bilayer/bending_fluctuations_'+c+'_n8_f1_to_f500.csv',index_col=0),pd.read_csv(c+'/'+tail_length+'/bilayer/4x_bilayer/bending_fluctuations_'+c+'_n16_f1_to_f500.csv',index_col=0)]
    #for df in df_list:
    #    df=df[1:]#get rid of q=0 divergence
    #    xdata=df[df.columns[0]]
    #    ydata=df[df.columns[1]]-np.mean(df[df.columns[1]].values[20:])
    #    plt.scatter(xdata, ydata)#,label=filename)
    #    print(c)
    #    popt,pcov=curve_fit(power,xdata,ydata)
    #    print("kappa (kT)={}".format(popt[0]*10**(-20)/kT)) #10**(-20) since <|hq|**2> is in terms of inverse Angstroms squared from bending fluctuations code
    #    plt.plot(xdata,power(xdata,*popt))



    df = pd.read_csv(c+'/'+tail_length+'/bilayer/4x_bilayer/bending_fluctuations_'+c+'_n16_f100_to_f300.csv',index_col=0)
    df=df[1:]#get rid of q=0 divergence
    xdata=df[df.columns[0]]
    #ydata=df[df.columns[1]]
    ydata=df[df.columns[1]]-np.mean(df[df.columns[1]].values[20:])
    #print(np.mean(df[df.columns[1]].values[20:]))
    plt.scatter(xdata, ydata)#,label=filename)
    #print(filename)
    print(c)
    popt,pcov=curve_fit(power,xdata,ydata)
    print("kappa (kT)={}".format(popt[0]*10**(-20)/kT)) #10**(-20) since <|hq|**2> is in terms of inverse Angstroms squared from bending fluctuations code
    plt.plot(xdata,power(xdata,*popt))

        

plt.show()
