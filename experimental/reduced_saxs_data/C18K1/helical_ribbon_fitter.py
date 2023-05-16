# see /home/joey/bedzyk_research_group/form_factors/terech_helical_ribbon.ipynb for more details


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('veusz')
import scipy.integrate as integrate
import scipy.special as sp
from scipy.optimize import curve_fit, leastsq
from scipy.signal import find_peaks
import sys
from tqdm import tqdm
from scipy.interpolate import UnivariateSpline as us




##l = ["6","7","8p7","11"]
#cg=[pd.read_csv("/home/joey/Documents/bedzyk_research_group/reduced_data/C18K1/sample134.csv",names=["q","I","error"],header=1,comment='#')]#for l in ls]
##cg2=[pd.DataFrame.drop_duplicates(c,['q']) for c in cg]
#cg[0]['I']=cg[0]['I']/np.max(cg[0]['I'].values)
#cg



def print_params(p,p_names):
    for i in range(0,len(p)):
        print("{} = {}".format(p_names[i],p[i]))


# define a sort key
def sort_key(company):
    return company[0]



def f1(q,w,n,eps_n):
    return eps_n*(np.sinc(n*w/2/np.pi))**2


def s_poly_bilayer_a2(qvals,R,delP,psi,A):
    intensity=[]
    print(R,delP,psi,A)
    lh=0.755
    rho_h=410 #for c18k
    #rho_h=410 #for c16k
    rho_t=303 #for c18k pH 6
    #rho_t=304 #for c18k high pH
    #rho_t=300 #for c16k
    rho_s=334
    
    a = 0.951 #for c18K pH 6
    #a=0.950 #for c18k high pH
    #a=0.951 #for c16K
    del_R = 0.1*R # 5% polydispersity
    Rs = np.arange(R-del_R,R+del_R,del_R/10)
    intensity=[]
    for q in qvals:    
        temp=0
        intensity1=[]
        for R in Rs:
            b=1/np.tan(psi)
            for n in [0,1,2,3]:
                if n==0:
                    epsilon=1
                    if q>=n*b/R:
                        qn=n*b/(q*R)
                        g = (rho_h-rho_s)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R, R)[0])+(rho_t-rho_h)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R+lh, R-lh)[0])
                    elif q<n*b/R:
                        qn=1
                        g = (rho_h-rho_s)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R, R)[0])+(rho_t-rho_h)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R+lh, R-lh)[0])
                    temp+=g**2*f1(q,2*np.pi*delP,n,epsilon)   
                else:
                    epsilon=2
                    if q>=n*b/R:
                        qn=n*b/(q*R)
                        g = (rho_h-rho_s)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R, R)[0])+(rho_t-rho_h)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R+lh, R-lh)[0])
                    elif q<n*b/R:
                        qn=1
                        g = (rho_h-rho_s)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R, R)[0])+(rho_t-rho_h)*(integrate.quad(lambda x: x*sp.jv(n,q*x*(1-qn**2)**(0.5)), a*R+lh, R-lh)[0])
                    temp+=g**2*f1(q,2*np.pi*delP,n,epsilon)
            intensity1.append(A*temp*(delP**2)*2*np.pi*R*np.tan(psi)/q)
        intensity.append(np.mean(intensity1,axis=0))
    s=us(qvals,intensity/np.max(intensity),s=0)
    return s(qvals)

#R,delP,psi,A; a=0.95, del_R = 20%R, lh=0.755
#C16K1, pH 6 : guess = [66.56594824997335,0.542053606368151,0.5102673464476227,1.3734777470701654e-10]
#
#C18K1, pH 6.2: guess = [78.0,0.53,0.61,1],10% polydispersity, sample108.csv
# fit: 75.00216490093072 0.4907516165706524 0.5954913797188431 0.9949999999992776
# fit (with appropriate MAXS region for bilayer): 67.66831951560755 0.4495890822346444 0.6256326165138937 0.9901668466408708
#
#
#C18K1, pH 7.1: guess = [78.0,0.53,0.61,1], 10% polydispersity, sample111.csv
# fit: 72.4806761548213 0.4826254639634803 0.6035767231382473 0.09494924092544955
# fit (with appropriate MAXS region for bilayer): 72.48067184881398 0.48262546395123296 0.603576723124116 0.09494924092544955
#
#
#C18K1, pH 7.7: guess = [78.0,0.53,0.61,1], 10% polydispersity, sample116.csv
# fit: 73.96786940450153 0.5144775182247749 0.591846581531895 0.44405839036212796
# fit (with appropriate MAXS region for bilayer):73.96786735201256 0.5167022114583343 0.5890575007670479 0.44405839033118527
#
#
#C18K1, pH 8.7: guess = [78.0,0.53,0.61,1], 10% polydispersity, sample123.csv
# fit: 77.7099684856166 0.5358074975187369 0.6086769519284287 0.2623805688206657
# fit (with appropriate MAXS region for bilayer):77.70991453534504 0.5425084008259806 0.602842598840674 0.2623805688550106
#
#
#C18K1, pH 9.2: guess = [78.0,0.53,0.61,1], 10% polydispersity, sample128.csv
# fit: 76.36809293501852 0.5284899771966044 0.5953462257982433 0.13087669959869866
#
#
#C18K1, pH 10.6: guess = [78.0,0.53,0.61,1], 10% polydispersity, sample131.csv
# fit: 79.82284702772685 0.5404406320006936 0.6553594472317784 0.8866709949167255
#
#
#C18K1, pH 11.5: guess = [78.0,0.53,0.61,1], 10% polydispersity, sample134.csv
# fit: 77.6080994437633 0.5325031965032093 0.6325363873650387 0.24540238652360022
#C18K1, pH 11.5: guess = [78.0,0.53,0.61,1], 5% polydispersity, sample134.csv
# fit: 90.1860416104088 0.5335120704888537 0.6802658651551056 0.9511262788350276 (basically a straight line)
# fit (with appropriate MAXS region for bilayer):77.09686768326718 0.5378321760228344 0.6285713257986 0.08582850946296867
#



##guess = [78.0,0.53,0.61,1]
#guess=[77.09686768326718, 0.5378321760228344, 0.6285713257986, 0.08582850946296867]
#xdata=10*cg[0][cg[0].columns[0]].values[0:100]
#ydata=cg[0][cg[0].columns[1]].values[0:100]
##xdata=10*cg[0][cg[0].columns[0]].values[0:710]
##ydata=cg[0][cg[0].columns[1]].values[0:710]
#w= [[k, v] for k, v in zip(xdata, ydata)]
#w.sort(key=sort_key, reverse=False)
#xdata = list(zip(*w))[0]
#ydata = list(zip(*w))[1]
#plt.loglog(xdata,ydata,label="C$_{18}$K")
#bkg=2e-4
#plt.loglog(xdata,s_poly_bilayer_a2(xdata,*guess)+bkg,label="sim")
##popt,pcov = curve_fit(s_poly_bilayer_a2,xdata,ydata,p0=guess,sigma=ydata,bounds=[(40,0,np.pi/12,0),(100,1,np.pi/3,1)])
##plt.loglog(xdata,s_poly_bilayer_a2(xdata,*popt),label="fit")
#plt.legend()
#plt.xlabel('q (nm$^{-1}$)')
##plt.ylabel('Intensity (cm$^{-1}$)')
#plt.ylabel('Intensity/I$_0$')
##plt.title('pH 6')
#print("sim")
#print_params(guess,["R","delta/P","psi","scaling"])
##print("fit")
##print_params(popt,["R","delta/P","psi","scaling"])
#plt.show()






file_list=['sample108.csv','sample111.csv','sample116.csv','sample123.csv','sample128.csv','sample131.csv','sample134.csv']
params=[[73.66831951560755, 0.4495890822346444, 0.6256326165138937, 0.9901668466408708],\
        [72.48067184881398, 0.48262546395123296, 0.603576723124116, 0.09494924092544955],\
        [73.96786735201256, 0.5167022114583343, 0.5890575007670479, 0.44405839033118527],\
        [75.70991453534504, 0.5425084008259806, 0.602842598840674, 0.2623805688550106],\
        [76.36809293501852, 0.5284899771966044, 0.5953462257982433, 0.13087669959869866],\
        [79.82284702772685, 0.5404406320006936, 0.6553594472317784, 0.8866709949167255],\
        [77.09686768326718, 0.5378321760228344, 0.6285713257986, 0.08582850946296867]]


choice=6


cg=[pd.read_csv("/media/joey/DATA/Office_OptiPlex_7050_Documents/bedzyk_research_group/reduced_data/C18K1/"+file_list[choice],names=["q","I","error"],header=1,comment='#')]
cg[0]['I']=cg[0]['I']/np.max(cg[0]['I'].values)

xdata=10*cg[0][cg[0].columns[0]].values[0:710]
ydata=cg[0][cg[0].columns[1]].values[0:710]
w= [[k, v] for k, v in zip(xdata, ydata)]
w.sort(key=sort_key, reverse=False)
xdata = list(zip(*w))[0]
ydata = list(zip(*w))[1]
plt.loglog(xdata,ydata,label="C$_{18}$K")



bkg=2.7e-4
sim = s_poly_bilayer_a2(xdata,*params[choice])+bkg
plt.loglog(xdata,sim,label="sim")
###############################################
#TO SAVE FILES
#df_sim=pd.DataFrame({"q":xdata, "Intensity":sim})
#df_sim.to_csv(file_list[choice]+".fit")
###############################################

plt.show()
