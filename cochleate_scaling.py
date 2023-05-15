import numpy as np
import matplotlib.pyplot as plt
plt.style.use('veusz')
from scipy.optimize import fsolve 
from scipy import optimize

l0=0.348*10**(-11) #chiral coupling strength (l/l0 = 1, maximum)
ls=[l0/x for x in np.linspace(1,6,20)] #varying chiral coupling
#extra points closer to l=0
ls.append(l0/10)
ls.append(l0/15)
ls.append(l0/20)
ls.append(l0/25)
ls.append(l0/50) 
radii=[]
ds=[]
for l in ls:
    #print(l/ls[0]) 
    kappa=160*10**(-21) # bending rigidity = 39kT
    d_const=20.9*10**(-9) # d spacing from SAXS power law (6.40 + 2.05*(c)^(-0.5), c = Molar salt concentration=0.02 M)
    R_const=173*10**(-9) # supposed radius of cochleate

    # solve self consistently for constant d spacing and free energy model
    def func(variables): 
        (d,R)=variables 
        eqn1=d*(R/d)*np.arctan(R/d)-(kappa/l)*(8*(1+(9*d**2/(8*R**2)))/(3*(1+(d**2/(R**2)))**(3/2))+np.arcsinh(R/d)) #Helfrich prost minimization from PNAS 2019 Gao et al.
        eqn2=d-d_const #constant d spacing
        #eqn1=d-((kappa/l)*(8*(1+(9*d**2/(8*R**2)))/(3*(1+(d**2/(R**2)))**(3/2))+np.arcsinh(R/d)))/((R/d)*np.arctan(R/d)) 
        #eqn2=R-R_const
        return [eqn1,eqn2] 
    d0,radius=fsolve(func,(d_const,R_const*100))/10**(-9)
    ds.append(d0)
    radii.append(radius)
    print(l/l0,radius)
    #print(radius)
#plt.plot([l/ls[0] for l in ls],radii, label='R (nm)')
plt.plot([l/ls[0] for l in ls],ds,label='D (nm)',color='k')
plt.xlabel("$\lambda/\lambda_0$")
#plt.legend()
plt.ylabel("D (nm)")
plt.show()                                               
 




# plot gradient line to match with cochleate to flat bilayer coloring
from matplotlib.collections import LineCollection

x=[l/ls[0] for l in ls]
y=np.array(radii)
cols = np.linspace(1,0,len(x))

points = np.array([x, y]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

fig, ax = plt.subplots()
lc = LineCollection(segments, cmap='PuBu',alpha=1)
lc.set_array(cols)
lc.set_linewidth(8)
line=ax.add_collection(lc)
plt.ylim(0,y[-4])
plt.xlabel("$\lambda/\lambda_0$")
plt.ylabel("R (nm)")
#fig.colorbar(line,ax=ax)
plt.show()
