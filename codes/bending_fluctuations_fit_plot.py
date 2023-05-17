import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

plt.style.use('veusz')

#bilayer
x8s=[40]
x12s=[0,10,20,30,40,50]                                                   
x16s=[0,10,20,30,40,50,60,70,80,90]                                     
x20s=[0,10,20,30,40,50,60,70,80,90,100]
y8s=[0.00663]
y12s=[0.01275,0.01515,0.01684,0.01562,0.01648,0.01458]                                  
y16s=[0.02426,0.02549,0.02573,0.02537,0.02323,0.02679,0.01656,0.01764,0.01430,0.01158]
y20s=[0.03021,0.03353,0.03409,0.03191,0.02502,0.02718,0.02386,0.02608,0.01998,0.01660,0.01360]
y8s_kT=[16.40]
y12s_kT=[31.52,37.46,41.63,38.61,40.75,36.05]
y16s_kT=[59.97,63.02,63.63,62.73,57.42,66.24,40.95,43.62,35.35,28.62] 
y20s_kT=[74.69,82.89,84.27,78.89,61.86,67.20,58.98,64.46,49.41,41.06,33.62]

#4x bilayer
x8s=[40]
x12s=[10,20,30,40,50]                                                   
x16s=[10,20,40,50,60,70,80,90]                                     
x20s=[10,20,30,40,50,60,70,80,90,100]
y8s_kT=[18.40]
y12s_kT=[32.84,48.30,46.32,49.50,49.83]
y16s_kT=[68.40,77.23,83.69,92.73,52.77,49.26,53.24,45.70] 
y20s_kT=[104.11,107.64,125.35,121.95,126.20,92.16,70.47,68.29,64.34,49.64]





plt.figure()
plt.legend(fontsize=20)
plt.xlabel("Charge Percent (%)")
plt.ylabel("$\kappa$ ($k_B T$)")
plt.scatter(x8s,y8s_kT,label="C8",s=100)   
plt.scatter(x12s,y12s_kT,label="C12",s=100)                                 
plt.scatter(x16s,y16s_kT,label="C16",s=100)                                 
plt.scatter(x20s,y20s_kT,label="C20",s=100)
plt.legend(fontsize=20)
plt.show()


plt.figure()
plt.legend(fontsize=20)
plt.xlabel("Charge Percent (%)")
plt.ylabel("$\kappa/\kappa_0$")
#plt.scatter(x8s,[y/y8s_kT[0] for y in y8s_kT],label="C8")
#plt.scatter(x12s,[y/y12s_kT[0] for y in y12s_kT],label="C12")
plt.scatter(x16s,[y/y16s_kT[0] for y in y16s_kT],label="C16",color='#4DAF4A',s=100)
plt.scatter(x20s,[y/y20s_kT[0] for y in y20s_kT],label="C20",color='#F884BC',s=100)


#ALL DATA
xs=x12s+x16s+x20s
ys=[y/y12s_kT[0] for y in y12s_kT]+[y/y16s_kT[0] for y in y16s_kT]+[y/y20s_kT[0] for y in y20s_kT]
data=zip(xs,ys)
data=list(set(data))
def sortFirst(val):
    #print(val[1])
    return val[0]
#print(data)
data.sort(key=sortFirst) #sort data by charge percent
#print(data)

xs=list(zip(*data))[0]
ys=list(zip(*data))[1]


#JUST C16 and C20

xs=x16s+x20s
ys=[y/y16s_kT[0] for y in y16s_kT]+[y/y20s_kT[0] for y in y20s_kT]
data=zip(xs,ys)
data=list(set(data))
def sortFirst(val):
    #print(val[1])
    return val[0]
#print(data)
data.sort(key=sortFirst) #sort data by charge percent
#print(data)

xs=list(zip(*data))[0]
ys=list(zip(*data))[1]


def power(x,A,B,m,n):
    b=1
    return A*(x)**m +B*x*n+ b

def line(x,A,m,b):
    return m*(x)**A+b

#guess=[-0.51245198,  0.55577823,  0.92815471,  0.90693002]
#popt,pcov=curve_fit(power,xs,ys,p0=guess)
crossover=8
#popt,pcov=curve_fit(line,xs[:crossover],ys[:crossover])
#popt2,pcov2=curve_fit(line,xs[crossover:],ys[crossover:])
m,b=np.polyfit(xs[:crossover],ys[:crossover],1)
m2,b2=np.polyfit(xs[crossover:],ys[crossover:],1)
m3,b3=np.polyfit(xs[crossover-2:crossover+2],ys[crossover-2:crossover+2],1)

print(m,b,m2,b2)
#plt.plot(xs,power(xs,*popt),label="Fit")

plt.plot(xs[:crossover],[m*x+b for x in xs[:crossover]],color='k',linestyle='-')
plt.plot(xs[crossover:],[m2*x+b2 for x in xs[crossover:]],color='k',linestyle='-')
plt.plot(xs[crossover-2:crossover+2],[m3*x+b3 for x in xs[crossover-2:crossover+2]],color='k',linestyle='--')

#print(popt)

plt.legend(fontsize=20)

plt.show()
