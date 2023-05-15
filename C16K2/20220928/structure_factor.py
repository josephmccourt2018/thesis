import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.style.use('veusz')
from scipy.interpolate import UnivariateSpline as us
from scipy.integrate import trapz


form=pd.read_csv("sample140.csv.fit",index_col=0) 
data=pd.read_csv("sample140.csv.scaled",index_col=0)
structure=pd.DataFrame({"q":data[data.columns[0]],"S":data[data.columns[1]].values/form[form.columns[1]].values})

plt.loglog(form[form.columns[0]],form[form.columns[1]])
plt.loglog(data[data.columns[0]],data[data.columns[1]]) 
plt.loglog(structure[structure.columns[0]],structure[structure.columns[1]])
plt.show()


xlim=[0,.05]

structure=structure[structure[structure.columns[0]]>xlim[0]]
structure=structure[structure[structure.columns[0]]<xlim[1]]


spl = us(structure[structure.columns[0]].values, structure[structure.columns[1]].values,s=0)

plt.loglog(structure[structure.columns[0]],structure[structure.columns[1]])
plt.loglog(structure[structure.columns[0]], spl(structure[structure.columns[0]]))
plt.show()

rs=np.linspace(2*np.pi/xlim[1],1000,100)
gr=[]
for r in rs:
    xs=np.linspace(0,0.05,10000)
    ys=np.array([q*(spl(q)-1)*np.sin(q*r) for q in xs])
    gr.append(1+trapz(xs,ys)/(2*np.pi**(2)*r))
gr=np.array(gr)
print(gr**2)
plt.plot(rs,abs(gr))
plt.show()
