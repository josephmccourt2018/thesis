# Command to use this script:
#
# python plotter_scaled.py sample{1}.csv sample{2}.csv ...
#
#
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('veusz')
import pandas as pd
import sys
import os

plt.figure()
basedir=os.getcwd()


scale=1
max_list=[]
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    max_list.append(np.max(df[df.columns[1]]))

maximum=np.max(max_list)
count=0
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    plt.loglog(df[df.columns[0]]*10,df[df.columns[1]]*scale/maximum,label=sys.argv[i])
    scale*=4
    count+=1

#plt.legend(fontsize=16)
plt.xlim(1,20)
#plt.xticks([np.round(x,0) for x in np.linspace(12,17,6)])#,labels=[str(x) for x in np.linspace(12,17,5)])
plt.ylim(10**(-4),0.02)
plt.xlabel('q [nm$^{-1}$]')
plt.ylabel('Intensity [arb. units]')
plt.show()
