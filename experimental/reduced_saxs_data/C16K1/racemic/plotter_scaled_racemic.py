# Command to use this script:
#
# python plotter_scaled_racemic.py sample353.csv sample370.csv 
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

salts=[0,100]

scale=1
max_list=[]
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    max_list.append(np.max(df[df.columns[1]]))

maximum=np.max(max_list)
count=0
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    plt.loglog(df[df.columns[0]]*10,df[df.columns[1]]*scale/maximum,label=salts[i-1],linewidth=2.5)
    scale*=10
    count+=1

plt.legend(fontsize=20,title='[NaCl] (mM)')
plt.xlabel('q [nm$^{-1}$]')
plt.ylabel('Intensity [arb. units]')
plt.title('Racemic-C$_{16}$K$_1$')
plt.show()
