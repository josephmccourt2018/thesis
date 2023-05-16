# Command to use this script:
#
# python plotter_scaled.py sample646.csv sample647.csv sample648.csv sample650.csv sample651.csv sample652.csv
# 100% L to 50%L/50%D C16K1 (racemic)
#
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('veusz')
import pandas as pd
import sys
import os

plt.figure()
basedir=os.getcwd()


ls=[100,88.8,81.7,71.8,55.4,48.7]
ds=[round(100-l,1) for l in ls]
colors=['#984ea3','#8e4da2', '#844ca1','#66489f','#49459d','#21409a']

scale=1
max_list=[]
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    max_list.append(np.max(df[df.columns[1]]))

maximum=np.max(max_list)
count=0
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    #plt.loglog(df[df.columns[0]]*10,df[df.columns[1]]*scale/maximum,label=sys.argv[i])
    plt.loglog(df[df.columns[0]]*10,df[df.columns[1]]*scale/maximum,label="{}  :  {}".format(ls[i-1],ds[i-1]),color=colors[i-1])
    scale*=10
    count+=1

plt.xlim(df[df.columns[0]].values[0]*10,5)
plt.legend(title="%L  :  %D",fontsize=16)
plt.ylim(10e-6,10e6)
plt.xlabel('q [nm$^{-1}$]')
plt.ylabel('Intensity [arb. units]')
plt.show()
