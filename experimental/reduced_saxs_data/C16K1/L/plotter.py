import numpy as np
import matplotlib.pyplot as plt
plt.style.use('presentation')
import pandas as pd
import sys
import os

plt.figure()
basedir=os.getcwd()
for i in range(1,len(sys.argv)):
    df = pd.read_csv(basedir+'/'+sys.argv[i],names=['q','I','error'],comment='#',header=1)
    plt.loglog(df[df.columns[0]],df[df.columns[1]],label=sys.argv[i])
    plt.legend()
    plt.xlabel('q [$\AA^{-1}$]')
    plt.ylabel('Intensity [cm$^{-1}$]')
plt.show()
