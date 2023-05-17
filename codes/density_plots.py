import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('veusz')
charge="10"
tail="4"
basedir="/home/joey/Documents/quest_b1021/CK1/MARTINI_2/bilayers/charge_variation/{}/{}/bilayer/4x_bilayer/".format(charge,tail)
tails = pd.read_csv(basedir+"tails.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);
heads = pd.read_csv(basedir+"heads.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);
ions = pd.read_csv(basedir+"ions.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);
water = pd.read_csv(basedir+"pw.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);

plt.figure()
#plt.plot(tails['z'],tails['density'],color='#6b6b6b',label='Tails')
#plt.plot(heads['z'],heads['density'],color='#2273b5',label='Heads')
#plt.plot(ions['z'],ions['density'],color="#2cff2c",label='Cl$^-$')
#plt.plot(water['z'],water['density'],label='PW')
plt.plot(tails['z'],tails['density'],color='#6b6b6b')
plt.plot(heads['z'],heads['density'],color='#2273b5')
plt.plot(ions['z'],ions['density'],color="#2cff2c")
plt.plot(water['z'],water['density'])
plt.plot(0,0,color='k',label=charge)
plt.xlabel("z (nm)")
plt.ylabel("$\\rho$ (kg/m$^3$)")
plt.xticks(np.arange(-10,12,2))
plt.legend()
data_list = [tails,heads,ions,water]

charge="90"
tail="4"
basedir="/home/joey/Documents/quest_b1021/CK1/MARTINI_2/bilayers/charge_variation/{}/{}/bilayer/4x_bilayer/".format(charge,tail)
tails = pd.read_csv(basedir+"tails.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);
heads = pd.read_csv(basedir+"heads.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);
ions = pd.read_csv(basedir+"ions.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);
water = pd.read_csv(basedir+"pw.xvg",sep = '\s+',skiprows=24,names = ["z", "density"],header = None);

#plt.plot(tails['z'],tails['density'],color='#6b6b6b',label='Tails',linestyle='--')
#plt.plot(heads['z'],heads['density'],color='#2273b5',label='Heads',linestyle='--')
#plt.plot(ions['z'],ions['density'],color="#2cff2c",label='Cl$^-$',linestyle='--')
#plt.plot(water['z'],water['density'],label='PW',color='r',linestyle='--')
plt.plot(tails['z'],tails['density'],color='#6b6b6b',linestyle='--')
plt.plot(heads['z'],heads['density'],color='#2273b5',linestyle='--')
plt.plot(ions['z'],ions['density'],color="#2cff2c",linestyle='--')
plt.plot(water['z'],water['density'],color='r',linestyle='--')
plt.plot(0,0,color='k',linestyle='--',label=charge)
plt.xlabel("z (nm)")
plt.ylabel("$\\rho$ (kg/m$^3$)")
plt.xticks(np.arange(-10,12,2))
plt.legend(title="Charge Percent")
plt.title("C16K1")
plt.show()
data_list = [tails,heads,ions,water]
