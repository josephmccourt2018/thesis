import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import os
import fnmatch
import re

#first have to load in Absolute.dat which was generated by Absolute_water_intensity.m\
absolute = pd.read_table('Absolute.dat',sep = '\s+',names = ["Temperature (K)", "Absolute Intensity (cm$^{-1}$)"],header = None)
#plt.plot(absolute[absolute.columns[0]],absolute[absolute.columns[1]])

folder = input("folder name (e.g. rack1): ");

choice = 1; # choice =1 implies water is the solvent. choice=2 implies water is not the solvent

Temp = 20; # Temperature in degree celcius. This is required for absolute intensity calibration

Vf = float(input("Volume fraction (e.g. 0.002): ")); # Volume fraction

smear = 2; # smear =1 implies create low resolution data. smear =2 implies no need for this step

region = 2; # This choice is only applicable if smear = 1 (smearing/low resolution data option is 'on')
            # region = 1 implies smearing will be performed in higher q
            # regions of SAXS data. This is useful when you don't want to
            # smear out the sharp diffraction peaks in the WAXS or low q
            # region SAXS. region = 2 is for smearing the whole of SAXS
            # data. This will work well if there are no diffraction peaks
            # in q-range covered by SAXS detector. region = 3 implies that
            # the smearing will be performed over whole of the
            # SAXS/MAXS/WAXS data. More smearing options can be added as
            # needed.
            
smear_range = 0.06; # This applies to smear = 1 && region = 1 choices. 
              # This implies that the smearing will be performed over a region q_saxs(n_saxs)-range, 
              #where n_saxs is the total number of points in the saxs data.

#read in data files

file_prefix = folder + "/plot_files/mean_plot_files/";

date = 0 # declare variable for date that will be extracted from *.gpt filename
         # for 'empt' files

if choice==1:
    
    empt = input('scan number for empty: ');
    
    wat = input('scan number for water: ');
    
    samp = input('scan number for sample: ');
      
    string=0 
    for file in os.listdir(file_prefix):
        if fnmatch.fnmatch(file, str(empt)+'*.gpt'): 
        #if fnmatch.fnmatch(file, '*'+str(59)+'*.gpt'): 
            string = file
            break;
    date= re.search(r'_(\d{8})_',string).group(1)
    print("Date of Experiment: " + date)
    
    
    filename1 = file_prefix + folder + '_hs104_' + empt + '.txt';
    
    filename2 = file_prefix + folder + '_hs104_' + wat + '.txt';
    
    filename3 = file_prefix + folder + '_hs104_' + samp + '.txt';
    
    filename4 = file_prefix + folder + '_hs103_' + empt + '.txt';
    
    filename5 = file_prefix + folder + '_hs103_' + wat + '.txt';
    
    filename6 = file_prefix + folder + '_hs103_' + samp + '.txt';
    
    filename7 = file_prefix + folder + '_hs102_' + empt + '.txt';
    
    filename8 = file_prefix + folder + '_hs102_' + wat + '.txt';
    
    filename9 = file_prefix + folder + '_hs102_' + samp + '.txt';
    
    empt_saxs = pd.read_csv(filename1,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    wat_saxs = pd.read_csv(filename2,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    samp_saxs = pd.read_csv(filename3,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    empt_maxs = pd.read_csv(filename4,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    wat_maxs = pd.read_csv(filename5,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    samp_maxs = pd.read_csv(filename6,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    empt_waxs = pd.read_csv(filename7,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    wat_waxs = pd.read_csv(filename8,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    samp_waxs = pd.read_csv(filename9,sep = '\s+',comment='#',names = ["q", "Intensity", "sigma"],header = None);
    
    sol_saxs = wat_saxs;
    
    sol_maxs = wat_maxs;
    
    sol_waxs = wat_waxs;


#________________________________________________________________________________________________________________________________________
#The absolute intensity for water at q=0. This requires the file Absolute.dat which is created by Absolute_intensity_water.m
#________________________________________________________________________________________________________________________________________

T = absolute[absolute.columns[0]];

Abs_wat = absolute[absolute.columns[1]];

m = len(T);

T_new = T-(Temp + 273.16);

T_new = abs(T_new);

I = min(T_new); #index

R = Temp%10;

if R==0:
    
    Abs = Abs_wat[I];
     
elif R<=5:
    
    Abs = (1-R/10)*Abs_wat(I) + R/10*Abs_wat(I+1);
    
else:
    Abs = (1-R/10)*Abs_wat(I) + R/10*Abs_wat(I-1);

# _____________________________________________________________________
# Other Inputs
#_____________________________________________________________________
      
q_saxs_i = 0.02; # This is the starting q value for fitting the water+cap - cap dat for intensity calibration

q_saxs_f = 0.06; # This is the end q value for fitting the water+cap - cap dat for intensity calibration
 
Energy = input('Energy in keV: ');

Lambda = 12.398/float(Energy); #wavelength. This is only required in the smearing option.

#_________________________________________________________________________
#Defining q
#_________________________________________________________________________

q_saxs = empt_saxs['q'];

q_maxs = empt_maxs['q'];
 
q_waxs = empt_waxs['q'];
 
n_saxs = len(q_saxs);
 
n_maxs = len(q_maxs);
 
n_waxs = len(q_waxs);

#  ____________________________________________________________________________________________________________________________________________
# This part of the script uses the SAXS data from water+capillary and empty
# capillary to calibrate the intensities
# ________________________________________________________________________________________________________________________________________________
 
water = wat_saxs['Intensity']-empt_saxs['Intensity'];

err = np.sqrt((wat_saxs['sigma']**2 + (empt_saxs['sigma']**2)));

# plot data

plt.figure(1)

plt.errorbar(q_saxs, water, err)

plt.title('scattering from water')

plt.xlabel('q ($\AA ^{-1}$)');

plt.ylabel('Int (cm$^{-1}$)');    

start = 0;

while q_saxs[start] <= q_saxs_i:
    start = start + 1;

finish = 0;

while q_saxs[finish] <= q_saxs_f:
    finish = finish + 1;

X=[q_saxs[r] for r in np.arange(start,finish)]
Y=[water[r] for r in np.arange(start,finish)]
error=[err[r] for r in np.arange(start,finish)]

#  plot to test if it is chosing the desired range and do the linear fitting

plt.figure(2);

plt.errorbar(X, Y, error);

P = np.polyfit(X,Y,1);

m = len(X); 

Yfit = [P[0]*x + P[1] for x in X]

plt.plot(X, Yfit)

plt.title('Intensity Calibration')

plt.xlabel('q ($\AA ^{-1}$)');

plt.ylabel('Int (cm$^{-1}$)');   

scale_factor = Abs/P[1] # this scale factor will be used in the next section

#________________________________________________________________________________________________________________________________________
#This part of the script is used to obtain solvent subtracted data
#________________________________________________________________________________________________________________________________________

saxs = (samp_saxs['Intensity'] - empt_saxs['Intensity']*Vf - sol_saxs['Intensity']*(1-Vf))*scale_factor; # scale factor is obtained in the firt part of the script above
err_saxs = (np.sqrt((samp_saxs['sigma'])**2 + (1-Vf)**2*(sol_saxs['sigma']**2) + Vf**2*empt_saxs['sigma']**2))*scale_factor;

maxs = (samp_maxs['Intensity'] - empt_maxs['Intensity']*Vf - sol_maxs['Intensity']*(1-Vf))*scale_factor;
err_maxs = (np.sqrt((samp_maxs['sigma'])**2 + (1-Vf)**2*(sol_maxs['sigma']**2) + Vf**2*empt_maxs['sigma']**2))*scale_factor;

waxs = (samp_waxs['Intensity'] - empt_waxs['Intensity']*Vf - sol_waxs['Intensity']*(1-Vf))*scale_factor;
err_waxs = (np.sqrt((samp_waxs['sigma'])**2 + (1-Vf)**2*(sol_waxs['sigma']**2) + Vf**2*empt_waxs['sigma']**2))*scale_factor;

# plot data figure 3, 4, 5 

plt.figure(3)
plt.errorbar(q_saxs,saxs,err_saxs)
plt.title('background subtracted SAXS data');
plt.xlabel('q ($\AA ^{-1}$)');
plt.ylabel('Int (cm$^{-1}$)');   

plt.figure(4)
plt.errorbar(q_maxs,maxs,err_maxs);
plt.title('background subtracted MAXS data');
plt.xlabel('q ($\AA ^{-1}$)');
plt.ylabel('Int (cm$^{-1}$)');   

plt.figure(5)
plt.errorbar(q_waxs,waxs,err_waxs)
plt.title('background subtracted WAXS data');
plt.xlabel('q ($\AA ^{-1}$)');
plt.ylabel('Int (cm$^{-1}$)');   

#________________________________________________________________________________________________________________________________________

# Here the offset between the SAXS data and the MAXS+WAXS data is estimated
#_________________________________________________________________________________________________________________________________________

q_saxs_end = q_saxs[n_saxs-1];

q_maxs_begin = q_maxs[0];

q_waxs_begin = q_waxs[0];

start = 0;

finish = 0;
for k in range(0,n_saxs):
    while q_saxs[start] < q_maxs_begin:  
          start = start+1;
            
for k in range(0,n_maxs):
    while q_maxs[finish] < q_saxs_end:      
          finish = finish+1;

X1 = [q_saxs[k] for k in range(start,n_saxs)]

Y1 = [saxs[k] for k in range(start,n_saxs)]

#trapz was giving incorrect values for the integral, so just did my own rectangular integration
#^fixed this issue

Int1 = trapz(Y1,X1)

#Int1 = 0;
#for i in range(0,len(X1)-1):
#    Int1 +=(X1[i+1]-X1[i])*Y1[i]

X2 = [q_maxs[k] for k in range(0,finish)]

Y2 = [maxs[k] for k in range(0,finish)]

Int2 = trapz(Y2,X2)
#Int2 = 0;
#for i in range(0,len(X2)-1):
#    Int2 +=(X2[i+1]-X2[i])*Y2[i]

offset= (Int1-Int2)/(X2[finish-1]-X2[0])

q = []
I = []
error = []
for k in range(0,n_saxs+n_maxs+n_waxs):

    if k < n_saxs:
        if q_saxs[k] < q_maxs_begin:
            q.append(q_saxs[k]);
            #I.append(saxs[k]-offset);
            I.append(saxs[k])
            error.append(err_saxs[k]);
        else:
            continue;
     
    elif k >= n_saxs and k < n_saxs + n_maxs:
        if q_maxs[k - n_saxs] < q_waxs_begin:
            q.append(q_maxs[k - n_saxs]);
            I.append(maxs[k - n_saxs]);
            error.append(err_maxs[k - n_saxs]);
        else:
            continue;
        
    elif k >= n_saxs + n_maxs and k < n_saxs + n_maxs + n_waxs:
        q.append(q_waxs[k - n_saxs - n_maxs]);
        I.append(waxs[k - n_saxs - n_maxs]);
        error.append(err_waxs[k - n_saxs - n_maxs]);

df = pd.DataFrame(data={'q':q, 'I':I, 'error':error})
plt.figure(6)
plt.xscale('log')
plt.yscale('log')
plt.plot(df['q'],df['I'])


#_____________________________________________________________________________________________________
# Filter out low intensity points and reset index numbering, then plot
#_____________________________________________________________________________________________________

df = df[df['I']>1e-4]
df.reset_index(drop=True)
plt.figure(7)
plt.xscale('log')
plt.yscale('log')
plt.plot(df['q'],df['I'])


#_____________________________________________________________________________________________________
#save
#_____________________________________________________________________________________________________

info = input("Sample Info: "); # information about sample added to header

filename='sample' + samp + '.csv';
print(filename)
with open(filename, 'w') as f:
    f.write('# Volume Fraction = {}\n'.format(Vf))
    f.write('# Energy [keV]  = {}\n'.format(Energy))
    f.write('# Sample Info: '+info+'\n')
    f.write('# Date of Experiment: '+date+'\n')
    df.to_csv(f)


plt.show()
