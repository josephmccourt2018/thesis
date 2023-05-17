import MDAnalysis as mda
import numpy as np
import pandas as pd
from tqdm import tqdm

u = mda.Universe('../md_semi.tpr','../combined.xtc')

molecules = 90

apl_list = []

for ts in tqdm(u.trajectory):
    apl = u.dimensions[0]*u.dimensions[1]/molecules
    apl_list.append([u.trajectory.time, apl])

df = pd.DataFrame(apl_list,columns=['Time (ps)','APL (Ang^2/lipid)'])

df.to_csv("apl_vs_time.csv")
