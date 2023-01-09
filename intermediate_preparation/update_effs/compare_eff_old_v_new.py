#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 20:52:18 2022

@author: nolangrieves
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

eff_v1 = pd.read_csv('v1_files_simulation/effs.txt',header=0,sep=r"\s+",skiprows=[1])

eff_dec2022 = pd.read_csv('june2022_files/updated_files/NIRPS_effs.txt',header=0,sep=r"\s+",skiprows=[1])

eff_new = pd.read_csv('updated_files_jan2023/NIRPS_effs.txt',header=0,sep=r"\s+",skiprows=[1])

for i in range(len(eff_v1.columns)):
    plt.figure()
    v1col = eff_v1.columns[i]
    eff_dec2022col = eff_dec2022.columns[i]
    newcol = eff_new.columns[i]
    print(v1col,eff_dec2022col,newcol)
    plt.title(v1col)
    print()
    plt.plot(eff_v1['#Lambda'],eff_v1[v1col],'r.',alpha=0.5,label='v1 simiulated')
    plt.plot(eff_dec2022['#wavelength'],eff_dec2022[eff_dec2022col],'k.',alpha=0.5,label='jun2022 grating1')
    plt.plot(eff_new['#wavelength'],eff_new[newcol],'b.',alpha=0.5,label='dec2022 grating2')
    plt.legend()                 
