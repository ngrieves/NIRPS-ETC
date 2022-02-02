#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  9 20:52:18 2022

@author: nolangrieves
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

eff_old = pd.read_csv('old_files/effs.txt',header=0,sep=r"\s+",skiprows=[1])

eff_new = pd.read_csv('updated_files/NIRPS_effs.txt',header=0,sep=r"\s+",skiprows=[1])


for i in range(len(eff_old.columns)):
    plt.figure()
    oldcol = eff_old.columns[i]
    newcol = eff_new.columns[i]
    plt.title(oldcol)
    plt.plot(eff_old['#Lambda'],eff_old[oldcol],'r.',label='old')
    plt.plot(eff_new['#wavelength'],eff_new[newcol],'k.',label='new')
    plt.legend()                 
