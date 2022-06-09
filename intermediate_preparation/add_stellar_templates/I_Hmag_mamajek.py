#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 17:00:47 2022

@author: nolangrieves
"""

import pandas as pd
import numpy as np


table = pd.read_csv('EEM_dwarf_UBVIJHK_colors_Teff.txt',sep=r"\s+",header=22)

table[table['#SpT'] == 'B3V']


#I = H+(Ic-H)
#Ic-H = (Mv – V-Ic) – (H-Ks + Ks)

SpTs = ['B3V','B8V','B9V','A1V']
for SpT in SpTs:
    Ic_H = float(table[table['#SpT'] == SpT]['Mv']) - float(table[table['#SpT'] == SpT]['V-Ic']) - (float(table[table['#SpT'] == SpT]['H-Ks']) + float(table[table['#SpT'] == SpT]['M_Ks']))
    print(SpT,': I-H =',Ic_H)



#print(table[table['#SpT'] == SpT]['Mv'])
#print(table[table['#SpT'] == SpT]['V-Ic'])
#print(table[table['#SpT'] == SpT]['H-Ks'])
#print(table[table['#SpT'] == SpT]['M_Ks'])

#For F0V:
#Mv = 2.51
#V-Ic = 0.339
#H-Ks = 0.045
#Ks = 1.76
#Ic-H = (2.51 – 0.339) – (0.045 + 1.76) = 0.366


#For G8V
#Mv = 5.32
#V-Ic = 0.786
#H-Ks = 0.082
#Ks = 3.55
#Ic-H = (5.32 – 0.786) – (0.082 + 3.55) = 0.902
