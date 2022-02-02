#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fit Spirou Template Quality factors to get a relation with spectral type
Created on Thu Nov 26 15:22:01 2020
@author: nolangrieves
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

vsini_cut = 2.5
templow = 2500
temphigh = 5000
polydegree = 2

output_temps = np.array([5000,4500,4000,3900,3700,3600,3400,3200,3000,2800,2700,2600,2500])
#output_temps = np.array([4500,4000,3900,3700,3600,3400,3200,3000,2800,2700,2600,2500])

#The Template_s1d_GJ9827_sc1d_v_file_AB.fits has the outlier Y-band Q value at 4340

bandpass = 'eniric' # 'CFHT' or 'eniric'
input_file = 'spirou_template_Qvalues_'+bandpass+'-bandpass.txt'
output_file = 'spirou_fit_Qvalues_'+bandpass+'-bandpass.txt'

spirou_results = pd.read_csv(input_file, sep=',', header=0)
#spirou_results = pd.read_csv('spriou_template_Qvalues.txt', sep=',', header=0)
#spirou_results.info()
spirou_results = spirou_results.sort_values(by = ['Teff'])
spirou_results = spirou_results[spirou_results['vsini'] < vsini_cut] 
print('mean spirou template vsinis: %.2f'%np.mean(spirou_results['vsini']))
print('# stars: ',len(spirou_results['vsini']))
      
plt.figure(figsize=(13,10))
plt.plot(spirou_results['Teff'],spirou_results['Q_Y'],'--b.',
         label='Y Spirou Templates')
plt.plot(spirou_results['Teff'],spirou_results['Q_J'],'--g.',
         label='J Spirou Templates')
plt.plot(spirou_results['Teff'],spirou_results['Q_H'],'--r.',
         label='H Spirou Templates')
plt.title('Spirou Template Quality Factors',fontsize=20)
plt.xlabel('Effective Temperature [K]',fontsize=20)
plt.ylabel('Quality Factor',fontsize=20)
plt.text(3500,6000,
         '%i stars with mean vsini=%.2f km/s'%(len(spirou_results['vsini']),np.mean(spirou_results['vsini'])),
         fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

y_fit = np.polyfit(spirou_results['Teff'],spirou_results['Q_Y'],polydegree)
py = np.poly1d(y_fit)
j_fit = np.polyfit(spirou_results['Teff'],spirou_results['Q_J'],polydegree)
pj = np.poly1d(j_fit)
h_fit = np.polyfit(spirou_results['Teff'],spirou_results['Q_H'],polydegree)
ph = np.poly1d(h_fit)
xp = np.linspace(templow,temphigh,100000)
plt.plot(xp,py(xp),'b')
plt.plot(xp,pj(xp),'g')
plt.plot(xp,ph(xp),'r')



y_fit_qt = []
j_fit_qt = []
h_fit_qt = []
for itemp in output_temps:
    #print(itemp)
    idx = (np.abs(xp - itemp)).argmin()
    #print(xp[idx])
    y_fit_qt.append(py(xp[idx]))
    j_fit_qt.append(pj(xp[idx]))
    h_fit_qt.append(ph(xp[idx]))

text_file = open(output_file, "w")
#text_file = open("spriou_fit_qvalues.txt", "w")
print('Teff,Q_Y,Q_J,Q_H')
text_file.write('Teff,Q_Y,Q_J,Q_H\n')
for i in np.arange(len(output_temps)):
    print('%i,%.1f,%.1f,%.1f'%(output_temps[i],y_fit_qt[i],j_fit_qt[i],h_fit_qt[i]))
    text_file.write('%i,%.1f,%.1f,%.1f\n'%(output_temps[i],y_fit_qt[i],j_fit_qt[i],h_fit_qt[i]))
text_file.close()