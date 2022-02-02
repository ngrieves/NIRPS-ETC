#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Phoenix Quality Factor comparison for different resolutions
Created on Thu Nov 26 16:14:48 2020
@author: nolangrieves
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

en_res = '70k'
vsini = 1.0
sampling = 4.3

c=299792458

output_temps = np.array([5000,4500,4000,3900,3700,3600,3400,3200,3000,2800,2700,2600,2500])

bandpass = 'CFHT' # 'CFHT' or 'eniric'


input_file = 'phoenix_eniric_Qfactors_'+bandpass+'-bandpass.csv'
output_file = 'phoenix_Q_conversions_'+bandpass+'-bandpass.txt'

eniric_results = pd.read_csv(input_file, sep=',', header=0)
#eniric_results.info()

if bandpass == 'eniric':
    y_label = 'Y'
    j_label = 'J'
    h_label = 'H'
if bandpass == 'CFHT':
    y_label = 'Y_CFHT'
    j_label = 'J_CFHT'
    h_label = 'H_CFHT'

en_y = ((eniric_results['sampling'] == sampling) & 
                  (eniric_results['vsini'] == vsini) &
                  (eniric_results['resolution'] == en_res) &
                  (eniric_results['band'] == y_label))

en_j = ((eniric_results['sampling'] == sampling) & 
                  (eniric_results['vsini'] == vsini) &
                  (eniric_results['resolution'] == en_res) &
                  (eniric_results['band'] == j_label))

en_h = ((eniric_results['sampling'] == sampling) & 
                  (eniric_results['vsini'] == vsini) &
                  (eniric_results['resolution'] == en_res) &
                  (eniric_results['band'] == h_label))

en_yval_70k = eniric_results[en_y]
en_jval_70k = eniric_results[en_j]
en_hval_70k = eniric_results[en_h]

plt.figure(figsize=(13,10))
plt.plot(en_yval_70k['temp'],en_yval_70k['quality'],'-b.',
         label='Y Eniric Phoenix '+en_res+' vsini='+str(vsini)+' sampling='+str(sampling))
plt.plot(en_jval_70k['temp'],en_jval_70k['quality'],'-g.',
         label='J Eniric Phoenix '+en_res+' vsini='+str(vsini)+' sampling='+str(sampling))
plt.plot(en_hval_70k['temp'],en_hval_70k['quality'],'-r.',
         label='H Eniric Phoenix '+en_res+' vsini='+str(vsini)+' sampling='+str(sampling))

en_res = '80k'
en_y = ((eniric_results['sampling'] == sampling) & 
                  (eniric_results['vsini'] == vsini) &
                  (eniric_results['resolution'] == en_res) &
                  (eniric_results['band'] == y_label))

en_j = ((eniric_results['sampling'] == sampling) & 
                  (eniric_results['vsini'] == vsini) &
                  (eniric_results['resolution'] == en_res) &
                  (eniric_results['band'] == j_label))

en_h = ((eniric_results['sampling'] == sampling) & 
                  (eniric_results['vsini'] == vsini) &
                  (eniric_results['resolution'] == en_res) &
                  (eniric_results['band'] == h_label))

en_yval_xk = eniric_results[en_y]
en_jval_xk = eniric_results[en_j]
en_hval_xk = eniric_results[en_h]
plt.plot(en_yval_xk['temp'],en_yval_xk['quality'],'--b*',
         label='Y Eniric Phoenix '+en_res+' vsini='+str(vsini)+' sampling='+str(sampling))
plt.plot(en_jval_xk['temp'],en_jval_xk['quality'],'--g*',
         label='J Eniric Phoenix '+en_res+' vsini='+str(vsini)+' sampling='+str(sampling))
plt.plot(en_hval_xk['temp'],en_hval_xk['quality'],'--r*',
         label='H Eniric Phoenix '+en_res+' vsini='+str(vsini)+' sampling='+str(sampling))

plt.title('Pheonix Quality Factor Comparison',fontsize=20)
plt.xlabel('Effective Temperature [K]',fontsize=20)
plt.ylabel('Quality Factor',fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.figure(figsize=(13,10))
plt.plot(en_yval_xk['temp'],np.array(en_yval_70k['quality'])/np.array(en_yval_xk['quality']),'b',label='Y')
plt.plot(en_yval_xk['temp'],np.array(en_jval_70k['quality'])/np.array(en_jval_xk['quality']),'g',label='J')
plt.plot(en_yval_xk['temp'],np.array(en_hval_70k['quality'])/np.array(en_hval_xk['quality']),'r',label='H')
plt.title('Quality Factor Comparison',fontsize=20)
plt.xlabel('Effective Temperature [K]',fontsize=20)
plt.ylabel(('Q 70k / Q '+en_res),fontsize=20)
plt.legend(fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


temps = en_yval_xk['temp'].values
y_conversion = np.array(en_yval_70k['quality'])/np.array(en_yval_xk['quality'])
j_conversion = np.array(en_jval_70k['quality'])/np.array(en_jval_xk['quality'])
h_conversion = np.array(en_hval_70k['quality'])/np.array(en_hval_xk['quality'])

#text_file = open("phoenix_Q_conversions.txt", "w")
text_file = open(output_file, "w")
print('Teff,Q_Y_conv,Q_J_conv,Q_H_conv')
text_file.write('Teff,Q_Y,Q_J,Q_H\n')
for i in np.arange(len(temps)):
    print('%i,%.4f,%.4f,%.4f'%(temps[i],y_conversion[i],j_conversion[i],h_conversion[i]))
    text_file.write('%i,%.4f,%.4f,%.4f\n'%(temps[i],y_conversion[i],j_conversion[i],h_conversion[i]))
text_file.close()

