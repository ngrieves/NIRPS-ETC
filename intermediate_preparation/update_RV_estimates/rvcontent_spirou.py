#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get RV content from high resolution combined Spirou spectra (templates)
Created on Thu Nov 26 14:38:14 2020
@author: nolangrieves
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.misc import derivative

def deltavrms(spe,wave,sigdet):

    "Compute photon noise uncertainty of one order"
  
    tot=0.
    nx=len(wave)
    for i in range(1,nx-1):
        arg1=wave[i]/(wave[i+1]-wave[i-1])
        arg2=(spe[i+1]-spe[i-1])/(spe[i]+sigdet**2) # different eniric
        #if np.isfinite(spe[i]):
        tot=tot+(spe[i]+sigdet**2)*((arg1*arg2)**2)
    dvrms=3.e8/np.sqrt(tot)
    qval = np.sqrt(tot)/np.sqrt(sum(spe))
    return dvrms,qval

#precision for the photon noise uncertainty equation
#v/c --> 100000/3e8 = 0.00033
#v/c --> 70000/3e8 = 0.00023
precision = 0.00023
c = 299792458.0
sampling = 4.3
snr_pix = 100.0

bandpass = 'eniric' #'eniric' OR 'CFHT'

showplots = True

#sampling of 1 km/s - similar 1 extracted pixel of nirps
#10000 photoelectron --> 100 snr per extracted pixel
#snr_res = snr_pix*(4.3**(1/2))
#snr_pix = snr_res/((number of pixels per bin of resolution)**(1/2))

snr_res = snr_pix*(4.3**(1.0/2.0))
flux_level = snr_pix**2.0
print('Assume SNR = %.1f/pixel (/S1D bin in the case of Spirou templates) or a flux level = %.1f/pixel'%(snr_pix,flux_level))
print('this corresponds to SNR = %.1f/resolution_element given a sampling of %.1f'%(snr_res,sampling))


if bandpass == 'CFHT':
    #### CFHT/WIRCAM bandpasses ###
    #Y-band
    wlnmin_y = 938.6
    wlnmax_y = 1113.4
    #J-band
    wlnmin_j = 1153.5
    wlnmax_j = 1354.4
    #H-band
    wlnmin_h = 1462.8
    wlnmax_h = 1808.5
    #K-band
    wlnmin_k = 1957.7
    wlnmax_k = 2343.1

elif bandpass == 'eniric':
    ### Eniric bandpasses ###
    #Y-band
    wlnmin_y = 1000.0
    wlnmax_y = 1100.0
    #J-band
    wlnmin_j = 1170.0
    wlnmax_j = 1330.0
    #H-band
    wlnmin_h = 1500.0
    wlnmax_h = 1750.0
    #K-band
    wlnmin_k = 2070.0
    wlnmax_k = 2350.0

else:
    print('!!!!!!!!! undefined wavelength bandpass !!!!!!!!!!!')
    wlnmin_y = 0
    wlnmax_y = 0
    wlnmin_j = 0
    wlnmax_j = 0
    wlnmin_h = 0
    wlnmax_h = 0
    wlnmin_k = 0
    wlnmax_k = 0

print('Y band: '+str(wlnmin_y)+' - '+str(wlnmax_y)+' nm')
print('J band: '+str(wlnmin_j)+' - '+str(wlnmax_j)+' nm')
print('H band: '+str(wlnmin_h)+' - '+str(wlnmax_h)+' nm')
print('K band: '+str(wlnmin_k)+' - '+str(wlnmax_k)+' nm')


template_list = pd.read_csv('list_templates_SPIROU.csv', sep=',', header=0)
#template_list.info()
filenames = template_list['File Name']
teffs = template_list['Teff']
spectypes = template_list['spectral type']
vsinis = template_list['vsini']
jmags = template_list['Jmag']


text_file = open('spirou_template_Qvalues_'+bandpass+'-bandpass.txt', "w")
#text_file = open("spriou_template_Qvalues.txt", "w")
print('File,Teff,SpectralType,vsini,Jmag,Q_FullSpec,Q_Y,Q_J,Q_H,Q_K')
text_file.write('File,Teff,SpectralType,vsini,Jmag,Q_FullSpec,Q_Y,Q_J,Q_H,Q_K\n')
        
for ifile in np.arange(len(filenames)):
    file = str(filenames[ifile])
    openfile = fits.open('spirou_templates/'+file)
    #openfile.info()
    objhead = openfile[0].header
    datahead = openfile[1].header
    data = openfile[1].data

    wln = data['wavelength']
    flx = data['flux']
    #get rid of nans
    wln = wln[np.isfinite(flx)]
    flx = flx[np.isfinite(flx)]
    wlno = wln.copy()
    flxo = flx.copy()
    npix = len(flxo)
    
    if showplots:
        #plot each template
        plt.figure(figsize=(13,8))
        plt.plot(wln,flx,'k')
        plt.title(file+' SPIRou Template',fontsize=20)
        plt.xlabel('Wavelength [nm]',fontsize=20)
        plt.ylabel('Normalized Flux',fontsize=20)
        plt.ylim(0,np.nanmax(flx)+0.02)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.axvline(wlnmin_y,color='blue')
        plt.axvline(wlnmax_y,color='blue')
        plt.text(1020,0.2,'Y',color='blue',fontsize=26)
        plt.axvline(wlnmin_j,color='green')
        plt.axvline(wlnmax_j,color='green')
        plt.text(1220,0.2,'J',color='green',fontsize=26)
        plt.axvline(wlnmin_h,color='orange')
        plt.axvline(wlnmax_h,color='orange')
        plt.text(1630,0.2,'H',color='orange',fontsize=26)
        plt.axvline(wlnmin_k,color='red')
        plt.axvline(wlnmax_k,color='red')
        plt.text(2190,0.2,'K',color='red',fontsize=26)
        plt.tight_layout()

    dvrms = deltavrms(flxo,wlno,precision)
    q_all = dvrms[1].copy()
    sigrv_all = c/(q_all*np.sqrt(npix*flux_level))

    #Y-band
    flx = flxo[(wlno > wlnmin_y) & (wlno < wlnmax_y)]
    wln = wlno[(wlno > wlnmin_y) & (wlno < wlnmax_y)]
    dvrms = deltavrms(flx,wln,precision)
    q_y = dvrms[1].copy()
    npix_y = len(wln)
    sigrv_y = c/(q_y*np.sqrt(npix_y*flux_level))

    #J-band
    flx = flxo[(wlno > wlnmin_j) & (wlno < wlnmax_j)]
    wln = wlno[(wlno > wlnmin_j) & (wlno < wlnmax_j)]
    dvrms = deltavrms(flx,wln,precision)
    q_j = dvrms[1].copy()
    npix_j = len(wln)
    sigrv_j = c/(q_j*np.sqrt(npix_j*flux_level))

    #H-band
    flx = flxo[(wlno > wlnmin_h) & (wlno < wlnmax_h)]
    wln = wlno[(wlno > wlnmin_h) & (wlno < wlnmax_h)]
    dvrms = deltavrms(flx,wln,precision)
    q_h = dvrms[1].copy()
    npix_h = len(wln)
    sigrv_h = c/(q_h*np.sqrt(npix_h*flux_level))

    #K-band
    flx = flxo[(wlno > wlnmin_k) & (wlno < wlnmax_k)]
    wln = wlno[(wlno > wlnmin_k) & (wlno < wlnmax_k)]
    dvrms = deltavrms(flx,wln,precision)
    q_k = dvrms[1].copy()
    npix_k = len(wln)
    sigrv_k = c/(q_k*np.sqrt(npix_k*flux_level))

    ##check if they all add up
    combine_yjh = 1/np.sqrt(sigrv_y**-2+sigrv_j**-2+sigrv_h**-2)
    combine_yjhk = 1/np.sqrt(sigrv_y**-2+sigrv_j**-2+sigrv_h**-2+sigrv_k**-2)
    combine_yjh_100 = combine_yjh/1.38 #from Artigaue 2018 
    combine_yjhk_100 = combine_yjhk/1.38 #from Artigaue 2018 

    print('%s,%.1f,%s,%.2f,%.2f,%.1f,%.1f,%.1f,%.1f,%.1f'%(file,teffs[ifile],spectypes[ifile],vsinis[ifile],jmags[ifile],q_all,q_y,q_j,q_h,q_k))
    text_file.write('%s,%.1f,%s,%.2f,%.2f,%.1f,%.1f,%.1f,%.1f,%.1f\n'%(file,teffs[ifile],spectypes[ifile],vsinis[ifile],jmags[ifile],q_all,q_y,q_j,q_h,q_k))

text_file.close()