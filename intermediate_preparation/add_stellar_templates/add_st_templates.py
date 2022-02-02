#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 22:30:27 2022

@author: nolangrieves
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from astropy.io import fits


sts = pd.read_csv('ST_templates.txt',skiprows=[1],sep=r"\s+",header=0)
print(sts.columns)
wlns = sts['#Lambda']

plt.figure(figsize=(14,12))
stars = ['f0v', 'f5v', 'g0v', 'g5v', 'k0v', 'k3v', 'k7v', 'm0v',
       'm1v', 'm2v', 'm3v', 'm4v', 'm5v', 'm6v', 'm7v', 'm8v', 'm9v', 'l1v',
       'l2v', 'l3v', 'l5v', 'l6v', 'l8v', 't2v']
for col in stars:           
    plt.plot(sts['#Lambda'],sts[col],label=col)

#wavelength in [nm]
#flux in [Wm-2um-1]

################################################################
################################################################
# Get G8V star from IRTF:
#http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/Data/G8V_HD75732.txt
g8v = pd.read_csv('G8V_HD75732.txt',skiprows=47,sep=r"\s+",header=None)

g8v_wln = g8v[0]
g8v_flx = g8v[1]
g8v_wln = g8v_wln[g8v_flx != -999]*1000 #um to nm
g8v_flx = g8v_flx[g8v_flx != -999]

#interpolate to wln range
fint = interp1d(g8v_wln, g8v_flx, kind='cubic',fill_value="extrapolate")
g8v_flxint = fint(wlns)

plt.plot(sts['#Lambda'],g8v_flxint,'--',color='red',label='g8v added')
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel('Flux [Wm-2um-1]',fontsize=18)
plt.xlabel('Wavelength [nm]',fontsize=18)
plt.legend()

#add column
sts['g8v'] = g8v_flxint
                       
#stso = sts.reindex(columns=['#Lambda', 'f0v', 'f5v', 'g0v', 'g5v','g8v','k0v', 'k3v', 'k7v', 'm0v',
#       'm1v', 'm2v', 'm3v', 'm4v', 'm5v', 'm6v', 'm7v', 'm8v', 'm9v', 'l1v',
#       'l2v', 'l3v', 'l5v', 'l6v', 'l8v', 't2v'])
#print(stso.columns)

print(sts.columns)

sts.to_csv('STAR_templates.txt', header=sts.columns,index=None, sep='\t')
################################################################
################################################################

             
             
######## A Star Template Work Below ########
             
#             
#################################################################
#################################################################          
##Get A star from CRIRES+ spectrophotometric standard stars
##https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar
#             
##A0V    9700
##A9V    7400 
#
##HIP104139 = HD104139, Teff=9.200K, SpT=A1V, K4.1mag
##HR7950 : A1 V
##HR3454 Teff=17000 K and logg=4.0
##HR4468 Teff=10500 K and logg=4.0
##HR5501 : B9.5V 
##HR7950 : A1 V 
##HR8634 : B8 V Teff = 11000 K and logg=4.0
#
#crires_stars = ['HIP104139','HR7950','HR5501','HR7950','HR8634']
#
#sed_HIP104139 = pd.read_csv('HIP104139.txt',comment='#',sep=r"\s+",header=None)
#wln_HIP104139 = sed_HIP104139[0]*10.0 # lambda [nm] to Angstrom
#flx_HIP104139 = sed_HIP104139[1] #F(nu) [Jy]
##convert units 
## F(nu) [Jy] to erg/s/cm2/A
## Flux Density conversion: F_lambda [erg/s/cm2/A] = 3.0x10^-5 * F_nu / wln[A]^2
#flx_lambda_HIP104139 = flx_HIP104139 * 3.0e-5 / wln_HIP104139**2
##divide by 10000 for Angstrom to um, multiply by 10000 for m-2 to cm-2 (so units the same)
##1 Watt = 10000000 erg/s       so / (1e7)
#flx_lambda_HIP104139 = flx_lambda_HIP104139 / 1e7
#wln_HIP104139 = wln_HIP104139 / 10.0 #Angstrom to nm
#
#
#
#sed_HR3454 = pd.read_csv('HR3454.txt',comment='#',sep=r"\s+",header=None)
#wln_HR3454 = sed_HR3454[0]*10.0 # lambda [nm] to Angstrom
#flx_HR3454 = sed_HR3454[1] #F(nu) [Jy]
##convert units 
## F(nu) [Jy] to erg/s/cm2/A
## Flux Density conversion: F_lambda [erg/s/cm2/A] = 3.0x10^-5 * F_nu / wln[A]^2
##From: https://www.stsci.edu/~strolger/docs/UNITS.txt
##[Y erg/cm^2/s/A]             = 2.99792458E-05 * [X1 Jy] / [X2 A]^2
#
#
#flx_lambda_HR3454 = flx_HR3454 * 3.0e-5 / wln_HR3454**2
##divide by 10000 for Angstrom to um, multiply by 10000 for m-2 to cm-2 (so units the same)
##1 Watt = 10000000 (1e7) erg/s          so / (1e7)
#flx_lambda_HR3454 = flx_lambda_HR3454 / 1e7
#wln_HR3454 = wln_HR3454 / 10.0 #Angstrom to nm
#
#plt.plot(wln_HIP104139,flx_lambda_HIP104139,color='black',label='HIP104139 A1V')
#plt.plot(wln_HR3454,flx_lambda_HR3454,color='green',label='HR3454')
#################################################################
#################################################################
#
#
#################################################################
#################################################################
##### HIGH RESOLUTION spectra from goettingen
##http://phoenix.astro.physik.uni-goettingen.de/
##Husser 2013:
##https://ui.adsabs.harvard.edu/abs/2013A%26A...553A...6H/abstract
#hdu_wln = fits.open('goettingen/highres/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')
#data_wln = hdu_wln[0].data /10.0
#
#hdu_8000 = fits.open('goettingen/highres/lte08000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')
#data_8000 = hdu_8000[0].data #'erg/s/cm^2/cm'  
##divide by 10000 for cm-2 to m-2        so / (1e4)
##divide by 10000 for cm-1 to um-1       so / (1e4)
##1 Watt = 10000000 (1e7) erg/s          so / (1e7)
#data_8000 = data_8000 / 1e15
#
#hdu_9800 = fits.open('goettingen/highres/lte09800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits')
#data_9800 = hdu_9800[0].data / 1e15
#
##plt.plot(data_wln,data_8000,label='Teff=8000')
##plt.plot(data_wln,data_9800,label='Teff=9800')
#
#################################################################
#################################################################
#
#plt.xlim(900,2000)
#plt.xticks(fontsize=16)
#plt.yticks(fontsize=16)
#plt.ylabel('Flux [Wm-2um-1]',fontsize=18)
#plt.xlabel('Wavelength [nm]',fontsize=18)
#plt.legend()