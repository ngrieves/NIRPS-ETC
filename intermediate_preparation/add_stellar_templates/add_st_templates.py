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
###### Add in G8V star ######
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

plt.plot(sts['#Lambda'],g8v_flxint,':',label='g8v')#,color='red'
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel('Flux [Wm-2um-1]',fontsize=18)
plt.xlabel('Wavelength [nm]',fontsize=18)


#add column
sts['g8v'] = g8v_flxint
                       
#stso = sts.reindex(columns=['#Lambda', 'f0v', 'f5v', 'g0v', 'g5v','g8v','k0v', 'k3v', 'k7v', 'm0v',
#       'm1v', 'm2v', 'm3v', 'm4v', 'm5v', 'm6v', 'm7v', 'm8v', 'm9v', 'l1v',
#       'l2v', 'l3v', 'l5v', 'l6v', 'l8v', 't2v'])
#print(stso.columns)

################################################################
###### Add in A and B stars ######
################################################################
##CRIRES+ spectrophotometric standard stars
##https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar

##HIP104139 = HD104139, Teff=9.200K, SpT=A1V, K4.1mag
##HR7950 : A1 V
##HR3454 Teff=17000 K and logg=4.0
##HR4468 Teff=10500 K and logg=4.0
##HR5501 : B9.5V 
##HR8634 : B8 V Teff = 11000 K and logg=4.0
####Use HR5501 for B9p5V


#Use HR3454 for B3V
#Use HR8634 for B8V
#Use HR4468 for B9V
#Use HR7950 for A1V

crires_stars = ['HR7950','HR3454','HR4468','HR5501','HR8634']
crires_stars_label = ['A1V_HR7950','B3V_HR3454','B9V_HR4468','B9V_HR5501','B8V_HR8634']

crires_stars = ['HR3454','HR8634','HR4468','HR7950']
crires_stars_label = ['b3v','b8v','b9v','a1v']

##### Convert Unites #####
#https://irsa.ipac.caltech.edu/data/SPITZER/docs/spitzermission/missionoverview/spitzertelescopehandbook/18/#_Toc60818984
#Infrared Flux Units
#The infrared flux density from a point source is most commonly given in units of Jansky (Jy) where:
#1 Jy = 10-23 erg s-1 cm-2 Hz-1 = 10-26 Watts m-2 Hz-1 = Fν
#The conversion between Janskys and flux density in Wm-2 per unit wavelength is given by:
#Fν x 10-26 x c/λ2= Fλ

c = 299792458 #m/s

for istar in range(len(crires_stars)):
    star = crires_stars[istar]
    slabel = crires_stars_label[istar]
    sed_star = pd.read_csv(star+'.txt',comment='#',sep=r"\s+",header=None)
    wln_star = sed_star[0] # lambda [nm]
    flx_star_o = sed_star[1] # F(nu) [Jy]
    wln_star_m = wln_star * 1e-9 #nm to m
    
    flx_star = flx_star_o * 1e-26 * (c/wln_star_m**2)
    flx_star = flx_star * 1e-6 #m-1 to micron-1
    
    #interpolate to wln range
    fint = interp1d(wln_star,flx_star, kind='cubic',fill_value="extrapolate")
    star_flxint = fint(wlns)

    plt.plot(wlns,star_flxint,':',label=slabel)
    sts[slabel] = star_flxint






plt.legend()
plt.savefig('NIRPS_ETC_stellar_templates.pdf')

print(sts.columns)

sts.to_csv('STAR_templates.txt', header=sts.columns,index=None, sep='\t')
################################################################
################################################################