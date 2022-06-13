 #!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Copyright (C) 2018, Bruno L. Canto Martins
# NIRPS ETC code python
# Last modification - March, 2018.
# 
# New edits 2020-2022 by Nolan Grieves:
# Updated readout noise calculation
# Keep object/template flux units in (W m-2 um-1) 
# New RV precision calculation
# include CFHT or Eniric bandpasses as options
# Update new instrument specifications
#  - transmission/efficiency
#  - resolutions
#  - wavelength/blaze values from in lab
# correct for difference between wavelength array size and pixel bin size for total number of electrons
# New WAVE and BLAZE from La Silla (June 2022)
# Add A adn B stellar templates (June 2022)

# Imports
import numpy
from os import system, chdir
import matplotlib.pyplot as plt
import matplotlib as mpl
#from astropy import constants
from math import *
from nirps_etc_lib import *
from decimal import Decimal
import pandas as pd
import numpy as np
import time
start_time = time.time()


################################################################
##### CHANGE INPUTS HERE ####
#
##Observation Mode (HA/HE): 
#obs_mode= 'HA'
##Seeing, in arcsec (range 0.7-1.2):
#seeing= 1.0
##Airmass (range 1.0-2.0):
#airmass= 1.0
##Object magnitude (H band):
#H = 9
##Exposure time (in sec):
#t_exp = 1200
#
##ST Models for IRFT
##Spectral type that only have SNR estimate:
##B3V/B8V/B9V/A1V/F0V/F5V/G0V/G5V/G8V/K0V/L1V/L2V/L3V/L5V/L6V/L8V/T2V
##Spectral type that have SNR and RV precision estimate: 
##K3V/K7V/M0V/M1V/M2V/M3V/M4V/M5V/M6V/M7V/M8V/M9V
#st = 'M4V'
#
##wavelength bandpasses for YJH ('CFHT' or 'Eniric')
#bandpass = 'CFHT'
#
##Vsini of star for Eniric Pheonix spectra RV calculations will be for [0.1,1.0,5.0,10.0] km/s
#

#Use all order information (TSR) or only FSR
waveselect = 'FSR'
#FSR = Free Spectral Range calculated from 50% of blaze
#TSR = Total Spectral Range
################################################################

input_targets_file = 'etc_targets_input_COMM5.txt'
output_targets_file = 'etc_targets_output_COMM5.txt'

input_targets = pd.read_csv(input_targets_file,sep=r"\s+",header=0)
targets = input_targets['target']
sts = input_targets['st']
obs_modes = input_targets['obs_mode']
seeings = input_targets['seeing']
airmasses = input_targets['airmass'] 
Hmags = input_targets['H']
t_exps = input_targets['t_exp']
bandpasses = input_targets['bandpass']



output_targets = open(output_targets_file, "w")
output_targets.write('target Mean_S/N_(ph/pxl) Mean_S/N_Y(ph/pxl) Mean_S/N_J(ph/pxl) Mean_S/N_H(ph/pxl) Hmax_1619nm(ph/pxl) spRV enRV_vsini01 enRV_vsini1 enRV_vsini5 enRV_vsini10 \n')


if waveselect == 'FSR':
    effs_file = 'NIRPS_effs_FSR.txt'
    wave_range_file = 'NIRPS_wave_range_FSR.txt'
    tapas_file = 'NIRPS_tapas_FSR.txt'
    st_templates_file = 'NIRPS_STAR_templates_FSR.txt'
else:
    effs_file = 'NIRPS_effs.txt'
    wave_range_file = 'NIRPS_wave_range.txt'
    tapas_file = 'NIRPS_tapas.txt'
    st_templates_file = 'NIRPS_STAR_templates.txt'


for itarget in range(len(targets)):
    target = targets[itarget]
    obs_mode= obs_modes[itarget]
    seeing = seeings[itarget]
    airmass = airmasses[itarget]
    H = Hmags[itarget]
    t_exp = t_exps[itarget]
    st = sts[itarget]
    bandpass = bandpasses[itarget]


    if bandpass == 'CFHT':
    	spriou_fit_qvalues_file = 'spirou_fit_Qvalues_CFHT-bandpass.txt'
    	phoenix_Q_conversions_file = 'phoenix_Q_conversions_CFHT-bandpass.txt'
    	phoenix_eniric_Qfactors_file = 'phoenix_eniric_Qfactors_CFHT-bandpass.csv'
    
    #use for Eniric bandpass selected below
    if bandpass == 'Eniric':
    	spriou_fit_qvalues_file = 'spirou_fit_Qvalues_eniric-bandpass.txt'
    	phoenix_Q_conversions_file = 'phoenix_Q_conversions_eniric-bandpass.txt'
    	phoenix_eniric_Qfactors_file = 'phoenix_eniric_Qfactors_eniric-bandpass.csv'
    
    # Constants definitions
    #
    
    ### Physical Constants ###
    h=6.62606957e-34 				# Plank's constant (J.s)
    c=299792458. 					# Light speed (m/s)
    
    ### Telescope parameters ###
    #airmass=1.0 					# Airmass
    d_t=3.57 						# Telescope diameter (m) [https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html]
    
    
    V_sky=0							# negligibleM for IJH bands (Email from E. Artigau)
    #z=1.22603e-09 					# zero point of the photometric system (I band) taken from: 
    z=1.14e-10 						# zero point of the photometric system (H band) taken from: 
    								# http://www.eso.org/observing/etc/doc/formulabook/node12.html
    								# (ergs/s/cmˆ2/A)
    
    ### Instruments Requirements ###
    pixel_size=15e-6 				# pixel size (m)
    N_dark=0 						# Lab tests in Dec/2017 (Email from E. Artigau)
    
    
    
    def calc_fiber_diameter(obs_mode):
    	if obs_mode == 'HA':
    		fiber_diameter=0.4 			# Fiber diameter for HA mode (arcsec)
    	else:
    		fiber_diameter=0.9 			# Fiber diameter for HE mode (arcsec)
    	return fiber_diameter
    
    def calc_resoving_power(obs_mode):		
    	if obs_mode == 'HA':
    		resolving_power=80000 		# Resolving power for HA mode
    	else:
    		resolving_power=70000 		# Resolving power for HE mode
    	return resolving_power
    
    def calc_sampling(obs_mode):
    	if obs_mode == 'HA':
    		sampling=3.8							# Sampling for HA mode (px)
    	else:
    		sampling=4.3 							# Sampling for HE mode (px)
    	return sampling
    
    def calc_N_pix_Y(obs_mode):
    	if obs_mode == 'HA':
    		N_pix_Y=3.6 								# Bin for HA mode
    	else:
    		N_pix_Y=16.7 								# Bin for HE mode
    	return N_pix_Y
    
    def calc_I_max(obs_mode):
    	if obs_mode == 'HA':
    		I_max=3.6*27300 								# Saturation limit for HA mode
    	else:
    		I_max=16.7*27300 								# Saturation limit for HE mode
    	return I_max
    
    
    def calc_flux_spectral_type(st,H):
    	flux_sts=[]
    
    	file_name='inputs/'+st_templates_file
    	lines=reading_table(file_name)
    
    	for line in lines:
    		if st == 'F0V':
    			flux_sts.append(line[1])
    			I=0.348+H
    			Ho=6.998
    		elif st == 'F5V':
    			flux_sts.append(line[2])
    			I=0.519+H
    			Ho=4.738
    		elif st == 'G0V':
    			flux_sts.append(line[3])
    			I=0.699+H
    			Ho=2.905
    		elif st == 'G5V':
    			flux_sts.append(line[4])
    			I=0.821+H
    			Ho=4.614
    		elif st == 'K0V':
    			flux_sts.append(line[5])
    			I=1.009+H
    			Ho=4.803
    		elif st == 'K3V':
    			flux_sts.append(line[6])
    			I=1.275+H
    			Ho=3.469
    		elif st == 'K7V':
    			flux_sts.append(line[7])
    			I=1.665+H
    			Ho=5.499
    		elif st == 'M0V':
    			flux_sts.append(line[8])
    			I=1.739+H
    			Ho=5.843
    		elif st == 'M1V':
    			flux_sts.append(line[9])
    			I=1.799+H
    			Ho=4.393
    		elif st == 'M2V':
    			flux_sts.append(line[10])
    			I=1.833+H
    			Ho=3.640
    		elif st == 'M3V':
    			flux_sts.append(line[11])
    			I=1.928+H
    			Ho=4.843
    		elif st == 'M4V':
    			flux_sts.append(line[12])
    			I=2.137+H
    			Ho=6.627
    		elif st == 'M5V':
    			flux_sts.append(line[13])
    			I=2.354+H
    			Ho=8.014
    		elif st == 'M6V':
    			flux_sts.append(line[14])
    			I=2.848+H
    			Ho=6.482
    		elif st == 'M7V':
    			flux_sts.append(line[15])
    			I=3.147+H
    			Ho=9.201
    		elif st == 'M8V':
    			flux_sts.append(line[16])
    			I=3.668+H
    			Ho=11.066
    		elif st == 'M9V':
    			flux_sts.append(line[17])
    			I=3.996+H
    			Ho=8.905
    		elif st == 'L1V':
    			flux_sts.append(line[18])
    			I=5.477+H
    			Ho=12.041
    		elif st == 'L2V':
    			flux_sts.append(line[19])
    			I=5.704+H
    			Ho=12.392
    		elif st == 'L3V':
    			flux_sts.append(line[20])
    			I=5.931+H
    			Ho=12.380
    		elif st == 'L5V':
    			flux_sts.append(line[21])
    			I=6.385+H
    			Ho=11.895
    		elif st == 'L6V':
    			flux_sts.append(line[22])
    			I=6.612+H
    			Ho=13.099
    		elif st == 'L8V':
    			flux_sts.append(line[23])
    			I=7.066+H
    			Ho=12.204
    		elif st == 'T2V':
    			flux_sts.append(line[24])
    			I=9.108+H
    			Ho=14.090
    		elif st == 'G8V':                 #NG add 06.17.2020
    			flux_sts.append(line[25])     
    			I=0.902+H                     #I-H calculated from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
    			Ho=4.265                      #H mag of star in IRFT spectrum header: http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/Data/G8V_HD75732.txt
    		elif st == 'B3V':                 #NG add 09.06.2022
    			flux_sts.append(line[26])     
    			I=H-0.279                     #I-H calculated from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt (intermediate_preparation/add_stellar_templates/I_Hmag_mamjek.py)
    			Ho=4.840                      #Spectra from HR3454 in CRIRES+ spectrophotometric standard stars: https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar
    		elif st == 'B8V':                 #NG add 09.06.2022
    			flux_sts.append(line[27])     
    			I=H-0.149                     #I-H calculated from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt (intermediate_preparation/add_stellar_templates/I_Hmag_mamjek.py)
    			Ho=3.53                       #Spectra from HR8634 in CRIRES+ spectrophotometric standard stars: https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar
    		elif st == 'B9V':                 #NG add 09.06.2022
    			flux_sts.append(line[28])     
    			I=H-0.076                     #I-H calculated from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt (intermediate_preparation/add_stellar_templates/I_Hmag_mamjek.py)
    			Ho=4.845                      #Spectra from HR4468 in CRIRES+ spectrophotometric standard stars: https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar
    		elif st == 'A1V':                 #NG add 09.06.2022
    			flux_sts.append(line[29])     
    			I=0.016+H                     #I-H calculated from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt (intermediate_preparation/add_stellar_templates/I_Hmag_mamjek.py)
    			Ho=3.71                       #Spectra from HR7950 in CRIRES+ spectrophotometric standard stars: https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar


    	if len(flux_sts) == 0:
    			print("\n################################################")
    			print("# There is no template for this spectral type. #")
    			print("################################################\n")
    			st=calc_spectral_type()
    
    	return (flux_sts,I,Ho)
    
    def calc_total_efficiency(obs_mode,seeing,airmass,I):
    	while seeing < 0.7 or seeing > 1.2:
    		print("\n#################################################")
    		print("# SEEING OUT OF RANGE. Please enter new seeing. #")
    		print("#################################################")
    		seeing=input_seeing()
    
    	while airmass < 1.0 or airmass > 2.0:
    		print("\n###################################################")
    		print("# AIRMASS OUT OF RANGE. Please enter new airmass. #")
    		print("###################################################")
    		airmass=input_airmass()
    
    	if I > 12 and obs_mode == 'HA':
    		print("\n#######################################################")
    		print("# THIS IS A FAINT STAR. Better to observe in HE mode. #")
    		print("#######################################################")
    		
    
    	file_name='inputs/'+wave_range_file
    	lines=reading_table(file_name)
     
    	central_wave=[]
    	beg_wave=[]
    	end_wave=[]
    	order_wave=[]
    
    	for line in lines:
    		central_wave.append(line[1])
    		beg_wave.append(line[3])
    		end_wave.append(line[4])
    		order_wave.append(int(line[0]))
    
    	file_name='inputs/'+tapas_file
    	lines=reading_table(file_name)
    
    	ATM_EFF=[]
    
    	for line in lines:
    
    		if airmass == 1.0:
    			eff=line[1]
    			ATM_EFF.append(eff)
    
    		elif airmass > 1.0 and airmass < 1.02:
    			f=line[1]+(line[2]-line[1])*(airmass-1.00)/0.02
    			eff=f
    			ATM_EFF.append(eff)
    
    		elif airmass == 1.02:
    			eff=line[2]
    			ATM_EFF.append(eff)
    
    		elif airmass > 1.02 and airmass < 1.06:
    			f=line[2]+(line[3]-line[2])*(airmass-1.02)/0.04
    			eff=f
    			ATM_EFF.append(eff)
    
    		elif airmass == 1.06:
    			eff=line[3]
    			ATM_EFF.append(eff)
    
    		elif airmass > 1.06 and airmass < 1.15:
    			f=line[3]+(line[4]-line[3])*(airmass-1.06)/0.09
    			eff=f
    			ATM_EFF.append(eff)
    
    		elif airmass == 1.15:
    			eff=line[4]
    			ATM_EFF.append(eff)
    
    		elif airmass > 1.15 and airmass < 1.31:
    			f=line[4]+(line[5]-line[4])*(airmass-1.15)/0.16
    			eff=f
    			ATM_EFF.append(eff)
    
    		elif airmass == 1.31:
    			eff=line[5]
    			ATM_EFF.append(eff)
    
    		elif airmass > 1.31 and airmass < 1.56:
    			f=line[5]+(line[6]-line[5])*(airmass-1.31)/0.25
    			eff=f
    			ATM_EFF.append(eff)
    
    		elif airmass == 1.56:
    			eff=line[6]
    			ATM_EFF.append(eff)
    
    		elif airmass > 1.56 and airmass < 2.00:
    			f=line[6]+(line[7]-line[6])*(airmass-1.56)/0.44
    			eff=f
    			ATM_EFF.append(eff)
    
    		elif airmass == 2.00:
    			eff=line[7]
    			ATM_EFF.append(eff)
    
    	file_name='inputs/'+effs_file
    	lines=reading_table(file_name)
    
    	EFF=[]
    	wavelengths=[]
    	wavelengths_nm=[]
    	orders=[]
    
    	i=0
    
    	for line in lines:
    		eff=line[1]*line[2]*line[3]*line[4]*line[31]*ATM_EFF[i]
    		wavelengths.append(line[0]*10**(-9))						# Transform wavelength from nm to m
    		wavelengths_nm.append(line[0])
    		orders.append(line[32])
    
    		i=i+1
    
    		if obs_mode == 'HA':
    			eff=eff*line[5]
    
    			if I <= 9:
    
    				if seeing == 0.7:
    					eff=eff*line[6]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[6]+(line[10]-line[6])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[10]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[10]+(line[14]-line[10])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[14]
    					EFF.append(eff)
    
    			elif I > 9 and I <10:
    
    				if seeing == 0.7:
    					f=line[6]+(line[7]-line[6])*(I-9.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f_07=line[6]+(line[7]-line[6])*(I-9.)/1.
    					f_09=line[10]+(line[11]-line[10])*(I-9.)/1.
    					f=f_07+(f_09-f_07)*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					f=line[10]+(line[11]-line[10])*(I-9.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f_09=line[10]+(line[11]-line[10])*(I-9.)/1.
    					f_12=line[14]+(line[15]-line[14])*(I-9.)/1.
    					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					f=line[14]+(line[15]-line[14])*(I-9.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    			elif I == 10:
    
    				if seeing == 0.7:
    					eff=eff*line[7]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[7]+(line[11]-line[7])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[11]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[11]+(line[15]-line[11])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[15]
    					EFF.append(eff)
    
    			elif I > 10 and I < 11:
    
    				if seeing == 0.7:
    					f=line[7]+(line[8]-line[7])*(I-10.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f_07=line[7]+(line[8]-line[7])*(I-10.)/1.
    					f_09=line[11]+(line[12]-line[11])*(I-10.)/1.
    					f=f_07+(f_09-f_07)*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					f=line[11]+(line[12]-line[11])*(I-10.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f_09=line[11]+(line[12]-line[11])*(I-10.)/1.
    					f_12=line[15]+(line[16]-line[15])*(I-10.)/1.
    					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					f=line[15]+(line[16]-line[15])*(I-10.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    			elif I == 11:
    
    				if seeing == 0.7:
    					eff=eff*line[8]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[8]+(line[12]-line[8])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[12]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[12]+(line[16]-line[12])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[16]
    					EFF.append(eff)
    
    			elif I > 11 and I < 12:
    
    				if seeing == 0.7:
    					f=line[8]+(line[9]-line[8])*(I-11.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f_07=line[8]+(line[9]-line[8])*(I-11.)/1.
    					f_09=line[12]+(line[13]-line[12])*(I-11.)/1.
    					f=f_07+(f_09-f_07)*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					f=line[12]+(line[13]-line[12])*(I-11.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f_09=line[12]+(line[13]-line[12])*(I-11.)/1.
    					f_12=line[16]+(line[17]-line[16])*(I-11.)/1.
    					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					f=line[16]+(line[17]-line[16])*(I-11.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    			elif I >= 12:
    
    				if seeing == 0.7:
    					eff=eff*line[9]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[9]+(line[13]-line[9])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[13]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[13]+(line[17]-line[13])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[17]
    					EFF.append(eff)
    
    		else:
    			eff=eff*line[18]
    
    			if I <= 9:
    
    				if seeing == 0.7:
    					eff=eff*line[19]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[19]+(line[23]-line[19])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[23]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[23]+(line[27]-line[23])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[27]
    					EFF.append(eff)
    
    			elif I > 9 and I <10:
    
    				if seeing == 0.7:
    					f=line[19]+(line[20]-line[19])*(I-9.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f_07=line[19]+(line[20]-line[19])*(I-9.)/1.
    					f_09=line[23]+(line[24]-line[23])*(I-9.)/1.
    					f=f_07+(f_09-f_07)*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					f=line[23]+(line[24]-line[23])*(I-9.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f_09=line[23]+(line[24]-line[23])*(I-9.)/1.
    					f_12=line[27]+(line[28]-line[27])*(I-9.)/1.
    					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					f=line[27]+(line[28]-line[27])*(I-9.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    			elif I == 10:
    
    				if seeing == 0.7:
    					eff=eff*line[20]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[20]+(line[24]-line[20])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[24]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[24]+(line[28]-line[24])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[28]
    					EFF.append(eff)
    
    			elif I > 10 and I < 11:
    
    				if seeing == 0.7:
    					f=line[20]+(line[21]-line[20])*(I-10.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f_07=line[20]+(line[21]-line[20])*(I-10.)/1.
    					f_09=line[24]+(line[25]-line[24])*(I-10.)/1.
    					f=f_07+(f_09-f_07)*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					f=line[24]+(line[25]-line[24])*(I-10.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f_09=line[24]+(line[25]-line[24])*(I-10.)/1.
    					f_12=line[28]+(line[29]-line[28])*(I-10.)/1.
    					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					f=line[28]+(line[29]-line[28])*(I-10.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    			elif I == 11:
    
    				if seeing == 0.7:
    					eff=eff*line[21]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[21]+(line[25]-line[21])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[25]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[25]+(line[29]-line[25])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[29]
    					EFF.append(eff)
    
    			elif I > 11 and I < 12:
    
    				if seeing == 0.7:
    					f=line[21]+(line[22]-line[21])*(I-11.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f_07=line[21]+(line[22]-line[21])*(I-11.)/1.
    					f_09=line[25]+(line[26]-line[25])*(I-11.)/1.
    					f=f_07+(f_09-f_07)*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					f=line[25]+(line[26]-line[25])*(I-11.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f_09=line[25]+(line[26]-line[25])*(I-11.)/1.
    					f_12=line[29]+(line[30]-line[29])*(I-11.)/1.
    					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					f=line[29]+(line[30]-line[29])*(I-11.)/1.
    					eff=eff*f
    					EFF.append(eff)
    
    			elif I >= 12:
    
    				if seeing == 0.7:
    					eff=eff*line[22]
    					EFF.append(eff)
    
    				elif seeing > 0.7 and seeing < 0.9:
    					f=line[22]+(line[26]-line[22])*(seeing-0.7)/0.2
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 0.9:
    					eff=eff*line[26]
    					EFF.append(eff)
    
    				elif seeing > 0.9 and seeing < 1.2:
    					f=line[26]+(line[30]-line[26])*(seeing-0.9)/0.3
    					eff=eff*f
    					EFF.append(eff)
    
    				elif seeing == 1.2:
    					eff=eff*line[30]
    					EFF.append(eff)
    
    	total_effs=EFF
    	return (total_effs,wavelengths,wavelengths_nm,orders,central_wave,beg_wave,end_wave,order_wave,obs_mode)
    
    
    # Functions definitions (science calculations)
    #
    
    def calc_readout_noise(t_exp):
        #N_read=82.16*t_exp**(-1) 					# Readout noise in sqrt(electrons/pixel) (Email from E. Artigau)
        N_read=numpy.sqrt(400./(t_exp/5.57 + 1)+36.) # Readout noise in sqrt(electrons/pixel)
        #N_read=0.1          ###### only for caluclation of SNR
        return N_read
    
    def calc_instrument_parameters(d_t,seeing,wavelength,resolving_power):
    #	A=numpy.pi*(d_t/2.)**2 						# Calculate the aperture (assuming no central_obscuration)
    	A=8.8564 									# M1 clear area (mˆ2) [From: https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html]
    	seeing_rads=seeing/3600.*numpy.pi/180. 		# Convert the seeing arcsecs into radians
    	#resolution=wavelength/resolving_power 		# Spectral resolution in meters
    	delta_lambda=wavelength/resolving_power 		# Spectral resolution in meters #NG change 1feb2022 as this is delta lambda not resolution
    	P=h*c/wavelength 							# Energy of one photon in Joules (J/ph)
    	return (A,seeing_rads,delta_lambda,P)
    
    def calc_electrons_rate(z,I,delta_lambda,sampling,A,P,t_exp,total_eff,flux_st,st,H,Ho):
    
    	flux=(10.**(0.4*(Ho-H)))*flux_st                    # Flux of object from IRTF Library template in (W m-2 um-1) scaled by magnitude
    	bin_size_pxl=delta_lambda*1e6/sampling 				# Convert bin size meter to microns
    	bin_size_bin=delta_lambda*1e6 				            # Convert bin size meter to microns
    	N_obj_pxl=flux*bin_size_pxl*t_exp*total_eff*A/P 	# Calculate electrons rate for the object (e-) 
    	N_obj_bin=flux*bin_size_bin*t_exp*total_eff*A/P 	# Calculate electrons rate for the object (e-) 
    	return (N_obj_pxl,N_obj_bin,bin_size_pxl,bin_size_bin,flux)
    
    def calc_signal_to_noise_ratio(N_obj_pxl,N_obj_bin,n_dark,n_bin_y,N_read,t_exp):
    	N_star_pxl=N_obj_pxl
    	N_star_bin=N_obj_bin
    	N_sky=0 #self.N_sky = self.calc_N_obj(self.V_sky, self.z)	Negligeble for IJH bands (not for K band)
    	S_N_pxl = N_star_pxl/numpy.sqrt(N_star_pxl+n_bin_y*N_sky+n_bin_y*N_read**2.+n_bin_y*n_dark*t_exp/3600.)
    	S_N_bin = N_star_bin/numpy.sqrt(N_star_bin+n_bin_y*N_sky+n_bin_y*N_read**2.+n_bin_y*n_dark*t_exp/3600.)
    	return (S_N_pxl,S_N_bin)
    
    def calc_S_N_order(obs_mode,t_exp,d_t,seeing,wavelength,resolving_power,z,I,total_eff,flux_st,N_dark,st,H,Ho):
    	fiber_diameter=calc_fiber_diameter(obs_mode)
    	sampling=calc_sampling(obs_mode)
    	N_read=calc_readout_noise(t_exp)
    	
    	### Calculating Instrument Parameters ###
    
    	(A,seeing_rads,delta_lambda,P)=calc_instrument_parameters(d_t,seeing,wavelength,resolving_power) #change resolution to delta_lambda 1feb2022
    
    	### Calculating rate of electrons per second per bin from the object ###
    	# Formulae from ESO exposure calculator
    	#
    
    	'''
    
    	INPUT
    
    	I = Object magnitude
    	z = zero point of the photometric system - taken from http://www.eso.org/observing/etc/doc/formulabook/node12.html
    
    	OUTPUT
    
    	N_obj = Object eletrons rate (ergs/s/cmˆ2/A)
    
    	'''
    
    	(N_obj_pxl,N_obj_bin,bin_size_pxl,bin_size_bin,flux)=calc_electrons_rate(
    		z,I,delta_lambda,sampling,A,P,t_exp,total_eff,flux_st,st,H,Ho) #change resolution to delta_lambda 1feb2022
    
    	### Calculating te Signal-to-Noise Rate ###
    	# Formulae from ESO exposure calculator
    	#
    
    	n_bin_y=calc_N_pix_Y(obs_mode)      					# Bins in y direction
    	
    	(S_N_pxl,S_N_bin)=calc_signal_to_noise_ratio(
    		N_obj_pxl,N_obj_bin,N_dark,n_bin_y,N_read,t_exp)
    
    	return (N_obj_pxl,N_obj_bin,bin_size_pxl,bin_size_bin,flux,S_N_pxl,S_N_bin)
    
    
    resolving_power=calc_resoving_power(obs_mode)
    I_max=calc_I_max(obs_mode)
    
    (flux_sts,I,Ho)=calc_flux_spectral_type(st,H)
    
    (total_effs,wavelengths,wavelengths_nm,orders,central_wave,beg_wave,end_wave,order_wave,obs_mode
    ) = calc_total_efficiency(obs_mode,seeing,airmass,I)
    
    
    S_N_pxl=[]
    S_N_bin=[]
    FLUX=[]
    N_OBJ=[]
    BIN_SIZE=[]
    S_N_pxl_mean=[]
    S_N_bin_mean=[]
    N_OBJ_mean=[]
    BIN_SIZE_mean=[]
    FLUX_mean=[]
    EFF_mean=[]
    SATURATION=[]
    
    
    i=0
    SNR_H_pxl=0
    SNR_H_bin=0
    order_now=orders[0]
    s_n_pxl_mean=0
    s_n_bin_mean=0
    n_obj_mean=0
    bin_size_mean=0
    eff_mean=0
    flux_mean=0
    divisor=0
    cont_saturation=0.0
    
    while i<len(wavelengths):
    	wavelength=wavelengths[i]
    	wavelength_nm=wavelengths_nm[i]
    	order=orders[i]
    	total_eff=total_effs[i]
    	flux_st=flux_sts[i]
    	if i == len(wavelengths)-1:
    		delta_wavelength = (wavelengths_nm[-1]-wavelengths_nm[-2])*0.001 #nm to microns
    	else:
    		delta_wavelength = (wavelengths_nm[i+1]-wavelengths_nm[i])*0.001
    
    	if order_now == order:
    		(n_obj_pxl,n_obj_bin,bin_size_pxl,bin_size_bin,flux,s_n_pxl,s_n_bin)=calc_S_N_order(obs_mode,t_exp,d_t,seeing,wavelength,resolving_power,z,I,total_eff,flux_st,N_dark,st,H,Ho)
    		S_N_pxl.append(s_n_pxl)
    		S_N_bin.append(s_n_bin)
    		FLUX.append(flux)
    		#N_OBJ.append(n_obj_pxl)
    		#For total number of electrons: multiply by ((delta wavlength of the array) / (bin size))
    		#to account for the array size of the effs/tapas/stellar template wavelength values
    		N_OBJ.append(n_obj_pxl*(delta_wavelength/bin_size_pxl)) 
    
    		if n_obj_pxl > I_max:
    			cont_saturation+=1.0
    
    		BIN_SIZE.append(bin_size_pxl)
    		bin_size_mean=bin_size_mean+bin_size_pxl
    		n_obj_mean=n_obj_mean+n_obj_pxl
    		s_n_pxl_mean=s_n_pxl_mean+s_n_pxl
    		s_n_bin_mean=s_n_bin_mean+s_n_bin
    		eff_mean=eff_mean+total_eff
    		flux_mean=flux_mean+flux
    		order_now=order
    		divisor=divisor+1
    	else:
    		sn_pxl=s_n_pxl_mean/divisor
    		sn_bin=s_n_bin_mean/divisor
    		nobj=n_obj_mean/divisor
    		binsize=bin_size_mean/divisor
    		effall=eff_mean/divisor
    		fluxobj=flux_mean/divisor
    		saturation=cont_saturation/divisor
    		S_N_pxl_mean.append(sn_pxl)
    		S_N_bin_mean.append(sn_bin)
    		N_OBJ_mean.append(nobj)
    		BIN_SIZE_mean.append(binsize)
    		EFF_mean.append(effall)
    		FLUX_mean.append(fluxobj)
    		SATURATION.append(saturation)
    		s_n_pxl_mean=0
    		s_n_bin_mean=0
    		n_obj_mean=0
    		bin_size_mean=0
    		eff_mean=0
    		flux_mean=0
    		divisor=0
    		cont_saturation=0.0
    		(n_obj_pxl,n_obj_bin,bin_size_pxl,bin_size_bin,flux,s_n_pxl,s_n_bin)=calc_S_N_order(obs_mode,t_exp,d_t,seeing,wavelength,resolving_power,z,I,total_eff,flux_st,N_dark,st,H,Ho)
    		S_N_pxl.append(s_n_pxl)
    		S_N_bin.append(s_n_bin)
    		FLUX.append(flux)
    		N_OBJ.append(n_obj_pxl)
    
    		if n_obj_pxl > I_max:
    			cont_saturation+=1.0
    
    		BIN_SIZE.append(bin_size_pxl)
    		s_n_pxl_mean=s_n_pxl_mean+s_n_pxl
    		s_n_bin_mean=s_n_bin_mean+s_n_bin
    		n_obj_mean=n_obj_mean+n_obj_pxl
    		bin_size_mean=bin_size_mean+bin_size_pxl
    		eff_mean=eff_mean+total_eff
    		flux_mean=flux_mean+flux
    		order_now=order
    		divisor=divisor+1
    
    	i=i+1
    
    #Hband SNR
    hband_wave = 1619
    hband_ind = np.argmin(np.abs(np.array(wavelengths_nm) - hband_wave))
    SNR_pxl_H = S_N_pxl[hband_ind]
    SNR_bin_H = S_N_bin[hband_ind]
    sn_h = SNR_bin_H
    
    
    sn_pxl=s_n_pxl_mean/divisor
    sn_bin=s_n_bin_mean/divisor
    nobj=n_obj_mean/divisor
    binsize=bin_size_mean/divisor
    effall=eff_mean/divisor
    fluxobj=flux_mean/divisor
    
    saturation=cont_saturation/divisor
    S_N_pxl_mean.append(sn_pxl)
    S_N_bin_mean.append(sn_bin)
    N_OBJ_mean.append(nobj)
    BIN_SIZE_mean.append(binsize)
    EFF_mean.append(effall)
    FLUX_mean.append(fluxobj)
    SATURATION.append(saturation)
    
    print ("=================================================================")
    
    
    if bandpass == 'Eniric':
        print('### Eniric bandpasses ###')
        ##Y-band
        wlnmin_y = 1000.0
        wlnmax_y = 1100.0
        #J-band
        wlnmin_j = 1170.0
        wlnmax_j = 1330.0
        #H-band
        wlnmin_h = 1500.0
        wlnmax_h = 1750.0
        #K-band
        #wlnmin_k = 2070.0
        #wlnmax_k = 2350.0
    
    if bandpass == 'CFHT':
        print('#### CFHT/WIRCAM bandpasses ###')
        ##Y-band
        wlnmin_y = 938.6
        wlnmax_y = 1113.4
        #J-band
        wlnmin_j = 1153.5
        wlnmax_j = 1354.4
        #H-band
        wlnmin_h = 1462.8
        wlnmax_h = 1808.5
        ##K-band
        ##wlnmin_k = 1957.7
        ##wlnmax_k = 2343.1
    
    ### OLD NIRPS ETC LIMITS ###
    #wlnmin_y = 980
    #wlnmax_y = 1110
    #wlnmin_j = 1200
    #wlnmax_j = 1330
    #wlnmin_h = 1510
    #wlnmax_h = 1735
    
    print('Y band: '+str(wlnmin_y)+' - '+str(wlnmax_y)+' nm')
    print('J band: '+str(wlnmin_j)+' - '+str(wlnmax_j)+' nm')
    print('H band: '+str(wlnmin_h)+' - '+str(wlnmax_h)+' nm')
    
    print ("\nTEMPLATE STAR: %s"%(st))
    print ("H (mag): %.3f"%(H))
    print ("Exposure time (seconds): %.1f\n"%(t_exp))
    
    print ("Observation Mode: %s"%(obs_mode))
    print ("seeing (arcsec): %.1f"%(seeing))
    print ("airmass: %.1f\n"%(airmass))
    
    print ("Saturation limit (e-/pxl): %d\n\n"%(I_max))
    
    print ("=================================================================")
    
    #makedirs("outputs", exist_ok=True)
    #text_file = open("outputs/order_snrs.txt", "w")
    #text_file.write("order central_wave (beg-end)    Eff.    object     snr      snr         sat \n")
    #text_file.write("          (um)                         (e-/pxl)  (ph/pxl) (ph/res elem) (%) \n") 
    
    i=0
    
    SN_pxl=0
    SN_pxl_Y=0
    SN_pxl_J=0
    SN_pxl_H=0
    mean_pxl_Y=0
    mean_pxl_J=0
    mean_pxl_H=0
    SN_bin=0
    SN_bin_Y=0
    SN_bin_J=0
    SN_bin_H=0
    mean_bin_Y=0
    mean_bin_J=0
    mean_bin_H=0
    
    #Efficiency averages
    mean_all_eff = 0
    mean_eff_Y = 0
    mean_eff_J = 0
    mean_eff_H = 0
    
    while i < len(order_wave):
    
    	#text_file.write('# %s %7.2f (%7.2f-%7.2f) '%(repr(order_wave[i]).rjust(3),central_wave[i],beg_wave[i],end_wave[i]))
    	#text_file.write('%5.3f %.5e %5.1f %8.1f %s \n'%(EFF_mean[i],Decimal(N_OBJ_mean[i]),S_N_pxl_mean[i],S_N_bin_mean[i],repr(int(SATURATION[i]*100.)).rjust(10)))
    	SN_pxl=SN_pxl+S_N_pxl_mean[i]
    	SN_bin=SN_bin+S_N_bin_mean[i]
    	mean_all_eff = mean_all_eff+EFF_mean[i]
    	if wlnmin_y <= central_wave[i] <= wlnmax_y:
    		SN_pxl_Y=SN_pxl_Y+S_N_pxl_mean[i]
    		mean_pxl_Y=mean_pxl_Y+1
    		SN_bin_Y=SN_bin_Y+S_N_bin_mean[i]
    		mean_bin_Y=mean_bin_Y+1
    		mean_eff_Y=mean_eff_Y+EFF_mean[i]
    	if wlnmin_j <= central_wave[i] <= wlnmax_j:
    		SN_pxl_J=SN_pxl_J+S_N_pxl_mean[i]
    		mean_pxl_J=mean_pxl_J+1
    		SN_bin_J=SN_bin_J+S_N_bin_mean[i]
    		mean_bin_J=mean_bin_J+1
    		mean_eff_J=mean_eff_J+EFF_mean[i]
    	if wlnmin_h <= central_wave[i] <= wlnmax_h:
    		SN_pxl_H=SN_pxl_H+S_N_pxl_mean[i]
    		mean_pxl_H=mean_pxl_H+1
    		SN_bin_H=SN_bin_H+S_N_bin_mean[i]
    		mean_bin_H=mean_bin_H+1
    		mean_eff_H=mean_eff_H+EFF_mean[i]
    	i=i+1
    
    #text_file.close()
    print ("-----------------------------------------------------------------")
    print("\n SIGNAL TO NOISE RATIO:\n")
    print ("Mean S/N: %5.1f (ph/pxl) | %5.1f (ph/res elem)"%(SN_pxl/len(order_wave),SN_bin/len(order_wave)))
    print ("Mean S/N (ph/pxl): Y=%5.1f | J=%5.1f | H=%5.1f"%(SN_pxl_Y/mean_pxl_Y,SN_pxl_J/mean_pxl_J,SN_pxl_H/mean_pxl_H))
    print ("         (ph/res elem): Y=%5.1f | J=%5.1f | H=%5.1f\n\n"%(SN_bin_Y/mean_bin_Y,SN_bin_J/mean_bin_J,SN_bin_H/mean_bin_H))
    print ("S/N in H (1625 nm): %5.1f (ph/pxl) | %5.1f (ph/res elem)\n"%(SNR_pxl_H,SNR_bin_H))
    print ("-----------------------------------------------------------------")
    print("\nMean Efficiency: %5.3f "%(mean_all_eff/len(EFF_mean)))
    print("Mean Efficiencies Y=%5.3f | J=%5.3f | H=%5.3f \n"%(mean_eff_Y/mean_pxl_Y,mean_eff_J/mean_pxl_J,mean_eff_H/mean_pxl_H))
    print ("-----------------------------------------------------------------")
    
    N_OBJ_arr = np.array(N_OBJ)
    wavelengths_nm = np.array(wavelengths_nm)
    
    
    def calc_rv_precision(st,obs_mode,sn_h):    
        #temps = np.array([5000,4500,4000,3900,3700,3600,3400,3200,3000,2800,2700,2600,2500]) #original Teff values available from Spirou Template interpolation
        #stypes = np.array(['K2V','K5V','K7V','M0V','M1V','M2V','M3V','M4V','M5V','M6V','M7V','M8V','M9V'])
        temps = np.array([5000,4000,3900,3700,3600,3400,3200,3000,2800,2700,2600,2500])
        stypes = np.array(['K3V','K7V','M0V','M1V','M2V','M3V','M4V','M5V','M6V','M7V','M8V','M9V'])
        spirou_qs = spirou_results = pd.read_csv('inputs/'+spriou_fit_qvalues_file, sep=',', header=0)
        q_resolution_conversion = pd.read_csv('inputs/'+phoenix_Q_conversions_file, sep=',', header=0)
        eniric_pheonix_results = pd.read_csv('inputs/'+phoenix_eniric_Qfactors_file, sep=',', header=0)
    
        index = [np.where(stypes == st)[0]][0]
        temp = temps[index][0]
        print('')
        print(' RADIAL VELOCITY PRECISION:')
        print('')
        print(st,'Teff =',temp)
        print('YJH total e-: %i %i %i'%(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))]),
                                  sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_j))]),
                                  sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_h))])))
        print('')
        print('-- RV precision from Spirou Template Qs --')
        print('')
        if obs_mode == 'HE':
            resolution = '70k'
            sampling = 4.3
            spirou_quality_y = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_Y'].values[0]
            spirou_quality_j = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_J'].values[0]
            spirou_quality_h = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_H'].values[0]
        if obs_mode == 'HA': 
            resolution = '80k'
            sampling = 3.8
            spirou_quality_y = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_Y'].values[0]/q_resolution_conversion[(spirou_qs['Teff'] == temp)]['Q_Y'].values[0]
            spirou_quality_j = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_J'].values[0]/q_resolution_conversion[(spirou_qs['Teff'] == temp)]['Q_J'].values[0]
            spirou_quality_h = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_H'].values[0]/q_resolution_conversion[(spirou_qs['Teff'] == temp)]['Q_H'].values[0]
    
        rvprec_y = c/(spirou_quality_y*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))])))
        rvprec_j = c/(spirou_quality_j*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_j) & (wavelengths_nm <= wlnmax_j))])))
        rvprec_h = c/(spirou_quality_h*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_h) & (wavelengths_nm <= wlnmax_h))])))
        combine_yjh_sp = 1/np.sqrt(rvprec_y**-2+rvprec_j**-2+rvprec_h**-2)
        print('YJH Quality factors: %.1f %.1f %.1f'%(spirou_quality_y,spirou_quality_j,spirou_quality_h))
        print('YJH RV precision (m/s): %.2f %.2f %.2f'%(rvprec_y,rvprec_j,rvprec_h))
        print('YJH combined RV precision (m/s): %.2f'%combine_yjh_sp)
           
        def eniric_rv_precision(vsini):
            print('')
            print('vsini = %.1f km/s'%vsini)
            match_ph = ((eniric_pheonix_results['temp'] == temp) &
                            (eniric_pheonix_results['sampling'] == sampling) & 
                            (eniric_pheonix_results['vsini'] == vsini) &
                            (eniric_pheonix_results['resolution'] == resolution))
    
            if bandpass == 'CFHT':
                ph_quality_y = eniric_pheonix_results[((match_ph) & (eniric_pheonix_results['band'] == 'Y_CFHT'))]['quality'].values[0]
                ph_quality_j = eniric_pheonix_results[((match_ph) & (eniric_pheonix_results['band'] == 'J_CFHT'))]['quality'].values[0]
                ph_quality_h = eniric_pheonix_results[((match_ph) & (eniric_pheonix_results['band'] == 'H_CFHT'))]['quality'].values[0]
            elif bandpass == 'Eniric':
                ph_quality_y = eniric_pheonix_results[((match_ph) & (eniric_pheonix_results['band'] == 'Y'))]['quality'].values[0]
                ph_quality_j = eniric_pheonix_results[((match_ph) & (eniric_pheonix_results['band'] == 'J'))]['quality'].values[0]
                ph_quality_h = eniric_pheonix_results[((match_ph) & (eniric_pheonix_results['band'] == 'H'))]['quality'].values[0]
            else:
                ph_quality_y = 'NaN'
                ph_quality_j = 'NaN'
                ph_quality_h = 'NaN'
            print('YJH Quality factors:',ph_quality_y,ph_quality_j,ph_quality_h)
            
            rvprec_y = c/(ph_quality_y*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))])))
            rvprec_j = c/(ph_quality_j*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_j) & (wavelengths_nm <= wlnmax_j))])))
            rvprec_h = c/(ph_quality_h*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_h) & (wavelengths_nm <= wlnmax_h))])))
            combine_yjh_en = 1/np.sqrt(rvprec_y**-2+rvprec_j**-2+rvprec_h**-2)
            
            print('YJH RV precision (m/s): %.2f %.2f %.2f'%(rvprec_y,rvprec_j,rvprec_h))
            print('YJH combined RV precision (m/s): %.2f'%combine_yjh_en)
            return(combine_yjh_en)
        
        print('')
        print('-- RV precision from Eniric Pheonix spectra Qs --')
        vsinis = [0.1,1.0,5.0,10.0]
        eniric_rv_vsinis = np.zeros(len(vsinis))
        for ivsini in range(len(vsinis)):
            eniric_rv_vsinis[ivsini] = eniric_rv_precision(vsinis[ivsini])
        return(combine_yjh_sp,eniric_rv_vsinis)
    
    stypes_rv = np.array(['K3V','K7V','M0V','M1V','M2V','M3V','M4V','M5V','M6V','M7V','M8V','M9V'])
    if len(np.where(stypes_rv == st)[0]) == 1:
        rv_total_spirou,rv_total_eniric = calc_rv_precision(st,obs_mode,sn_h)
        spRV = '{:.2f}'.format(rv_total_spirou)
        enRV_vsini01 = '{:.2f}'.format(float(rv_total_eniric[0]))
        enRV_vsini1 = '{:.2f}'.format(float(rv_total_eniric[1]))
        enRV_vsini5 = '{:.2f}'.format(float(rv_total_eniric[2]))
        enRV_vsini10 = '{:.2f}'.format(float(rv_total_eniric[3]))
    else:
        print("# It's not possible to calculate the RV precision of this spectral type with the NIRPS ETC currently #\n")
        spRV = '--'
        enRV_vsini01 = '--'
        enRV_vsini1 = '--'
        enRV_vsini5 = '--'
        enRV_vsini10 = '--'
        
    
    print ("=================================================================\n\n")
    output_targets.write(target+' %.2f %.2f %.2f %.2f %.2f '%(SN_pxl/len(order_wave),SN_pxl_Y/mean_pxl_Y,SN_pxl_J/mean_pxl_J,SN_pxl_H/mean_pxl_H,SNR_pxl_H)+spRV+' '+enRV_vsini01+' '+enRV_vsini1+' '+enRV_vsini5+' '+enRV_vsini10+' \n')

output_targets.close()
