#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (C) 2018, Bruno L. Canto Martins

# NIRPS ETC code python
# Last modification - March, 2018.
# 

# Imports
#

import numpy
from os import system, chdir
import matplotlib.pyplot as plt
import matplotlib as mpl
#from astropy import constants
from math import *
from nirps_etc_lib import *
from decimal import Decimal


# Constants definitions
#

### Physical Constants ###

h=6.62606957e-34 				# Plank's constant (J.s)
c=299792458. 					# Light speed (m/s)

### Telescope parameters ###

#airmass=1.0 					# Airmass
d_t=3.57 						# Telescope diameter (m) [https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html]


V_sky=0							# negligible for IJH bands (Email from E. Artigau)
#z=1.22603e-09 					# zero point of the photometric system (I band) taken from: 
z=1.14e-10 						# zero point of the photometric system (H band) taken from: 
								# http://www.eso.org/observing/etc/doc/formulabook/node12.html
								# (ergs/s/cmˆ2/A)

### Instruments Requirements ###

pixel_size=15e-6 				# pixel size (m)
N_dark=0 						# Lab tests in Dec/2017 (Email from E. Artigau)

# Functions definitions (user interactions)
#

def input_observation_mode():
	obs_mode=input("Observation Mode (HA/HE): ") # Observation mode
	obs_mode=obs_mode.rstrip()
	return obs_mode 

def input_seeing():
	seeing=float(input("Seeing, in arcsec (range 0.7-2.0): ")) 	# Seeing of the telescope (FWHM) in arcseconds
	return seeing 													# Range from 0.7 to 2.0 (2019_01_21) --> From 0.7 - 1.2 to 0.7 - 2.0

def input_airmass():
	airmass=float(input("Airmass (range 1.00-5.50): ")) 	# Airmass
	return airmass 													

def input_magnitude():
	H=float(input("Object magnitude (H band): ")) # Object Magnitude
	return H

def input_exposure_time():
	t_exp=float(input("Exposure time (in sec): ")) # Exposure time in seconds
	return t_exp

def input_spectral_type():
	st=input("Spectral type (F0V/F5V/G0V/G5V/K0V/K3V/K7V/M0V/M1V/M2V/M3V/M4V/M5V/M6V/M7V/M8V/M9V/L1V/L2V/L3V/L4V/L5V/L6V/L8V/T2V): ") # ST Models for IRFT
	st=st.rstrip()

	return (st) 

# Functions calculations
#

def calc_fiber_diameter(obs_mode):
	if obs_mode == 'HA':
		fiber_diameter=0.4 			# Fiber diameter for HA mode (arcsec)
	else:
		fiber_diameter=0.9 			# Fiber diameter for HE mode (arcsec)
	return fiber_diameter

def calc_resoving_power(obs_mode):		
	if obs_mode == 'HA':
		resolving_power=100000 		# Resolving power for HA mode
	else:
		resolving_power=75000 		# Resolving power for HE mode
	return resolving_power

def calc_sampling(obs_mode):
	if obs_mode == 'HA':
		sampling=3.1 							# Sampling for HA mode (px)
	else:
		sampling=4.1 							# Sampling for HE mode (px)
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

def calc_delta_rv(st,obs_mode,sn_h):
	file_name='sigma_RV.txt'
	TABLE_ARCHIVE=open(file_name,'rt')
	lines=TABLE_ARCHIVE.readlines()

	ok=0
	y1_cond2=9999
	y1_cond3=9999
	j1_cond2=9999
	j1_cond3=9999
	h1_cond2=9999
	h1_cond3=9999
	y5_cond2=9999
	y5_cond3=9999
	j5_cond2=9999
	j5_cond3=9999
	h5_cond2=9999
	h5_cond3=9999
	y10_cond2=9999
	y10_cond3=9999
	j10_cond2=9999
	j10_cond3=9999
	h10_cond2=9999
	h10_cond3=9999

	for line in lines:
		colunas=line.split()
		if colunas[0]==st and colunas[1]==obs_mode:
			y1_cond2=float(colunas[2])
			y1_cond3=float(colunas[3])
			j1_cond2=float(colunas[4])
			j1_cond3=float(colunas[5])
			h1_cond2=float(colunas[6])
			h1_cond3=float(colunas[7])
			y5_cond2=float(colunas[8])
			y5_cond3=float(colunas[9])
			j5_cond2=float(colunas[10])
			j5_cond3=float(colunas[11])
			h5_cond2=float(colunas[12])
			h5_cond3=float(colunas[13])
			y10_cond2=float(colunas[14])
			y10_cond3=float(colunas[15])
			j10_cond2=float(colunas[16])
			j10_cond3=float(colunas[17])
			h10_cond2=float(colunas[18])
			h10_cond3=float(colunas[19])
			ok=1

	if y1_cond2==9999:
		print("\n###################################################")
		print("# It's not possible to calculate de RV precision! #")
		print("###################################################\n")
	else:
		deltaRV1_y2=y1_cond2*100/sn_h
		deltaRV1_y3=y1_cond3*100/sn_h
		deltaRV1_j2=j1_cond2*100/sn_h
		deltaRV1_j3=j1_cond3*100/sn_h
		deltaRV1_h2=h1_cond2*100/sn_h
		deltaRV1_h3=h1_cond3*100/sn_h
		delta_RV1_2=1/sqrt((1/deltaRV1_y2)**2+(1/deltaRV1_j2)**2+(1/deltaRV1_h2)**2)
		delta_RV1_3=1/sqrt((1/deltaRV1_y3)**2+(1/deltaRV1_j3)**2+(1/deltaRV1_h3)**2)
		deltaRV5_y2=y5_cond2*100/sn_h
		deltaRV5_y3=y5_cond3*100/sn_h
		deltaRV5_j2=j5_cond2*100/sn_h
		deltaRV5_j3=j5_cond3*100/sn_h
		deltaRV5_h2=h5_cond2*100/sn_h
		deltaRV5_h3=h5_cond3*100/sn_h
		delta_RV5_2=1/sqrt((1/deltaRV5_y2)**2+(1/deltaRV5_j2)**2+(1/deltaRV5_h2)**2)
		delta_RV5_3=1/sqrt((1/deltaRV5_y3)**2+(1/deltaRV5_j3)**2+(1/deltaRV5_h3)**2)
		deltaRV10_y2=y10_cond2*100/sn_h
		deltaRV10_y3=y10_cond3*100/sn_h
		deltaRV10_j2=j10_cond2*100/sn_h
		deltaRV10_j3=j10_cond3*100/sn_h
		deltaRV10_h2=h10_cond2*100/sn_h
		deltaRV10_h3=h10_cond3*100/sn_h
		delta_RV10_2=1/sqrt((1/deltaRV10_y2)**2+(1/deltaRV10_j2)**2+(1/deltaRV10_h2)**2)
		delta_RV10_3=1/sqrt((1/deltaRV10_y3)**2+(1/deltaRV10_j3)**2+(1/deltaRV10_h3)**2)


		print("\n RADIAL VELOCITY PRECISION: \n")

		print("vsini = 1.0 km/s:")
		print("RV precision (m/s) for YJH bands (min): %5.1f %5.1f %5.1f"%(deltaRV1_y3,deltaRV1_j3,deltaRV1_h3))
		print("RV precision (m/s) for YJH bands (max): %5.1f %5.1f %5.1f"%(deltaRV1_y2,deltaRV1_j2,deltaRV1_h2))
		print("Total RV precision (m/s) (min): %5.1f"%delta_RV1_3)
		print("Total RV precision (m/s) (max): %5.1f"%delta_RV1_2)
		
		print("\nvsini = 5.0 km/s:")
		print("RV precision (m/s) for YJH bands (min): %5.1f %5.1f %5.1f"%(deltaRV5_y3,deltaRV5_j3,deltaRV5_h3))
		print("RV precision (m/s) for YJH bands (max): %5.1f %5.1f %5.1f"%(deltaRV5_y2,deltaRV5_j2,deltaRV5_h2))
		print("Total RV precision (m/s) (min): %5.1f"%delta_RV5_3)
		print("Total RV precision (m/s) (max): %5.1f"%delta_RV5_2)
		
		print("\nvsini = 10.0 km/s:")
		print("RV precision (m/s) for YJH bands (min): %5.1f %5.1f %5.1f"%(deltaRV10_y3,deltaRV10_j3,deltaRV10_h3))
		print("RV precision (m/s) for YJH bands (max): %5.1f %5.1f %5.1f"%(deltaRV10_y2,deltaRV10_j2,deltaRV10_h2))
		print("Total RV precision (m/s) (min): %5.1f"%delta_RV10_3)
		print("Total RV precision (m/s) (max): %5.1f"%delta_RV10_2)

def calc_flux_spectral_type(st,H):
	flux_sts=[]

	file_name='ST_templates.txt'
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

	if len(flux_sts) == 0:
			print("\n################################################")
			print("# There is no template for this spectral type. #")
			print("################################################\n")
			st=calc_spectral_type()

	return (flux_sts,I,Ho)

def calc_total_efficiency(obs_mode,seeing,airmass,I):
	while seeing < 0.7 or seeing > 2.0:
		print("\n#################################################")
		print("# SEEING OUT OF RANGE. Please enter new seeing. #")
		print("#################################################")
		seeing=input_seeing()

	while airmass < 1.00 or airmass > 5.55:
		print("\n###################################################")
		print("# AIRMASS OUT OF RANGE. Please enter new airmass. #")
		print("###################################################")
		airmass=input_airmass()

	if I > 12 and obs_mode == 'HA':
		print("\n#######################################################")
		print("# THIS IS A FAINT STAR. Better to observe in HE mode. #")
		print("#######################################################")
		
#		obs_mode='HE'

	file_name='wave_range.txt'
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

	file_name='tapas.txt'
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

		elif airmass > 2.00 and airmass < 2.90:
			f=line[7]+(line[8]-line[7])*(airmass-2.00)/0.90
			eff=f
			ATM_EFF.append(eff)

		elif airmass == 2.90:
			eff=line[8]
			ATM_EFF.append(eff)

		elif airmass > 2.90 and airmass < 5.55:
			f=line[8]+(line[9]-line[8])*(airmass-2.90)/2.65
			eff=f
			ATM_EFF.append(eff)

		elif airmass == 5.55:
			eff=line[9]
			ATM_EFF.append(eff)

	file_name='effs.txt'
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[10]+(line[14]-line[10])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f_09=line[10]+(line[11]-line[10])*(I-9.)/1.
					f_12=line[14]+(line[15]-line[14])*(I-9.)/1.
					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[11]+(line[15]-line[11])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f_09=line[11]+(line[12]-line[11])*(I-10.)/1.
					f_12=line[15]+(line[16]-line[15])*(I-10.)/1.
					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[12]+(line[16]-line[12])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f_09=line[12]+(line[13]-line[12])*(I-11.)/1.
					f_12=line[16]+(line[17]-line[16])*(I-11.)/1.
					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[13]+(line[17]-line[13])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[23]+(line[27]-line[23])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f_09=line[23]+(line[24]-line[23])*(I-9.)/1.
					f_12=line[27]+(line[28]-line[27])*(I-9.)/1.
					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[24]+(line[28]-line[24])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f_09=line[24]+(line[25]-line[24])*(I-10.)/1.
					f_12=line[28]+(line[29]-line[28])*(I-10.)/1.
					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[25]+(line[29]-line[25])*(seeing-0.9)/0.3
					eff=eff*f
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

				elif seeing > 1.2 and seeing <= 2.0:
					f_09=line[25]+(line[26]-line[25])*(I-11.)/1.
					f_12=line[29]+(line[30]-line[29])*(I-11.)/1.
					f=f_09+(f_12-f_09)*(seeing-0.9)/0.3
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

				elif seeing > 1.2 and seeing <= 2.0:
					f=line[26]+(line[30]-line[26])*(seeing-0.9)/0.3
					eff=eff*f
					EFF.append(eff)

	total_effs=EFF

	return (total_effs,wavelengths,wavelengths_nm,orders,central_wave,beg_wave,end_wave,order_wave,obs_mode)


# Functions definitions (science calculations)
#

def calc_readout_noise(t_exp):
	N_read=82.16*t_exp**(-1/2) 					# Readout noise in sqrt(electrons/pixel) (Email from E. Artigau)
	return N_read

def calc_instrument_parameters(d_t,seeing,wavelength,resolving_power):
#	A=numpy.pi*(d_t/2.)**2 						# Calculate the aperture (assuming no central_obscuration)
	A=8.8564 									# M1 clear area (mˆ2) [From: https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html]
	seeing_rads=seeing/3600.*numpy.pi/180. 		# Convert the seeing arcsecs into radians
	resolution=wavelength/resolving_power 		# Spectral resolution in meters
	P=h*c/wavelength 							# Energy of one photon in Joules (J/ph)
	return (A,seeing_rads,resolution,P)

def calc_electrons_rate(z,I,resolution,sampling,A,P,t_exp,total_eff,flux_st,st,H,Ho):
	
	flux=(10.**(-1-.4*(H-Ho)))*flux_st 					# Flux of the object (ergs/s/cmˆ2/A) for IRTF Library
	bin_size_pxl=resolution*1e10/sampling 				# Convert bin size meter to Angstroms
	bin_size_bin=resolution*1e10 					    # Convert bin size meter to Angstroms
	A_cm=A*1e4 											# Convert area m to cmˆ2
	P_ergs=P*1e7										# Convert energy Joules to ergs
	N_obj_pxl=flux*bin_size_pxl*t_exp*total_eff*A_cm/P_ergs 	# Calculate electrons rate for the object (e-) 
	N_obj_bin=flux*bin_size_bin*t_exp*total_eff*A_cm/P_ergs 	# Calculate electrons rate for the object (e-) 
	return (N_obj_pxl,N_obj_bin,bin_size_pxl,bin_size_bin,flux)

def calc_signal_to_noise_ratio(N_obj_pxl,N_obj_bin,n_dark,n_bin_y,N_read,t_exp):
	N_star_pxl=N_obj_pxl
	N_star_bin=N_obj_bin
	N_sky=0 #self.N_sky = self.calc_N_obj(self.V_sky, self.z)	Negligeble for IJH bands (not for K band)
	#n_dark=N_dark    									# Dark counts per hour
	#n_bin_y=N_pix_Y(obs_mode)      					# Bins in y direction
	#n_ron=N_read 		    							# Read out noise
	S_N_pxl = N_star_pxl/numpy.sqrt(N_star_pxl+n_bin_y*N_sky+n_bin_y*N_read**2.+n_bin_y*n_dark*t_exp/3600.)
	S_N_bin = N_star_bin/numpy.sqrt(N_star_bin+n_bin_y*N_sky+n_bin_y*N_read**2.+n_bin_y*n_dark*t_exp/3600.)
	return (S_N_pxl,S_N_bin)

def calc_S_N_order(obs_mode,t_exp,d_t,seeing,wavelength,resolving_power,z,I,total_eff,flux_st,N_dark,st,H,Ho):
	fiber_diameter=calc_fiber_diameter(obs_mode)
	sampling=calc_sampling(obs_mode)
	N_read=calc_readout_noise(t_exp)
	
	### Calculating Instrument Parameters ###

	(A,seeing_rads,resolution,P)=calc_instrument_parameters(d_t,seeing,wavelength,resolving_power)

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

	(N_obj_pxl,N_obj_bin,bin_size_pxl,bin_size_bin,flux)=calc_electrons_rate(z,I,resolution,sampling,A,P,t_exp,total_eff,flux_st,st,H,Ho)

	### Calculating te Signal-to-Noise Rate ###
	# Formulae from ESO exposure calculator
	#

	n_bin_y=calc_N_pix_Y(obs_mode)      					# Bins in y direction
	
	(S_N_pxl,S_N_bin)=calc_signal_to_noise_ratio(N_obj_pxl,N_obj_bin,N_dark,n_bin_y,N_read,t_exp)

	return (N_obj_pxl,N_obj_bin,bin_size_pxl,bin_size_bin,flux,S_N_pxl,S_N_bin)

# Programme
#

###################################################################################################
# Informations about the program. Format of the input parameters and your utility.

system("clear")

print("NIRPS ETC - Version (May/2018)")

print("\n\n#######################################################################")
print("#                                                                     #")
print("#     INPUT parameters:                                               #")
print("#                                                                     #")
print("#     1 - Spectral Type Template                                      #")
print("#         IRTF Spectral Library.                                      #")
print("#         (http://irtfweb.ifa.hawaii.edu/~spex/                       #")
print("#          IRTF_Spectral_Library/index.html)]:                        #")
print("#         F0V/F5V/G0V/G5V/K0V/K3V/K7V/M0V/M1V/M2V/M3V/M4V/M5V         #")
print("#         M6V/M7V/M8V/M9V/L1V/L2V/L3V/L4V/L5V/L6V/L8V/T2V             #")
print("#                                                                     #")
print("#     2 - Observation mode:                                           #")
print("#         HA - High Accurace mode                                     #")
print("#         HE - High Efficience mode (whitout AO)                      #")
print("#                                                                     #")
print("#     3 - Seeing (in arcsec) @ 500nm:                                 #")
print("#         Range from 0.7 to 2.0                                       #")
print("#                                                                     #")
print("#     4 - Airmass:                                                    #")
print("#         Range from 1.00 to 5.55                                     #")
print("#                                                                     #")
print("#     5 - Object magnitude (H band)                                   #")
print("#                                                                     #")
print("#     6 - Exposure time (in seconds)                                  #")
print("#                                                                     #")
print("#     ATTENTION:                                                      #")
print("#         - This simulation is to derive the S/N at the central       #")
print("#           wavelegth of each order for the H4RG.                     #")
print("#         - The relation between S/N and RV precision is              #")
print("#           calculate only for M stars.                               #")
print("#         - RV precision is based on Figueira et al. (2016), with     #")
print("#           updated values from Jason et al., in preparation, and     #")
print("#           Artigau et al. (2018).                                    #")
print("#                                                                     #")
print("#######################################################################\n\n")

st=input_spectral_type()
obs_mode=input_observation_mode()
seeing=input_seeing()
airmass=input_airmass()
H=input_magnitude()
t_exp=input_exposure_time()

resolving_power=calc_resoving_power(obs_mode)
(flux_sts,I,Ho)=calc_flux_spectral_type(st,H)
(total_effs,wavelengths,wavelengths_nm,orders,central_wave,beg_wave,end_wave,order_wave,obs_mode)=calc_total_efficiency(obs_mode,seeing,airmass,I)
I_max=calc_I_max(obs_mode)


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

	if order_now == order:
		(n_obj_pxl,n_obj_bin,bin_size_pxl,bin_size_bin,flux,s_n_pxl,s_n_bin)=calc_S_N_order(obs_mode,t_exp,d_t,seeing,wavelength,resolving_power,z,I,total_eff,flux_st,N_dark,st,H,Ho)
		S_N_pxl.append(s_n_pxl)
		S_N_bin.append(s_n_bin)
		FLUX.append(flux)
		N_OBJ.append(n_obj_pxl)

		if wavelength_nm == 1625.003:
			SNR_pxl_H=s_n_pxl
			SNR_bin_H=s_n_bin
			sn_h=SNR_bin_H

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

system("clear")

print ("=================================================================")
print ("Y band: 980 - 1110 nm" )
print ("J band: 1200 - 1330 nm" )
print ("H band: 1510 - 1735 nm" )
print ("\nTEMPLATE STAR: %s"%(st))
print ("H (mag): %.1f"%(H))
print ("Exposure time (seconds): %.1f\n"%(t_exp))

print ("Observation Mode: %s"%(obs_mode))
print ("seeing (arcsec): %.1f"%(seeing))
print ("airmass: %.2f\n"%(airmass))

print ("Saturation limit (e-/pxl): %d\n\n"%(I_max))

print ("=================================================================")
print ("order central wave (beg-end)    Eff.    object     snr      snr         sat")
print ("              (um)                     (e-/pxl)  (ph/pxl) (ph/res elem) (%)") 

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

while i < len(order_wave):

	print('#'+repr(order_wave[i]).rjust(3),'%7.2f (%7.2f-%7.2f)'%(central_wave[i],beg_wave[i],end_wave[i]))
	print('%5.3f'%(EFF_mean[i]), '%.5e'%Decimal(N_OBJ_mean[i]), '%5.1f'%(S_N_pxl_mean[i]), '%8.1f'%(S_N_bin_mean[i]),   repr(int(SATURATION[i]*100.)).rjust(10))
	SN_pxl=SN_pxl+S_N_pxl_mean[i]
	SN_bin=SN_bin+S_N_bin_mean[i]
	if 980 <= central_wave[i] <= 1110:
		SN_pxl_Y=SN_pxl_Y+S_N_pxl_mean[i]
		mean_pxl_Y=mean_pxl_Y+1
		SN_bin_Y=SN_bin_Y+S_N_bin_mean[i]
		mean_bin_Y=mean_bin_Y+1
	if 1200 <= central_wave[i] <= 1330:
		SN_pxl_J=SN_pxl_J+S_N_pxl_mean[i]
		mean_pxl_J=mean_pxl_J+1
		SN_bin_J=SN_bin_J+S_N_bin_mean[i]
		mean_bin_J=mean_bin_J+1
	if 1510 <= central_wave[i] <= 1735:
		SN_pxl_H=SN_pxl_H+S_N_pxl_mean[i]
		mean_pxl_H=mean_pxl_H+1
		SN_bin_H=SN_bin_H+S_N_bin_mean[i]
		mean_bin_H=mean_bin_H+1
	i=i+1

print ("-----------------------------------------------------------------")
print("\n SIGNAL TO NOISE RATIO:\n")
print ("Mean S/N: %5.1f (ph/pxl) | %5.1f (ph/res elem)"%(SN_pxl/68,SN_bin/68))
print ("Mean S/N (ph/pxl): Y=%5.1f | J=%5.1f | H=%5.1f"%(SN_pxl_Y/mean_pxl_Y,SN_pxl_J/mean_pxl_J,SN_pxl_H/mean_pxl_H))
print ("         (ph/res elem): Y=%5.1f | J=%5.1f | H=%5.1f\n\n"%(SN_bin_Y/mean_bin_Y,SN_bin_J/mean_bin_J,SN_bin_H/mean_bin_H))
print ("S/N in H (1625 nm): %5.1f (ph/pxl) | %5.1f (ph/res elem)\n"%(SNR_pxl_H,SNR_bin_H))
print ("-----------------------------------------------------------------")

calc_delta_rv(st,obs_mode,sn_h)

print ("=================================================================\n\n")

################################################
#
#
#        Plot SNR vs Wavelength
#
#
################################################

fig = plt.figure(figsize=(15, 8)) 								# Size of the figure plot on the monitor

plt.plot(wavelengths_nm,S_N_pxl,label="All wavelength",color='gainsboro',zorder=1)
plt.axvspan(980, 1110, facecolor='#2ca02c', alpha=0.1)
plt.text(1045,-5,'Y')
plt.axvspan(1200, 1330, facecolor='#2ca02c', alpha=0.1)
plt.text(1265,-5,'J')
plt.axvspan(1510, 1735, facecolor='#2ca02c', alpha=0.1)
plt.text(1622.5,-5,'H')

cmap=plt.get_cmap('jet')
for (cw,snr,sat) in zip(central_wave,S_N_pxl_mean,SATURATION):

	plt.scatter(cw,snr,color=cmap(sat),zorder=2)


plt.title('S/N vs Wavelength graph', size='small')
plt.ylabel('S/N (ph/pxl)', size='small')
plt.xlabel('Wavelength (nm)', size='small')

# Make a figure and axes with dimensions as desired.
ax1 = fig.add_axes([0.91, 0.1, 0.01, 0.8])

# Set the colormap and norm to correspond to the data for which
# the colorbar will be used.
norm = mpl.colors.Normalize(vmin=0, vmax=100)

# ColorbarBase derives from ScalarMappable and puts a colorbar
# in a specified axes, so it has everything needed for a
# standalone colorbar.  There are many more kwargs, but the
# following gives a basic continuous colorbar with ticks
# and labels.
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='vertical')
cb1.set_label('Saturation (%)')

plt.xlim=(900,1900)


#plt.show()