#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 10:09:01 2023
Create NIRPS efficiency file in format that can be used for ESO ETC
a) separate the order-independent columns into one, similar file
(but with no order column) and monotonely increasing wavelengths and
b) keep the order-dependent column(s) and the order column in a second
file (with wavelengths monotonely increasing for each order)
@author: nolangrieves
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy.interpolate
from scipy.signal import medfilt


showplots = 1
order_array_size = 1000 #sample size that each order will be interpolated into
wave_range = 'TSR' #select 'TSR' or 'FSR'
#FSR = Free Spectral Range calculated from 50% of blaze
#TSR = Total Spectral Range

#set the physical orders - make sure to match the wave matrix and blaze files
new_orders_first = np.arange(107,150)[::-1]
new_orders_second = np.arange(77,105)[::-1]
new_orders = np.append(new_orders_first,new_orders_second)

#load new wave matrix and blaze data
wave_matrix = fits.open('dec2022_files/r.NIRPS.2022-12-08T22:56:24.084_WAVE_MATRIX_THAR_FP_A.fits') #HA mode from Dec 2022
wave_blaze = fits.open('dec2022_files/r.NIRPS.2022-12-08T22:06:25.499_BLAZE_A.fits') #HA mode from Dec 2022 

new_wave = wave_matrix[1].data
new_blaze = wave_blaze[1].data

#load other new tranmsission data
new_direc = 'transmission_values_v34/'
new_tel = pd.read_csv(new_direc+'telescope_trans.txt',header=None,sep=r"\s+",skiprows=[0])
new_H4RG = pd.read_csv(new_direc+'H4RG_QE.txt',header=None,sep=r"\s+",skiprows=[0])         
new_fiber_HA = pd.read_csv(new_direc+'fiber_throughput_HA.txt',header=None,sep=r"\s+",skiprows=[0])
new_atm = pd.read_csv(new_direc+'atm_trans_tapas.txt',header=None,sep=r"\s+",skiprows=[0])
new_fiber_HE = pd.read_csv(new_direc+'fiber_throughput_HE.txt',header=None,sep=r"\s+",skiprows=[0])
new_backend = pd.read_csv(new_direc+'backend_spectrometer_transmission.txt',header=None,sep=r"\s+",skiprows=[0])
new_frontend = pd.read_csv(new_direc+'frontend_trans.txt',header=None,sep=r"\s+",skiprows=[0])
new_coupling_see0p7 = pd.read_csv(new_direc+'coupling_see0p7.txt',header=None,sep=r"\s+",skiprows=[0])
new_coupling_see0p9 = pd.read_csv(new_direc+'coupling_see0p9.txt',header=None,sep=r"\s+",skiprows=[0])
new_coupling_see1p2 = pd.read_csv(new_direc+'coupling_see1p2.txt',header=None,sep=r"\s+",skiprows=[0])

########## FOR GRATING 2 DECEMBER 2022 ##########
new_backend[1] = new_backend[1]*1.40
#need to do at interpolation step for fiber HE new_fiber_HE*0.88 
#################################################

###############################################
############ wavelength range file ############
###############################################
#wave_range_array = np.empty([len(new_orders),5])
#wave_range_array[:,0] = new_orders

#ncols = 33 #number of columns for the effs table
#new_transmission_array = np.empty([order_array_size*len(new_orders),ncols])
new_wln_ord_eff = np.empty([order_array_size*len(new_orders),3])
for iord in range(len(new_orders)):
    #Get wavelengths from MATRIX file and Blaze from file and normalize each order
    #for each order interpolate to desired array size
    
    if wave_range == 'TSR':
        ############## TSR ############## 
        iord_newwave_init = new_wave[iord]/10.0
        iord_newwave_interp = scipy.interpolate.interp1d(np.arange(iord_newwave_init.size),iord_newwave_init)
        iord_newwave = iord_newwave_interp(np.linspace(0,iord_newwave_init.size-1,order_array_size))
    
        iord_blazenorm_init =  new_blaze[iord] / np.nanmax(new_blaze[iord]) #(new_blaze[iord] - np.nanmin(new_blaze[iord])) / (np.nanmax(new_blaze[iord]) - np.nanmin(new_blaze[iord])) #new_blaze[iord]/max(new_blaze[iord])
        iord_blazenorm_interp = scipy.interpolate.interp1d(np.arange(iord_blazenorm_init.size),iord_blazenorm_init)
        iord_blazenorm = iord_blazenorm_interp(np.linspace(0,iord_blazenorm_init.size-1,order_array_size))
    else:
        print('FSR not updated below')

    #don't use FSR for now
    # if wave_range == 'FSR':
    #     ############## FSR ############## 
    #     FSRcut = 0.5 #define blaze cut of FSR
    #     blaze_med_iord = medfilt(new_blaze[iord],kernel_size=101)
    #     gt50_med = (blaze_med_iord > FSRcut*max(blaze_med_iord))
    
    #     iord_newwave_init = new_wave[iord,gt50_med]/10.0
    #     iord_newwave_interp = scipy.interpolate.interp1d(np.arange(iord_newwave_init.size),iord_newwave_init)
    #     iord_newwave = iord_newwave_interp(np.linspace(0,iord_newwave_init.size-1,order_array_size))
        
    #     iord_blazenorm_init = (blaze_med_iord[gt50_med] - np.nanmin(blaze_med_iord[gt50_med])) / (np.nanmax(blaze_med_iord[gt50_med]) - np.nanmin(blaze_med_iord[gt50_med])) #blaze_med_iord[gt50_med]/max(blaze_med_iord[gt50_med])
    #     iord_blazenorm_interp = scipy.interpolate.interp1d(np.arange(iord_blazenorm_init.size),iord_blazenorm_init)
    #     iord_blazenorm = iord_blazenorm_interp(np.linspace(0,iord_blazenorm_init.size-1,order_array_size))
        
    for ipix in range(order_array_size):
        ### wln array for separate files
        new_wln_ord_eff[(iord*order_array_size)+(ipix),0] = iord_newwave[ipix]
        new_wln_ord_eff[(iord*order_array_size)+(ipix),1] = int(new_orders[iord])
        new_wln_ord_eff[(iord*order_array_size)+(ipix),2] = iord_blazenorm[ipix]

##################################################################
################ SAVE ORDER DEPENDENT INFORMATION ################
##################################################################
new_wln_ord_eff_columns = ['wavelength','order','blaze']

#transfer to data frame
new_wln_ord_eff_df = pd.DataFrame(data=new_wln_ord_eff,index=None,columns=new_wln_ord_eff_columns)
new_wln_ord_eff_df = new_wln_ord_eff_df[['wavelength','blaze','order']]
#adjust format
for icol in range(len(new_wln_ord_eff_df.columns)-1):
    #print(trans_out.columns[icol])
    new_wln_ord_eff_df[new_wln_ord_eff_df.columns[icol]] = new_wln_ord_eff_df[new_wln_ord_eff_df.columns[icol]].map(lambda x: '{0:.10f}'.format(x)) 
new_wln_ord_eff_df['order'] = new_wln_ord_eff_df['order'].apply(np.float64)
new_wln_ord_eff_df['order'] = new_wln_ord_eff_df['order'].apply(np.int64)
               
#save dataframe
new_wln_ord_eff_df.to_csv('mono_jan2023/effs_order.txt',index=None,sep='\t')
##################################################################
##################################################################


########### Make non wavelength dependent effs file ###########

ncols = 31 #number of columns for the effs table
new_transmission_array = np.empty([order_array_size*len(new_orders),ncols])

wave_mono = np.linspace(np.nanmin(new_wave)/10.0,np.nanmax(new_wave)/10.0,(order_array_size*len(new_orders)))
#         iord_newwave_init = new_wave[iord]/10.0
#         iord_newwave_interp = scipy.interpolate.interp1d(np.arange(iord_newwave_init.size),iord_newwave_init)
#         iord_newwave = iord_newwave_interp(np.linspace(0,iord_newwave_init.size-1,order_array_size))

#set mono wavelength
new_transmission_array[:,0] = wave_mono


#Telescope     
newtel_interp = scipy.interpolate.interp1d(new_tel[0]*1000.0,new_tel[1],fill_value="extrapolate")
new_transmission_array[:,1] = newtel_interp(new_transmission_array[:,0])

#Front-end
new_frontend_interp = scipy.interpolate.interp1d(new_frontend[0],new_frontend[4],fill_value="extrapolate")
new_transmission_array[:,2] = new_frontend_interp(new_transmission_array[:,0])

#Back-end
new_backend_interp = scipy.interpolate.interp1d(new_backend[0],new_backend[1],fill_value="extrapolate")
new_transmission_array[:,3] = new_backend_interp(new_transmission_array[:,0])

#H4RG
new_H4RG_interp = scipy.interpolate.interp1d(new_H4RG[0],new_H4RG[1],fill_value="extrapolate")
new_transmission_array[:,4] = new_H4RG_interp(new_transmission_array[:,0])

#HA mode
new_fiber_HA_interp = scipy.interpolate.interp1d(new_fiber_HA[0],(new_fiber_HA[2]*new_fiber_HA[7]),fill_value="extrapolate")
new_transmission_array[:,5] = new_fiber_HA_interp(new_transmission_array[:,0])

#HE mode
new_fiber_HE_interp = scipy.interpolate.interp1d(new_fiber_HE[0],(new_fiber_HE[2]*new_fiber_HE[7]*0.88),fill_value="extrapolate")
new_transmission_array[:,6] = new_fiber_HE_interp(new_transmission_array[:,0])


#Fiber Coupling with Imag & Seeing
#HA
new_fiber_HA_see07_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[1],fill_value="extrapolate")
new_transmission_array[:,7] = new_fiber_HA_see07_I9_interp(new_transmission_array[:,0])
new_fiber_HA_see07_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[2],fill_value="extrapolate")
new_transmission_array[:,8] = new_fiber_HA_see07_I10_interp(new_transmission_array[:,0])
new_fiber_HA_see07_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[3],fill_value="extrapolate")
new_transmission_array[:,9] = new_fiber_HA_see07_I11_interp(new_transmission_array[:,0])
new_fiber_HA_see07_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[4],fill_value="extrapolate")
new_transmission_array[:,10] = new_fiber_HA_see07_I12_interp(new_transmission_array[:,0])

new_fiber_HA_see09_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[1],fill_value="extrapolate")
new_transmission_array[:,11] = new_fiber_HA_see09_I9_interp(new_transmission_array[:,0])
new_fiber_HA_see09_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[2],fill_value="extrapolate")
new_transmission_array[:,12] = new_fiber_HA_see09_I10_interp(new_transmission_array[:,0])
new_fiber_HA_see09_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[3],fill_value="extrapolate")
new_transmission_array[:,13] = new_fiber_HA_see09_I11_interp(new_transmission_array[:,0])
new_fiber_HA_see09_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[4],fill_value="extrapolate")
new_transmission_array[:,14] = new_fiber_HA_see09_I12_interp(new_transmission_array[:,0])

new_fiber_HA_see12_I9_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[1],fill_value="extrapolate")
new_transmission_array[:,15] = new_fiber_HA_see12_I9_interp(new_transmission_array[:,0])
new_fiber_HA_see12_I10_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[2],fill_value="extrapolate")
new_transmission_array[:,16] = new_fiber_HA_see12_I10_interp(new_transmission_array[:,0])
new_fiber_HA_see12_I11_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[3],fill_value="extrapolate")
new_transmission_array[:,17] = new_fiber_HA_see12_I11_interp(new_transmission_array[:,0])
new_fiber_HA_see12_I12_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[4],fill_value="extrapolate")
new_transmission_array[:,18] = new_fiber_HA_see12_I12_interp(new_transmission_array[:,0])

#HE
new_fiber_HE_see07_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[5],fill_value="extrapolate")
new_transmission_array[:,19] = new_fiber_HE_see07_I9_interp(new_transmission_array[:,0])
new_fiber_HE_see07_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[6],fill_value="extrapolate")
new_transmission_array[:,20] = new_fiber_HE_see07_I10_interp(new_transmission_array[:,0])
new_fiber_HE_see07_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[7],fill_value="extrapolate")
new_transmission_array[:,21] = new_fiber_HE_see07_I11_interp(new_transmission_array[:,0])
new_fiber_HE_see07_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[8],fill_value="extrapolate")
new_transmission_array[:,22] = new_fiber_HE_see07_I12_interp(new_transmission_array[:,0])

new_fiber_HE_see09_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[5],fill_value="extrapolate")
new_transmission_array[:,23] = new_fiber_HE_see09_I9_interp(new_transmission_array[:,0])
new_fiber_HE_see09_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[6],fill_value="extrapolate")
new_transmission_array[:,24] = new_fiber_HE_see09_I10_interp(new_transmission_array[:,0])
new_fiber_HE_see09_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[7],fill_value="extrapolate")
new_transmission_array[:,25] = new_fiber_HE_see09_I11_interp(new_transmission_array[:,0])
new_fiber_HE_see09_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[8],fill_value="extrapolate")
new_transmission_array[:,26] = new_fiber_HE_see09_I12_interp(new_transmission_array[:,0])

new_fiber_HE_see12_I9_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[5],fill_value="extrapolate")
new_transmission_array[:,27] = new_fiber_HE_see12_I9_interp(new_transmission_array[:,0])
new_fiber_HE_see12_I10_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[6],fill_value="extrapolate")
new_transmission_array[:,28] = new_fiber_HE_see12_I10_interp(new_transmission_array[:,0])
new_fiber_HE_see12_I11_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[7],fill_value="extrapolate")
new_transmission_array[:,29] = new_fiber_HE_see12_I11_interp(new_transmission_array[:,0])
new_fiber_HE_see12_I12_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[8],fill_value="extrapolate")
new_transmission_array[:,30] = new_fiber_HE_see12_I12_interp(new_transmission_array[:,0])



final_columns = ['wavelength','telescope','frontend','backend',
                 'H4RG','HA_mode','HE_mode',
                 'HA_see07_I9','HA_see07_I10','HA_see07_I11','HA_see07_I12',
                 'HA_see09_I9','HA_see09_I10','HA_see09_I11','HA_see09_I12',
                 'HA_see12_I9','HA_see12_I10','HA_see12_I11','HA_see12_I12',
                 'HE_see07_I9','HE_see07_I10','HE_see07_I11','HE_see07_I12',
                 'HE_see09_I9','HE_see09_I10','HE_see09_I11','HE_see09_I12',
                 'HE_see12_I9','HE_see12_I10','HE_see12_I11','HE_see12_I12']

#transfer to data frame
trans_df = pd.DataFrame(data=new_transmission_array,index=None,columns=final_columns)

#make columns the same order as the original effs.txt document
trans_out = trans_df[['wavelength','telescope','frontend','backend','H4RG',
                     'HA_mode',
                     'HA_see07_I9','HA_see07_I10','HA_see07_I11','HA_see07_I12',
                     'HA_see09_I9','HA_see09_I10','HA_see09_I11','HA_see09_I12',
                     'HA_see12_I9','HA_see12_I10','HA_see12_I11','HA_see12_I12',
                     'HE_mode',
                     'HE_see07_I9','HE_see07_I10','HE_see07_I11','HE_see07_I12',
                     'HE_see09_I9','HE_see09_I10','HE_see09_I11','HE_see09_I12',
                     'HE_see12_I9','HE_see12_I10','HE_see12_I11','HE_see12_I12']]

#adjust format
for icol in range(len(trans_out.columns)-1):
    #print(trans_out.columns[icol])
    trans_out[trans_out.columns[icol]] = trans_out[trans_out.columns[icol]].map(lambda x: '{0:.10f}'.format(x)) 

#save dataframe
trans_out.to_csv('mono_jan2023/effs_wln.txt',index=None,sep='\t')

plt.figure()
plt.plot(trans_df['wavelength'],trans_df['telescope'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['frontend'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['backend'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['H4RG'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['HA_mode'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['HE_mode'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['HA_see09_I10'])
plt.figure()
plt.plot(trans_df['wavelength'],trans_df['HE_see12_I11'])

def isMonotonic(A):
    return (all(A[i] < A[i + 1] for i in range(len(A) - 1)))# or
            #all(A[i] >= A[i + 1] for i in range(len(A) - 1)))

print(isMonotonic(trans_df['wavelength']))
