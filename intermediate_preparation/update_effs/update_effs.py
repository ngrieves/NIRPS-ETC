#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 22 16:41:10 2021
Create new NIRPS transmission/efficiency table for Exposure Time Calculator
@author: nolangrieves@gmail.com
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
new_orders_first = np.arange(106,150)[::-1]
new_orders_second = np.arange(80,105)[::-1]
#new_orders_second = np.arange(81,105)[::-1]
#new_orders_second = np.arange(78,105)[::-1]
new_orders = np.append(new_orders_first,new_orders_second)

new_effs_filename = 'NIRPS_effs.txt'
new_waverange_filename = 'NIRPS_wave_range.txt'
new_tapas_filename = 'NIRPS_tapas.txt'
new_st_templates_filename = 'NIRPS_STAR_templates.txt'


#load new wave matrix and blaze data
wave_matrix = fits.open('NIRPS_WAVE_MATRIX_A.fits')
wave_blaze = fits.open('r.NIRPS_HA_FLAT166_0002_BLAZE_A.fits') #wave_blaze = fits.open('NIRPS_S2D_BLAZE_A.fits')

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

###############################################
############ wavelength range file ############
###############################################
wave_range_array = np.empty([len(new_orders),5])
wave_range_array[:,0] = new_orders

ncols = 33 #number of columns for the effs table
new_transmission_array = np.empty([order_array_size*len(new_orders),ncols])
for iord in range(len(new_orders)):
    #Get wavelengths from MATRIX file and Blaze from file and normalize each order
    #for each order interpolate to desired array size
    
    if wave_range == 'TSR':
        ############## TSR ############## 
        iord_newwave_init = new_wave[iord]/10.0
        iord_newwave_interp = scipy.interpolate.interp1d(np.arange(iord_newwave_init.size),iord_newwave_init)
        iord_newwave = iord_newwave_interp(np.linspace(0,iord_newwave_init.size-1,order_array_size))
    
        iord_blazenorm_init = new_blaze[iord]/max(new_blaze[iord])
        iord_blazenorm_interp = scipy.interpolate.interp1d(np.arange(iord_blazenorm_init.size),iord_blazenorm_init)
        iord_blazenorm = iord_blazenorm_interp(np.linspace(0,iord_blazenorm_init.size-1,order_array_size))
        
    if wave_range == 'FSR':
        ############## FSR ############## 
        FSRcut = 0.5 #define blaze cut of FSR
        blaze_med_iord = medfilt(new_blaze[iord],kernel_size=101)
        gt50_med = (blaze_med_iord > FSRcut*max(blaze_med_iord))
    
        iord_newwave_init = new_wave[iord,gt50_med]/10.0
        iord_newwave_interp = scipy.interpolate.interp1d(np.arange(iord_newwave_init.size),iord_newwave_init)
        iord_newwave = iord_newwave_interp(np.linspace(0,iord_newwave_init.size-1,order_array_size))
        
        iord_blazenorm_init = blaze_med_iord[gt50_med]/max(blaze_med_iord[gt50_med])
        iord_blazenorm_interp = scipy.interpolate.interp1d(np.arange(iord_blazenorm_init.size),iord_blazenorm_init)
        iord_blazenorm = iord_blazenorm_interp(np.linspace(0,iord_blazenorm_init.size-1,order_array_size))
        
    for ipix in range(order_array_size):
        new_transmission_array[(iord*order_array_size)+(ipix),0] = iord_newwave[ipix]
        new_transmission_array[(iord*order_array_size)+(ipix),1] = int(new_orders[iord])
        new_transmission_array[(iord*order_array_size)+(ipix),2] = iord_blazenorm[ipix]



    ###############################################
    ############ wavelength range file ############
    ###############################################
    wave_range_array[iord,1] = np.median(new_wave[iord]/10.0)
    wave_range_array[iord,2] = np.max(new_wave[iord]/10.0) - np.min(new_wave[iord]/10.0)    
    wave_range_array[iord,3] = np.min(new_wave[iord]/10.0)
    wave_range_array[iord,4] = np.max(new_wave[iord]/10.0)   

wave_range_columns = ['#order','central_wave','range','wave_min','wave_max']
wave_range_out = pd.DataFrame(data=wave_range_array,index=None,columns=wave_range_columns)

wave_range_out['#order'] = wave_range_out['#order'].apply(np.int64)
wave_range_out['central_wave'] = wave_range_out['central_wave'].map(lambda x: '{0:.3f}'.format(x)) 
wave_range_out['range'] = wave_range_out['range'].map(lambda x: '{0:.3f}'.format(x)) 
wave_range_out['wave_min'] = wave_range_out['wave_min'].map(lambda x: '{0:.3f}'.format(x)) 
wave_range_out['wave_max'] = wave_range_out['wave_max'].map(lambda x: '{0:.3f}'.format(x)) 
wave_range_out.to_csv('updated_files/'+new_waverange_filename,index=None,sep='\t')#,format=('%s','%.8f','%.8f','%.8f','%.8f'))


#Telescope     
newtel_interp = scipy.interpolate.interp1d(new_tel[0]*1000.0,new_tel[1],fill_value="extrapolate")
new_transmission_array[:,3] = newtel_interp(new_transmission_array[:,0])

#Front-end
new_frontend_interp = scipy.interpolate.interp1d(new_frontend[0],new_frontend[4],fill_value="extrapolate")
new_transmission_array[:,4] = new_frontend_interp(new_transmission_array[:,0])

#Back-end
new_backend_interp = scipy.interpolate.interp1d(new_backend[0],new_backend[1],fill_value="extrapolate")
new_transmission_array[:,5] = new_backend_interp(new_transmission_array[:,0])

#H4RG
new_H4RG_interp = scipy.interpolate.interp1d(new_H4RG[0],new_H4RG[1],fill_value="extrapolate")
new_transmission_array[:,6] = new_H4RG_interp(new_transmission_array[:,0])

#HA mode
new_fiber_HA_interp = scipy.interpolate.interp1d(new_fiber_HA[0],(new_fiber_HA[2]*new_fiber_HA[7]),fill_value="extrapolate")
new_transmission_array[:,7] = new_fiber_HA_interp(new_transmission_array[:,0])

#HE mode
new_fiber_HE_interp = scipy.interpolate.interp1d(new_fiber_HE[0],(new_fiber_HE[2]*new_fiber_HE[7]),fill_value="extrapolate")
new_transmission_array[:,8] = new_fiber_HE_interp(new_transmission_array[:,0])


#Fiber Coupling with Imag & Seeing
#HA
new_fiber_HA_see07_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[1],fill_value="extrapolate")
new_transmission_array[:,9] = new_fiber_HA_see07_I9_interp(new_transmission_array[:,0])
new_fiber_HA_see07_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[2],fill_value="extrapolate")
new_transmission_array[:,10] = new_fiber_HA_see07_I10_interp(new_transmission_array[:,0])
new_fiber_HA_see07_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[3],fill_value="extrapolate")
new_transmission_array[:,11] = new_fiber_HA_see07_I11_interp(new_transmission_array[:,0])
new_fiber_HA_see07_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[4],fill_value="extrapolate")
new_transmission_array[:,12] = new_fiber_HA_see07_I12_interp(new_transmission_array[:,0])

new_fiber_HA_see09_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[1],fill_value="extrapolate")
new_transmission_array[:,13] = new_fiber_HA_see09_I9_interp(new_transmission_array[:,0])
new_fiber_HA_see09_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[2],fill_value="extrapolate")
new_transmission_array[:,14] = new_fiber_HA_see09_I10_interp(new_transmission_array[:,0])
new_fiber_HA_see09_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[3],fill_value="extrapolate")
new_transmission_array[:,15] = new_fiber_HA_see09_I11_interp(new_transmission_array[:,0])
new_fiber_HA_see09_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[4],fill_value="extrapolate")
new_transmission_array[:,16] = new_fiber_HA_see09_I12_interp(new_transmission_array[:,0])

new_fiber_HA_see12_I9_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[1],fill_value="extrapolate")
new_transmission_array[:,17] = new_fiber_HA_see12_I9_interp(new_transmission_array[:,0])
new_fiber_HA_see12_I10_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[2],fill_value="extrapolate")
new_transmission_array[:,18] = new_fiber_HA_see12_I10_interp(new_transmission_array[:,0])
new_fiber_HA_see12_I11_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[3],fill_value="extrapolate")
new_transmission_array[:,19] = new_fiber_HA_see12_I11_interp(new_transmission_array[:,0])
new_fiber_HA_see12_I12_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[4],fill_value="extrapolate")
new_transmission_array[:,20] = new_fiber_HA_see12_I12_interp(new_transmission_array[:,0])

#HE
new_fiber_HE_see07_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[5],fill_value="extrapolate")
new_transmission_array[:,21] = new_fiber_HE_see07_I9_interp(new_transmission_array[:,0])
new_fiber_HE_see07_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[6],fill_value="extrapolate")
new_transmission_array[:,22] = new_fiber_HE_see07_I10_interp(new_transmission_array[:,0])
new_fiber_HE_see07_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[7],fill_value="extrapolate")
new_transmission_array[:,23] = new_fiber_HE_see07_I11_interp(new_transmission_array[:,0])
new_fiber_HE_see07_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p7[0],new_coupling_see0p7[8],fill_value="extrapolate")
new_transmission_array[:,24] = new_fiber_HE_see07_I12_interp(new_transmission_array[:,0])

new_fiber_HE_see09_I9_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[5],fill_value="extrapolate")
new_transmission_array[:,25] = new_fiber_HE_see09_I9_interp(new_transmission_array[:,0])
new_fiber_HE_see09_I10_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[6],fill_value="extrapolate")
new_transmission_array[:,26] = new_fiber_HE_see09_I10_interp(new_transmission_array[:,0])
new_fiber_HE_see09_I11_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[7],fill_value="extrapolate")
new_transmission_array[:,27] = new_fiber_HE_see09_I11_interp(new_transmission_array[:,0])
new_fiber_HE_see09_I12_interp = scipy.interpolate.interp1d(new_coupling_see0p9[0],new_coupling_see0p9[8],fill_value="extrapolate")
new_transmission_array[:,28] = new_fiber_HE_see09_I12_interp(new_transmission_array[:,0])

new_fiber_HE_see12_I9_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[5],fill_value="extrapolate")
new_transmission_array[:,29] = new_fiber_HE_see12_I9_interp(new_transmission_array[:,0])
new_fiber_HE_see12_I10_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[6],fill_value="extrapolate")
new_transmission_array[:,30] = new_fiber_HE_see12_I10_interp(new_transmission_array[:,0])
new_fiber_HE_see12_I11_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[7],fill_value="extrapolate")
new_transmission_array[:,31] = new_fiber_HE_see12_I11_interp(new_transmission_array[:,0])
new_fiber_HE_see12_I12_interp = scipy.interpolate.interp1d(new_coupling_see1p2[0],new_coupling_see1p2[8],fill_value="extrapolate")
new_transmission_array[:,32] = new_fiber_HE_see12_I12_interp(new_transmission_array[:,0])


final_columns = ['#wavelength','order','blaze','telescope','frontend','backend',
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
trans_out = trans_df[['#wavelength','telescope','frontend','backend','H4RG',
                     'HA_mode',
                     'HA_see07_I9','HA_see07_I10','HA_see07_I11','HA_see07_I12',
                     'HA_see09_I9','HA_see09_I10','HA_see09_I11','HA_see09_I12',
                     'HA_see12_I9','HA_see12_I10','HA_see12_I11','HA_see12_I12',
                     'HE_mode',
                     'HE_see07_I9','HE_see07_I10','HE_see07_I11','HE_see07_I12',
                     'HE_see09_I9','HE_see09_I10','HE_see09_I11','HE_see09_I12',
                     'HE_see12_I9','HE_see12_I10','HE_see12_I11','HE_see12_I12',
                     'blaze','order']]

#adjust format
for icol in range(len(trans_out.columns)-1):
    #print(trans_out.columns[icol])
    trans_out[trans_out.columns[icol]] = trans_out[trans_out.columns[icol]].map(lambda x: '{0:.10f}'.format(x)) 
trans_out['order'] = trans_out['order'].apply(np.int64)
               
#save dataframe
trans_out.to_csv('updated_files/'+new_effs_filename,index=None,sep='\t')


###### plots to compare old vs new effs #####
if showplots == 1:
    eff = pd.read_csv('old_files/effs.txt',header=0,sep=r"\s+",skiprows=[1])
    
    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['Order'],'.',color='red',alpha=0.7,markersize=1,label='Order (original)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,1],'.',color='blue',alpha=0.7,markersize=2,label='Order (new & interpolated)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Order',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_order_definition.png')
    
    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['Blaze'],'.',color='red',alpha=0.7,markersize=1,label='Blaze (original)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,2],'.',color='blue',alpha=0.7,markersize=2,label='Blaze (new & interploated)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Blaze',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_blaze_definition.png')

    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['Tel_clean'],'.',color='red',alpha=0.7,markersize=1,label='Tel_clean (original)')
    plt.plot(new_tel[0]*1000.0,new_tel[1],'.',color='blue',alpha=0.7,markersize=4,label='Tel_clean (new)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,3],'.',color='green',alpha=0.7,markersize=2,label='Tel_clean (new interpolate)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_telescope.png')
    
    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['Fiber-end'],'.',color='red',alpha=0.7,markersize=1,label='Fiber-End (original)')
    plt.plot(new_frontend[0],new_frontend[4],'.',color='blue',alpha=0.7,markersize=4,label='Front-End (new - from FP1)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,4],'.',color='green',alpha=0.7,markersize=2,label='Front-End (new - from FP1 interpolate)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_frontend.png')

    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['Back-End'],'.',color='red',alpha=0.7,markersize=1,label='Back-End (original)')
    plt.plot(new_backend[0],new_backend[1],'.',color='blue',alpha=0.7,markersize=4,label='Back-End (new)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,5],'.',color='green',alpha=0.7,markersize=2,label='Back-End (new interpolate)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_backend.png')

    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['H4RG'],'.',color='red',alpha=0.7,markersize=1,label='H4RG (original)')
    plt.plot(new_H4RG[0],new_H4RG[1],'.',color='blue',alpha=0.7,markersize=4,label='H4RG (new)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,6],'.',color='green',alpha=0.7,markersize=2,label='H4RG (new interpolate)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_H4RG.png')


    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['HA'],'.',color='red',alpha=0.7,markersize=1,label='HA (original)')
    plt.plot(new_fiber_HA[0],new_fiber_HA[2]*new_fiber_HA[7],'.',color='blue',alpha=0.7,markersize=4,label='Fiber throughput HA: Absorption in fiber T * Total HA measured (new)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,7],'.',color='green',alpha=0.7,markersize=2,label='Fiber throughput HA: Absorption in fiber T * Total HA measured (new interpolate)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.ylim(0.35,0.6)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=16)
    plt.tight_layout()
    plt.savefig('plots/transmission_HA.png')


    plt.figure(figsize=(12,8))
    plt.plot(eff['#Lambda'],eff['HE'],'.',color='red',alpha=0.7,markersize=1,label='HE (original)')
    plt.plot(new_fiber_HE[0],new_fiber_HE[2]*new_fiber_HE[7],'.',color='blue',alpha=0.7,markersize=4,label='Fiber throughput HE: Absorption in fiber T * Total HE measured (new)')
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,8],'.',color='green',alpha=0.7,markersize=2,label='Fiber throughput HE: Absorption in fiber T * Total HE measured (new interpolate)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.ylim(0.42,0.65)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=16)
    plt.tight_layout()
    plt.savefig('plots/transmission_HE.png')
        

    HA_Iseeing = ['0.7-I=9', '0.7-I=10', '0.7-I=11', '0.7-I=12', '0.9-I=9',
                 '0.9-I=10','0.9-I=11', '0.9-I=12', '1.2-I=9', '1.2-I=10', 
                 '1.2-I=11', '1.2-I=12']
    plt.figure(figsize=(12,8))
    plt.title('HA Fiber Coupling Transmission with Seeing & Imag',fontsize=20)
    for i in range(len(HA_Iseeing)):
        plt.plot(eff['#Lambda'],eff[HA_Iseeing[i]],'.',alpha=0.7,label=HA_Iseeing[i])
        plt.plot((new_transmission_array[:,0]),new_transmission_array[:,9+i],'.',color='black',alpha=0.7,markersize=1)
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,9+i],'.',color='black',alpha=0.7,markersize=1,label='(new & interpolated)')
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig('plots/transmission_HA_fiber_coupling.png')
    
    
    HE_Iseeing = ['0.7-I=9.1', '0.7-I=10.1', '0.7-I=11.1', '0.7-I=12.1',
                  '0.9-I=9.1', '0.9-I=10.1', '0.9-I=11.1', '0.9-I=12.1',
                  '1.2-I=9.1','1.2-I=10.1', '1.2-I=11.1', '1.2-I=12.1']
    plt.figure(figsize=(12,8))
    plt.title('HE Fiber Coupling Transmission with Seeing & Imag',fontsize=20)
    for i in range(len(HA_Iseeing)):
        plt.plot(eff['#Lambda'],eff[HE_Iseeing[i]],'.',alpha=0.7,label=HE_Iseeing[i])
        plt.plot((new_transmission_array[:,0]),new_transmission_array[:,21+i],'.',color='black',alpha=0.7,markersize=1)
    plt.plot((new_transmission_array[:,0]),new_transmission_array[:,21+i],'.',color='black',alpha=0.7,markersize=1,label='(new & interpolated)')                 
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(900,2000)
    plt.xlabel('Wavelength [nm]',fontsize=20)
    plt.ylabel('Transmission',fontsize=20)
    plt.legend(fontsize=18)
    plt.tight_layout()
    plt.savefig('plots/transmission_HE_fiber_coupling.png')

##########################################################################  
#### resample tapas file to match new length of effs #####
##########################################################################  
tapas = pd.read_csv('old_files/tapas.txt',header=0,sep=r"\s+",skiprows=[1])
ncols_tapas = 8
tapas_new = np.empty([order_array_size*len(new_orders),ncols_tapas])
tapas_new[:,0] = new_transmission_array[:,0]
tapas_new1_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_1.00'],fill_value="extrapolate")
tapas_new[:,1] = tapas_new1_interp(new_transmission_array[:,0])
tapas_new2_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_1.02'],fill_value="extrapolate")
tapas_new[:,2] = tapas_new2_interp(new_transmission_array[:,0])
tapas_new3_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_1.06'],fill_value="extrapolate")
tapas_new[:,3] = tapas_new3_interp(new_transmission_array[:,0])
tapas_new4_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_1.15'],fill_value="extrapolate")
tapas_new[:,4] = tapas_new4_interp(new_transmission_array[:,0])
tapas_new5_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_1.31'],fill_value="extrapolate")
tapas_new[:,5] = tapas_new5_interp(new_transmission_array[:,0])
tapas_new6_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_1.56'],fill_value="extrapolate")
tapas_new[:,6] = tapas_new6_interp(new_transmission_array[:,0])
tapas_new7_interp = scipy.interpolate.interp1d(tapas['#Lambda'],tapas['AM_2.00'],fill_value="extrapolate")
tapas_new[:,7] = tapas_new7_interp(new_transmission_array[:,0])
mask0 = (tapas_new < 0)
tapas_new[mask0] = 0

tapas_df = pd.DataFrame(data=tapas_new,index=None,columns=tapas.columns)
tapas_df.to_csv('updated_files/'+new_tapas_filename,index=None,sep='\t')
##########################################################################  
########################################################################## 

##########################################################################  
#### resample stellar templates file to match new length of effs #####
########################################################################## 
st_templates = pd.read_csv('old_files/STAR_templates.txt',header=0,sep=r"\s+")
st_templates_new = np.empty([order_array_size*len(new_orders),26])
st_templates_new[:,0] = new_transmission_array[:,0]

st_templates_columns = ['f0v', 'f5v', 'g0v', 'g5v', 'k0v', 'k3v', 'k7v', 
                        'm0v','m1v', 'm2v', 'm3v', 'm4v', 'm5v', 'm6v', 'm7v', 'm8v', 
                        'm9v', 'l1v','l2v', 'l3v', 'l5v', 'l6v', 'l8v', 't2v', 'g8v']
for icolst in range(len(st_templates_columns)):
    st_templates_interp_icol = scipy.interpolate.interp1d(st_templates['#Lambda'],st_templates[st_templates_columns[icolst]],fill_value="extrapolate")
    st_templates_new[:,icolst+1] = st_templates_interp_icol(new_transmission_array[:,0])

st_templates_new_df = pd.DataFrame(data=st_templates_new,index=None,columns=st_templates.columns)
st_templates_new_df.to_csv('updated_files/'+new_st_templates_filename,index=None, sep='\t')

#plt.figure()
#plt.plot(st_templates['#Lambda'],st_templates['m4v'],'.',label='original')
#plt.plot(st_templates_new[:,0],st_templates_new[:,12],'.',label='interp')
#plt.legend()
#plt.xlim(1350,1400)
##plt.ylim(0.8e-11,1e-11)


#plt.plot(tapas['#Lambda'],tapas['AM_1.02'],'.',label='original')
#plt.plot(tapas_new[:,0],tapas_new[:,1],'.',label='interp')
#plt.legend()
#plt.xlim(960,1000)
#plt.ylim(0.8e-11,1e-11)

##########################################################################  
########################################################################## 
    