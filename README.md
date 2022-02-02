# NIRPS-ETC
### Exposure Time Calculator (ETC) and radial velocity precision estimator for the Near InfraRed Planet Searcher (NIRPS) spectrograph

February 2022 - Before NIRPS on sky

Original NIRPS ETC code by Bruno L. Canto Martins in 2018-2019 (NIRPS-ETC/intermediate_preparation/ETC_v3.0_CantoMartins/)

Additional edits by Nolan Grieves (University of Geneva) from 2020-2022

### Overview 
* The NIRPS ETC uses spectra from the NASA Infrared Telescope Facility (IRTF) as SEDs to get estimated flux values for different spectral types: http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/
* The ETC calculates efficiency at different wavelengths using seeing, atmospheric efficiency from TAPAS (http://cds-espri.ipsl.fr/tapas/), and the measured global efficiency of the instrument
* The signal to noise ratio (SNR) at each pixel or bin is calculated from the fiber diameter, sampling, readout noise, resolution, efficiency, and flux in the pixel or bin from the IRTF template (flux=(10.\**(0.4\*(Ho-H)))\*flux_st)
* RV precisions are calculated using equation 12 of Bouchy et al. (2001: https://ui.adsabs.harvard.edu/abs/2001A%26A...374..733B/abstract), dRV=c/(Q\*sqrt(Ne-). The quality factors Q for spectra are calculated with ENIRIC from Phoenix simulated spectra or from spectral templates from the Spirou spectrograph 
  - -> see: NIRPS-ETC/intermediate_preparation/update_RV_estimates/README_update_RV_estimates

### Files
* inputs/
* intermediate_preparation/
* outputs/
* NIRPS_ETC.py
* NIRPS_ETC_script.py
* etc_targets_input.txt
* etc_targets_output.txt
* nirps_etc_lib.py

### Use
$ python NIRPS_ETC.py
* change observing options within the code at the top
  - Observation Mode (HA/HE)
  - Seeing, in arcsec (range 0.7-1.2)
  - Airmass (range 1.0-2.0)
  - Object magnitude (H band)
  - Exposure time (in sec)
  - Spectral type (F0V/F5V/G0V/G5V/G8V/K0V/K3V/K7V/M0V/M1V/M2V/M3V/M4V/M5V/M6V/M7V/M8V/M9V/L1V/L2V/L3V/L4V/L5V/L6V/L8V/T2V)
  - bandpass ('CFHT' or 'Eniric') #YJH bandpasses that will affect the range of the spectra used to calculate RV precision

OR use script version:
$ python NIRPS_ETC_script.py
* change inputs for each target in a space separated text file with columns:
  - target st obs_mode seeing airmass H t_exp bandpass
* change input and output text files within code to desired option


