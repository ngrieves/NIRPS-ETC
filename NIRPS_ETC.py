from nirps_etc_driver import run_nirps_etc
import time

start_time = time.time()


###############################################################
#### CHANGE INPUTS HERE ####

# Observation Mode (HA/HE):
obs_mode = 'HA'
# Seeing, in arcsec (range 0.7-1.2):
seeing = 1.2
# Airmass (range 1.0-2.0):
airmass = 1.2
# Object magnitude (H band):
H = 3.5
# Exposure time (in sec):
t_exp = 60

# ST Models for IRFT
# Spectral type that only have SNR estimate:
# B3V/B8V/B9V/A1V/F0V/F5V/G0V/G5V/G8V/K0V/L1V/L2V/L3V/L5V/L6V/L8V/T2V
# Spectral type that have SNR and RV precision estimate:
# K3V/K7V/M0V/M1V/M2V/M3V/M4V/M5V/M6V/M7V/M8V/M9V
st = 'A1V'

# wavelength bandpasses for YJH ('CFHT' or 'Eniric')
bandpass = 'CFHT'

# Use all order information (TSR) or only FSR
waveselect = 'FSR'
# FSR = Free Spectral Range calculated from 50% of blaze
# TSR = Total Spectral Range

# Vsini of star for Eniric Pheonix spectra RV calculations will be for [0.1,1.0,5.0,10.0] km/s

###############################################################
###############################################################
###############################################################

run_nirps_etc(obs_mode, st, H, seeing, airmass, t_exp,
              bandpass, waveselect, plot=True)

print("--- %s seconds ---" % (time.time() - start_time))
