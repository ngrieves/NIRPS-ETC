# Operations related to ETC calculation execution
import numpy as np
import pandas as pd
from nirps_etc_lib import reading_table

### Physical Constants ###
h = 6.62606957e-34                 # Plank's constant (J.s)
c = 299792458.                     # Light speed (m/s)

### Telescope parameters ###
#airmass=1.0                     # Airmass
# WARN: d_t is defined here but also not used in function that uses it I think, to double check
d_t = 3.57                         # Telescope diameter (m) [https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html]


### Functions for science calculations ###
def calc_fiber_diameter(obs_mode):
    if obs_mode == "HA":
        fiber_diameter = 0.4  # Fiber diameter for HA mode (arcsec)
    else:
        fiber_diameter = 0.9  # Fiber diameter for HE mode (arcsec)
    return fiber_diameter


def calc_resoving_power(obs_mode):
    if obs_mode == "HA":
        resolving_power = 80000  # Resolving power for HA mode
    else:
        resolving_power = 70000  # Resolving power for HE mode
    return resolving_power


def calc_I_max(obs_mode):
    if obs_mode == "HA":
        I_max = 3.6 * 27300  # Saturation limit for HA mode
    else:
        I_max = 16.7 * 27300  # Saturation limit for HE mode
    return I_max


def calc_sampling(obs_mode):
    if obs_mode == "HA":
        sampling = 3.8  # Sampling for HA mode (px)
    else:
        sampling = 4.3  # Sampling for HE mode (px)
    return sampling


def calc_N_pix_Y(obs_mode):
    if obs_mode == "HA":
        N_pix_Y = 3.6  # Bin for HA mode
    else:
        N_pix_Y = 16.7  # Bin for HE mode
    return N_pix_Y


def calc_flux_spectral_type(st, H, st_templates_file):
    st_templates_in = pd.read_csv('inputs/'+st_templates_file,
                                  header=0, sep=r"\s+")
    # Spectral types B3V-A1V,G8V added by NG 09.06.2022
    # I-H calculated from: http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
    # (intermediate_preparation/add_stellar_templates/I_Hmag_mamjek.py)
    st_template_mag_conversions = {
        # B3V-A1V spectra in CRIRES+ spectrophotometric standard stars
        # See star names next to dict entry.
        # (https://github.com/ivh/cr2rep/tree/master/catalogs/stdstar)
        'B3V': {'I-H': -0.279, 'Ho': 4.840},  # HR3454
        'B8V': {'I-H': -0.149, 'Ho': 3.53},  # HR8634
        'B9V': {'I-H': -0.076, 'Ho': 4.845},  # HR4468
        'A1V': {'I-H': 0.016, 'Ho': 3.71},  # HR4468
        'F0V': {'I-H': 0.348, 'Ho': 6.998},
        'F5V': {'I-H': 0.519, 'Ho': 4.738},
        'G0V': {'I-H': 0.699, 'Ho': 2.905},
        'G5V': {'I-H': 0.821, 'Ho': 4.614},
        # G8V: H mag of star in IRFT spectrum header:
        # http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/Data/G8V_HD75732.txt
        'G8V': {'I-H': 0.902, 'Ho': 4.265},
        'K0V': {'I-H': 1.009, 'Ho': 4.803},
        'K3V': {'I-H': 1.275, 'Ho': 3.469},
        'K7V': {'I-H': 1.665, 'Ho': 5.499},
        'M0V': {'I-H': 1.739, 'Ho': 5.843},
        'M1V': {'I-H': 1.799, 'Ho': 4.393},
        'M2V': {'I-H': 1.833, 'Ho': 3.640},
        'M3V': {'I-H': 1.928, 'Ho': 4.843},
        'M4V': {'I-H': 2.137, 'Ho': 6.627},
        'M5V': {'I-H': 2.354, 'Ho': 8.014},
        'M6V': {'I-H': 2.848, 'Ho': 6.482},
        'M7V': {'I-H': 3.147, 'Ho': 9.201},
        'M8V': {'I-H': 3.668, 'Ho': 11.066},
        'M9V': {'I-H': 3.996, 'Ho': 8.905},
        'L1V': {'I-H': 5.477, 'Ho': 12.041},
        'L2V': {'I-H': 5.704, 'Ho': 12.392},
        'L3V': {'I-H': 5.931, 'Ho': 12.380},
        'L5V': {'I-H': 6.385, 'Ho': 11.895},
        'L6V': {'I-H': 6.612, 'Ho': 13.099},
        'L8V': {'I-H': 7.066, 'Ho': 12.204},
        'T2V': {'I-H': 9.108, 'Ho': 14.090}
    }

    if st in st_template_mag_conversions:
        flux_sts = st_templates_in[st]
        I = H + st_template_mag_conversions[st]['I-H']
        Ho = st_template_mag_conversions[st]['Ho']

    else:
        # TODO: This function was called but not implemented. Raising an error instead for now
        # st = calc_spectral_type()
        raise ValueError(f"There is no error for spectral type {st}")

    return (flux_sts, I, Ho)


def calc_total_efficiency(obs_mode, seeing, airmass, I, wave_range_file,
                          tapas_file, effs_file):
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

    central_wave = []
    beg_wave = []
    end_wave = []
    order_wave = []

    for line in lines:
        central_wave.append(line[1])
        beg_wave.append(line[3])
        end_wave.append(line[4])
        order_wave.append(int(line[0]))

    file_name = 'inputs/' + tapas_file
    lines = reading_table(file_name)

    ATM_EFF = []

    for line in lines:

        if airmass == 1.0:
            eff = line[1]
            ATM_EFF.append(eff)

        elif airmass > 1.0 and airmass < 1.02:
            f = line[1]+(line[2]-line[1])*(airmass-1.00)/0.02
            eff = f
            ATM_EFF.append(eff)

        elif airmass == 1.02:
            eff = line[2]
            ATM_EFF.append(eff)

        elif airmass > 1.02 and airmass < 1.06:
            f = line[2]+(line[3]-line[2])*(airmass-1.02)/0.04
            eff = f
            ATM_EFF.append(eff)

        elif airmass == 1.06:
            eff = line[3]
            ATM_EFF.append(eff)

        elif airmass > 1.06 and airmass < 1.15:
            f = line[3]+(line[4]-line[3])*(airmass-1.06)/0.09
            eff = f
            ATM_EFF.append(eff)

        elif airmass == 1.15:
            eff = line[4]
            ATM_EFF.append(eff)

        elif airmass > 1.15 and airmass < 1.31:
            f = line[4]+(line[5]-line[4])*(airmass-1.15)/0.16
            eff = f
            ATM_EFF.append(eff)

        elif airmass == 1.31:
            eff = line[5]
            ATM_EFF.append(eff)

        elif airmass > 1.31 and airmass < 1.56:
            f = line[5]+(line[6]-line[5])*(airmass-1.31)/0.25
            eff = f
            ATM_EFF.append(eff)

        elif airmass == 1.56:
            eff = line[6]
            ATM_EFF.append(eff)

        elif airmass > 1.56 and airmass < 2.00:
            f = line[6]+(line[7]-line[6])*(airmass-1.56)/0.44
            eff = f
            ATM_EFF.append(eff)

        elif airmass == 2.00:
            eff = line[7]
            ATM_EFF.append(eff)

    file_name = 'inputs/'+effs_file
    lines = reading_table(file_name)

    EFF = []
    wavelengths = []
    wavelengths_nm = []
    orders = []

    i = 0

    for line in lines:
        eff = line[1]*line[2]*line[3]*line[4]*line[31]*ATM_EFF[i]
        wavelengths.append(line[0]*10**(-9))  # Transform wavelength from nm to m
        wavelengths_nm.append(line[0])
        orders.append(line[32])

        i = i+1

        if obs_mode == 'HA':
            eff = eff*line[5]

            if I <= 9:

                if seeing == 0.7:
                    eff = eff*line[6]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[6]+(line[10]-line[6])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[10]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[10]+(line[14]-line[10])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[14]
                    EFF.append(eff)

            elif I > 9 and I <10:

                if seeing == 0.7:
                    f = line[6]+(line[7]-line[6])*(I-9.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f_07 = line[6]+(line[7]-line[6])*(I-9.)/1.
                    f_09 = line[10]+(line[11]-line[10])*(I-9.)/1.
                    f = f_07+(f_09-f_07)*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    f = line[10]+(line[11]-line[10])*(I-9.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f_09 = line[10]+(line[11]-line[10])*(I-9.)/1.
                    f_12 = line[14]+(line[15]-line[14])*(I-9.)/1.
                    f = f_09+(f_12-f_09)*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    f = line[14]+(line[15]-line[14])*(I-9.)/1.
                    eff = eff*f
                    EFF.append(eff)

            elif I == 10:

                if seeing == 0.7:
                    eff = eff*line[7]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[7]+(line[11]-line[7])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[11]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[11]+(line[15]-line[11])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[15]
                    EFF.append(eff)

            elif I > 10 and I < 11:

                if seeing == 0.7:
                    f = line[7]+(line[8]-line[7])*(I-10.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f_07 = line[7]+(line[8]-line[7])*(I-10.)/1.
                    f_09 = line[11]+(line[12]-line[11])*(I-10.)/1.
                    f = f_07+(f_09-f_07)*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    f = line[11]+(line[12]-line[11])*(I-10.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f_09 = line[11]+(line[12]-line[11])*(I-10.)/1.
                    f_12 = line[15]+(line[16]-line[15])*(I-10.)/1.
                    f = f_09+(f_12-f_09)*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    f = line[15]+(line[16]-line[15])*(I-10.)/1.
                    eff = eff*f
                    EFF.append(eff)

            elif I == 11:

                if seeing == 0.7:
                    eff = eff*line[8]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[8]+(line[12]-line[8])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[12]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[12]+(line[16]-line[12])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[16]
                    EFF.append(eff)

            elif I > 11 and I < 12:

                if seeing == 0.7:
                    f = line[8]+(line[9]-line[8])*(I-11.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f_07 = line[8]+(line[9]-line[8])*(I-11.)/1.
                    f_09 = line[12]+(line[13]-line[12])*(I-11.)/1.
                    f = f_07+(f_09-f_07)*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    f = line[12]+(line[13]-line[12])*(I-11.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f_09 = line[12]+(line[13]-line[12])*(I-11.)/1.
                    f_12 = line[16]+(line[17]-line[16])*(I-11.)/1.
                    f = f_09+(f_12-f_09)*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    f = line[16]+(line[17]-line[16])*(I-11.)/1.
                    eff = eff*f
                    EFF.append(eff)

            elif I >= 12:

                if seeing == 0.7:
                    eff = eff*line[9]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[9]+(line[13]-line[9])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[13]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[13]+(line[17]-line[13])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[17]
                    EFF.append(eff)

        else:
            eff = eff*line[18]

            if I <= 9:

                if seeing == 0.7:
                    eff = eff*line[19]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[19]+(line[23]-line[19])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[23]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[23]+(line[27]-line[23])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[27]
                    EFF.append(eff)

            elif I > 9 and I < 10:

                if seeing == 0.7:
                    f = line[19]+(line[20]-line[19])*(I-9.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f_07 = line[19]+(line[20]-line[19])*(I-9.)/1.
                    f_09 = line[23]+(line[24]-line[23])*(I-9.)/1.
                    f = f_07+(f_09-f_07)*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    f = line[23]+(line[24]-line[23])*(I-9.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f_09 = line[23]+(line[24]-line[23])*(I-9.)/1.
                    f_12 = line[27]+(line[28]-line[27])*(I-9.)/1.
                    f = f_09+(f_12-f_09)*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    f = line[27]+(line[28]-line[27])*(I-9.)/1.
                    eff = eff*f
                    EFF.append(eff)

            elif I == 10:

                if seeing == 0.7:
                    eff = eff*line[20]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[20]+(line[24]-line[20])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[24]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[24]+(line[28]-line[24])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[28]
                    EFF.append(eff)

            elif I > 10 and I < 11:

                if seeing == 0.7:
                    f = line[20]+(line[21]-line[20])*(I-10.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f_07 = line[20]+(line[21]-line[20])*(I-10.)/1.
                    f_09 = line[24]+(line[25]-line[24])*(I-10.)/1.
                    f = f_07+(f_09-f_07)*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    f = line[24]+(line[25]-line[24])*(I-10.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f_09 = line[24]+(line[25]-line[24])*(I-10.)/1.
                    f_12 = line[28]+(line[29]-line[28])*(I-10.)/1.
                    f = f_09+(f_12-f_09)*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    f = line[28]+(line[29]-line[28])*(I-10.)/1.
                    eff = eff*f
                    EFF.append(eff)

            elif I == 11:

                if seeing == 0.7:
                    eff = eff*line[21]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[21]+(line[25]-line[21])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[25]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[25]+(line[29]-line[25])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[29]
                    EFF.append(eff)

            elif I > 11 and I < 12:

                if seeing == 0.7:
                    f = line[21]+(line[22]-line[21])*(I-11.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f_07 = line[21]+(line[22]-line[21])*(I-11.)/1.
                    f_09 = line[25]+(line[26]-line[25])*(I-11.)/1.
                    f = f_07+(f_09-f_07)*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    f = line[25]+(line[26]-line[25])*(I-11.)/1.
                    eff = eff*f
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f_09 = line[25]+(line[26]-line[25])*(I-11.)/1.
                    f_12 = line[29]+(line[30]-line[29])*(I-11.)/1.
                    f = f_09+(f_12-f_09)*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    f = line[29]+(line[30]-line[29])*(I-11.)/1.
                    eff = eff*f
                    EFF.append(eff)

            elif I >= 12:

                if seeing == 0.7:
                    eff = eff*line[22]
                    EFF.append(eff)

                elif seeing > 0.7 and seeing < 0.9:
                    f = line[22]+(line[26]-line[22])*(seeing-0.7)/0.2
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 0.9:
                    eff = eff*line[26]
                    EFF.append(eff)

                elif seeing > 0.9 and seeing < 1.2:
                    f = line[26]+(line[30]-line[26])*(seeing-0.9)/0.3
                    eff = eff*f
                    EFF.append(eff)

                elif seeing == 1.2:
                    eff = eff*line[30]
                    EFF.append(eff)

    total_effs = EFF

    return (total_effs, wavelengths, wavelengths_nm, orders, central_wave,
            beg_wave, end_wave, order_wave, obs_mode)


def calc_readout_noise(t_exp):
    # N_read=82.16*t_exp**(-1)  # Readout noise in sqrt(electrons/pixel) (Email from E. Artigau)
    N_read = np.sqrt(400./(t_exp/5.57 + 1)+36.)  # Readout noise in sqrt(electrons/pixel)
    # N_read=0.1  # only for caluclation of SNR
    return N_read


def calc_instrument_parameters(d_t, seeing, wavelength, resolving_power):
    # A=numpy.pi*(d_t/2.)**2  # Calculate the aperture (assuming no central_obscuration)
    # A from: https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html
    A = 8.8564  # M1 clear area (mˆ2)
    seeing_rads = seeing/3600.*np.pi/180.  # Convert the seeing arcsecs into radians
    # resolution = wavelength/resolving_power  # Spectral resolution in meters
    delta_lambda = wavelength/resolving_power  # Spectral resolution in meters #NG change 1feb2022 as this is delta lambda not resolution
    P = h*c/wavelength                             # Energy of one photon in Joules (J/ph)
    return (A, seeing_rads, delta_lambda, P)


def calc_electrons_rate(z, I, delta_lambda, sampling, A, P, t_exp, total_eff,
                        flux_st, st, H, Ho):

    # Flux of object from IRTF Library template in (W m-2 um-1) scaled by magnitude
    flux=(10.**(0.4*(Ho-H)))*flux_st 

    # Convert bin size meters to microns
    bin_size_pxl=delta_lambda*1e6/sampling  # /sampling to get per pixel
    bin_size_bin=delta_lambda*1e6

    # Calculate electrons rate for the object (e-)
    N_obj_pxl=flux*bin_size_pxl*t_exp*total_eff*A/P
    N_obj_bin=flux*bin_size_bin*t_exp*total_eff*A/P     # Calculate electrons rate for the object (e-)

    return (N_obj_pxl, N_obj_bin, bin_size_pxl, bin_size_bin, flux)


def calc_signal_to_noise_ratio(N_obj_pxl, N_obj_bin, n_dark, n_bin_y, N_read, t_exp):
    N_star_pxl = N_obj_pxl
    N_star_bin = N_obj_bin
    N_sky = 0  # self.N_sky = self.calc_N_obj(self.V_sky, self.z)  Negligeble for IJH bands (not for K band)
    S_N_pxl = N_star_pxl/np.sqrt(N_star_pxl+n_bin_y*N_sky+n_bin_y*N_read**2.+n_bin_y*n_dark*t_exp/3600.)
    S_N_bin = N_star_bin/np.sqrt(N_star_bin+n_bin_y*N_sky+n_bin_y*N_read**2.+n_bin_y*n_dark*t_exp/3600.)
    return (S_N_pxl, S_N_bin)


def calc_S_N_order(obs_mode, t_exp, d_t, seeing, wavelength, resolving_power,
                   z, I, total_eff, flux_st, N_dark, st, H, Ho):
    # WARN: Fiber diameter unused
    fiber_diameter = calc_fiber_diameter(obs_mode)
    sampling = calc_sampling(obs_mode)
    N_read = calc_readout_noise(t_exp)

    ### Calculating Instrument Parameters ###
    # change resolution to delta_lambda 1feb2022
    (A, seeing_rads, delta_lambda, P) = calc_instrument_parameters(
        d_t, seeing, wavelength, resolving_power
    )

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
    # change resolution to delta_lambda 1feb2022
    (N_obj_pxl, N_obj_bin, bin_size_pxl, bin_size_bin, flux) = calc_electrons_rate(
        z, I, delta_lambda, sampling, A, P, t_exp, total_eff, flux_st, st, H, Ho
    )

    ### Calculating te Signal-to-Noise Rate ###
    # Formulae from ESO exposure calculator
    #
    n_bin_y = calc_N_pix_Y(obs_mode)                          # Bins in y direction

    (S_N_pxl, S_N_bin) = calc_signal_to_noise_ratio(
        N_obj_pxl, N_obj_bin, N_dark, n_bin_y, N_read, t_exp
    )

    return (N_obj_pxl, N_obj_bin, bin_size_pxl, bin_size_bin, flux, S_N_pxl, S_N_bin)


def calc_rv_precision(st, obs_mode, sn_h, spirou_fit_qvalues_file,
                      phoenix_Q_conversions_file, phoenix_eniric_Qfactors_file,
                      N_OBJ_arr, wavelengths_nm, bandpass,
                      wlnmax_y, wlnmax_j, wlnmax_h,
                      wlnmin_y, wlnmin_j, wlnmin_h):
    # temps = np.array([5000,4500,4000,3900,3700,3600,3400,3200,3000,2800,2700,2600,2500]) #original Teff values available from Spirou Template interpolation
    # stypes = np.array(['K2V','K5V','K7V','M0V','M1V','M2V','M3V','M4V','M5V','M6V','M7V','M8V','M9V'])
    temps = np.array(
        [5000, 4000, 3900, 3700, 3600, 3400, 3200, 3000, 2800, 2700, 2600, 2500]
    )
    stypes = np.array(
        ['K3V', 'K7V', 'M0V', 'M1V', 'M2V', 'M3V', 'M4V', 'M5V', 'M6V', 'M7V', 'M8V', 'M9V']
    )
    spirou_qs = spirou_results = pd.read_csv(
        'inputs/'+spirou_fit_qvalues_file, sep=',', header=0
    )
    q_resolution_conversion = pd.read_csv(
        'inputs/'+phoenix_Q_conversions_file, sep=',', header=0
    )
    eniric_pheonix_results = pd.read_csv(
        'inputs/'+phoenix_eniric_Qfactors_file, sep=',', header=0
    )

    index = [np.where(stypes == st)[0]][0]
    temp = temps[index][0]
    print('')
    print(' RADIAL VELOCITY PRECISION:')
    print('')
    print(st, 'Teff =', temp)
    print('YJH total e-: %i %i %i' % (
        sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))]),
        sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_j))]),
        sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_h))])
    ))
    print('')
    print('-- RV precision from Spirou Template Qs --')
    print('')

    spirou_quality_y = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_Y'].values[0]
    spirou_quality_j = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_J'].values[0]
    spirou_quality_h = spirou_qs[(spirou_qs['Teff'] == temp)]['Q_H'].values[0]
    if obs_mode == 'HE':
        resolution = '70k'
        sampling = 4.3
    if obs_mode == 'HA':
        resolution = '80k'
        sampling = 3.8
        q_conv_y = q_resolution_conversion[(spirou_qs['Teff'] == temp)]['Q_Y'].values[0]
        q_conv_j = q_resolution_conversion[(spirou_qs['Teff'] == temp)]['Q_J'].values[0]
        q_conv_h = q_resolution_conversion[(spirou_qs['Teff'] == temp)]['Q_H'].values[0]
        spirou_quality_y = spirou_quality_y / q_conv_y
        spirou_quality_j = spirou_quality_j / q_conv_j
        spirou_quality_h = spirou_quality_h / q_conv_h

    rvprec_y = c/(spirou_quality_y*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))])))
    rvprec_j = c/(spirou_quality_j*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_j) & (wavelengths_nm <= wlnmax_j))])))
    rvprec_h = c/(spirou_quality_h*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_h) & (wavelengths_nm <= wlnmax_h))])))
    combine_yjh_sp = 1/np.sqrt(rvprec_y**-2+rvprec_j**-2+rvprec_h**-2)
    print('YJH Quality factors: %.1f %.1f %.1f' % (
        spirou_quality_y, spirou_quality_j, spirou_quality_h
    ))
    print('YJH RV precision (m/s): %.2f %.2f %.2f' % (rvprec_y, rvprec_j, rvprec_h))
    print('YJH combined RV precision (m/s): %.2f' % combine_yjh_sp)

    def eniric_rv_precision(vsini):
        print('')
        print('vsini = %.1f km/s' % vsini)
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
        print('YJH Quality factors:', ph_quality_y, ph_quality_j, ph_quality_h)

        rvprec_y = c/(ph_quality_y*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))])))
        rvprec_j = c/(ph_quality_j*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_j) & (wavelengths_nm <= wlnmax_j))])))
        rvprec_h = c/(ph_quality_h*np.sqrt(sum(N_OBJ_arr[((wavelengths_nm >= wlnmin_h) & (wavelengths_nm <= wlnmax_h))])))
        combine_yjh_en = 1/np.sqrt(rvprec_y**-2+rvprec_j**-2+rvprec_h**-2)

        print('YJH RV precision (m/s): %.2f %.2f %.2f' % (rvprec_y, rvprec_j, rvprec_h))
        print('YJH combined RV precision (m/s): %.2f' % combine_yjh_en)
        return(combine_yjh_en)

    print('')
    print('-- RV precision from Eniric Pheonix spectra Qs --')
    vsinis = [0.1, 1.0, 5.0, 10.0]
    eniric_rv_vsinis = np.zeros(len(vsinis))
    for ivsini in range(len(vsinis)):
        eniric_rv_vsinis[ivsini] = eniric_rv_precision(vsinis[ivsini])

    return(combine_yjh_sp, eniric_rv_vsinis)
