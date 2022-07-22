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
# Get SNR of orders from the central 5% (June 2022)

from os import makedirs

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import nirps_etc_calc as nec


def run_nirps_etc(obs_mode, st, H, seeing, airmass, t_exp, bandpass,
                  waveselect, plot=True, name=None, show=True,
                  save_details=True, save_plot=True):
    # Telescope diameter (m) [https://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview.html]
    # Area is input directly as of July 2022, so d_t not used.
    d_t = 3.57

    pixel_size = 15e-6                 # pixel size (m)
    N_dark = 0                         # Lab tests in Dec/2017 (Email from E. Artigau)


    name_suffix = f"_{name}" if name is not None else ""

    # zero point of the photometric system taken from:
    # http://www.eso.org/observing/etc/doc/formulabook/node12.html
    # (ergs/s/cmˆ2/A)
    # z = 1.22603e-09  # I band
    z = 1.14e-10  # H  band

    if bandpass == 'CFHT':
        spirou_fit_qvalues_file = 'spirou_fit_Qvalues_CFHT-bandpass.txt'
        phoenix_Q_conversions_file = 'phoenix_Q_conversions_CFHT-bandpass.txt'
        phoenix_eniric_Qfactors_file = 'phoenix_eniric_Qfactors_CFHT-bandpass.csv'

    # use for Eniric bandpass selected below
    if bandpass == 'Eniric':
        spirou_fit_qvalues_file = 'spirou_fit_Qvalues_eniric-bandpass.txt'
        phoenix_Q_conversions_file = 'phoenix_Q_conversions_eniric-bandpass.txt'
        phoenix_eniric_Qfactors_file = 'phoenix_eniric_Qfactors_eniric-bandpass.csv'

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

    resolving_power = nec.calc_resoving_power(obs_mode)
    I_max = nec.calc_I_max(obs_mode)

    (flux_sts, I, Ho) = nec.calc_flux_spectral_type(st, H, st_templates_file)

    (total_effs, wavelengths, wavelengths_nm, orders, central_wave, beg_wave,
     end_wave, order_wave, obs_mode) = nec.calc_total_efficiency(
        obs_mode, seeing, airmass, I, wave_range_file, tapas_file, effs_file,
    )

    S_N_pxl = []
    S_N_bin = []
    FLUX = []
    N_OBJ = []
    BIN_SIZE = []
    S_N_pxl_mean = []
    S_N_bin_mean = []
    N_OBJ_mean = []
    BIN_SIZE_mean = []
    FLUX_mean = []
    EFF_mean = []
    SATURATION = []
    N_OBJ_PIXEL = []  # NG add in 13 July 2022 to account for juste-/pxl not used in RV caluclation

    i = 0
    SNR_H_pxl = 0
    SNR_H_bin = 0
    order_now = orders[0]
    s_n_pxl_mean = 0
    s_n_bin_mean = 0
    n_obj_mean = 0
    bin_size_mean = 0
    eff_mean = 0
    flux_mean = 0
    divisor = 0
    cont_saturation = 0.0

    while i < len(wavelengths):
        wavelength = wavelengths[i]
        wavelength_nm = wavelengths_nm[i]
        order = orders[i]
        total_eff = total_effs[i]
        flux_st = flux_sts[i]

        if i == len(wavelengths)-1:
            delta_wavelength = (wavelengths_nm[-1]-wavelengths_nm[-2])*0.001  # nm to microns
        else:
            delta_wavelength = (wavelengths_nm[i+1]-wavelengths_nm[i])*0.001

        if order_now == order:
            (n_obj_pxl, n_obj_bin, bin_size_pxl, bin_size_bin, flux, s_n_pxl,
             s_n_bin) = nec.calc_S_N_order(obs_mode, t_exp, d_t, seeing,
                                           wavelength, resolving_power, z, I,
                                           total_eff, flux_st, N_dark, st, H,
                                           Ho)
            S_N_pxl.append(s_n_pxl)
            S_N_bin.append(s_n_bin)
            FLUX.append(flux)
            # N_OBJ.append(n_obj_pxl)
            # For total number of electrons: multiply by ((delta wavlength of the array) / (bin size))
            # to account for the array size of the effs/tapas/stellar template wavelength values
            N_OBJ.append(n_obj_pxl*(delta_wavelength/bin_size_pxl))
            N_OBJ_PIXEL.append(n_obj_pxl)  # NG add in 13 July 2022 to account for juste-/pxl not used in RV caluclation
            if n_obj_pxl > I_max:
                cont_saturation += 1.0

            BIN_SIZE.append(bin_size_pxl)
            bin_size_mean = bin_size_mean+bin_size_pxl
            n_obj_mean = n_obj_mean+n_obj_pxl
            s_n_pxl_mean = s_n_pxl_mean+s_n_pxl
            s_n_bin_mean = s_n_bin_mean+s_n_bin
            eff_mean = eff_mean+total_eff
            flux_mean = flux_mean+flux
            order_now = order
            divisor = divisor+1
        else:
            sn_pxl = s_n_pxl_mean/divisor
            sn_bin = s_n_bin_mean/divisor
            nobj = n_obj_mean/divisor
            binsize = bin_size_mean/divisor
            effall = eff_mean/divisor
            fluxobj = flux_mean/divisor
            saturation = cont_saturation/divisor
            S_N_pxl_mean.append(sn_pxl)
            S_N_bin_mean.append(sn_bin)
            N_OBJ_mean.append(nobj)
            BIN_SIZE_mean.append(binsize)
            EFF_mean.append(effall)
            FLUX_mean.append(fluxobj)
            SATURATION.append(saturation)
            s_n_pxl_mean = 0
            s_n_bin_mean = 0
            n_obj_mean = 0
            bin_size_mean = 0
            eff_mean = 0
            flux_mean = 0
            divisor = 0
            cont_saturation = 0.0
            (n_obj_pxl, n_obj_bin, bin_size_pxl, bin_size_bin, flux, s_n_pxl,
             s_n_bin) = nec.calc_S_N_order(obs_mode, t_exp, d_t, seeing,
                                           wavelength, resolving_power, z, I,
                                           total_eff, flux_st, N_dark, st, H,
                                           Ho)
            S_N_pxl.append(s_n_pxl)
            S_N_bin.append(s_n_bin)
            FLUX.append(flux)
            N_OBJ.append(n_obj_pxl*(delta_wavelength/bin_size_pxl))  # NG fix 13 July 2022
            N_OBJ_PIXEL.append(n_obj_pxl)  # NG add in 13 July 2022 to account for juste-/pxl not used in RV caluclation

            if n_obj_pxl > I_max:
                cont_saturation += 1.0

            BIN_SIZE.append(bin_size_pxl)
            s_n_pxl_mean = s_n_pxl_mean+s_n_pxl
            s_n_bin_mean = s_n_bin_mean+s_n_bin
            n_obj_mean = n_obj_mean+n_obj_pxl
            bin_size_mean = bin_size_mean+bin_size_pxl
            eff_mean = eff_mean+total_eff
            flux_mean = flux_mean+flux
            order_now = order
            divisor = divisor+1

        i = i+1

    # Hband SNR
    hband_wave = 1619
    hband_ind = np.argmin(np.abs(np.array(wavelengths_nm) - hband_wave))
    SNR_pxl_H = S_N_pxl[hband_ind]
    SNR_bin_H = S_N_bin[hband_ind]
    sn_h = SNR_bin_H

    sn_pxl = s_n_pxl_mean/divisor
    sn_bin = s_n_bin_mean/divisor
    nobj = n_obj_mean/divisor
    binsize = bin_size_mean/divisor
    effall = eff_mean/divisor
    fluxobj = flux_mean/divisor

    saturation = cont_saturation/divisor
    S_N_pxl_mean.append(sn_pxl)
    S_N_bin_mean.append(sn_bin)
    N_OBJ_mean.append(nobj)
    BIN_SIZE_mean.append(binsize)
    EFF_mean.append(effall)
    FLUX_mean.append(fluxobj)
    SATURATION.append(saturation)

    print("=================================================================")

    if bandpass == 'Eniric':
        print('### Eniric bandpasses ###')
        # Y-band
        wlnmin_y = 1000.0
        wlnmax_y = 1100.0
        # J-band
        wlnmin_j = 1170.0
        wlnmax_j = 1330.0
        # H-band
        wlnmin_h = 1500.0
        wlnmax_h = 1750.0
        # K-band
        # wlnmin_k = 2070.0
        # wlnmax_k = 2350.0

    if bandpass == 'CFHT':
        print('#### CFHT/WIRCAM bandpasses ###')
        # Y-band
        wlnmin_y = 938.6
        wlnmax_y = 1113.4
        # J-band
        wlnmin_j = 1153.5
        wlnmax_j = 1354.4
        # H-band
        wlnmin_h = 1462.8
        wlnmax_h = 1808.5
        # K-band
        # wlnmin_k = 1957.7
        # wlnmax_k = 2343.1

    ### OLD NIRPS ETC LIMITS ###
    # wlnmin_y = 980
    # wlnmax_y = 1110
    # wlnmin_j = 1200
    # wlnmax_j = 1330
    # wlnmin_h = 1510
    # wlnmax_h = 1735

    print('Y band: '+str(wlnmin_y)+' - '+str(wlnmax_y)+' nm')
    print('J band: '+str(wlnmin_j)+' - '+str(wlnmax_j)+' nm')
    print('H band: '+str(wlnmin_h)+' - '+str(wlnmax_h)+' nm')

    print("\nTEMPLATE STAR: %s" % (st))
    print("H (mag): %.3f" % (H))
    print("Exposure time (seconds): %.1f\n" % (t_exp))

    print("Observation Mode: %s" % (obs_mode))
    print("seeing (arcsec): %.1f" % (seeing))
    print("airmass: %.1f\n" % (airmass))

    print("Saturation limit (e-/pxl): %d\n\n" % (I_max))

    print("=================================================================")

    if save_details or save_plot:
        makedirs("outputs", exist_ok=True)
    if save_details:
        text_file = open(f"outputs/order_snrs{name_suffix}.txt", "w")
        text_file.write(
            "##--  calculated from the central 5% of each order --## \n"
        )
        text_file.write(
            "order central_wave (beg-end)    Eff.    object     snr      snr         sat \n"
        )
        text_file.write(
            "          (um)                         (e-/pxl)  (ph/pxl) (ph/res elem) (%) \n"
        )

    ######## UPDATE 13 June 2022 ########
    # get central 200 pixels (central 5%) for oder S/N instead of all wavelengths
    #####################################
    SN_pxl_order_central = np.zeros(len(order_wave))
    SN_bin_order_central = np.zeros(len(order_wave))
    N_OBJ_order_central = np.zeros(len(order_wave))
    EFF_order_central = np.zeros(len(order_wave))
    # len(S_N_pxl)
    order_sampling_size = int(len(S_N_pxl)/len(order_wave))
    ind1_center = int(len(S_N_pxl)/len(order_wave)*0.475)-1
    ind2_center = int(len(S_N_pxl)/len(order_wave)*0.525)-1

    for i in range(len(order_wave)):
        SN_pxl_order_central[i] = np.nanmean(
            S_N_pxl[(i*order_sampling_size+ind1_center):(i*order_sampling_size+ind2_center)]
        )
        SN_bin_order_central[i] = np.nanmean(
            S_N_bin[(i*order_sampling_size+ind1_center):(i*order_sampling_size+ind2_center)]
        )  # [(i*order_sampling_size+ind1_center):(i*order_sampling_size+ind2_center)])
        N_OBJ_order_central[i] = np.nanmean(
            N_OBJ_PIXEL[(i*order_sampling_size+ind1_center):(i*order_sampling_size+ind2_center)]
        )
        EFF_order_central[i] = np.nanmean(
            total_effs[(i*order_sampling_size+ind1_center):(i*order_sampling_size+ind2_center)]
        )
        if save_details:
            text_file.write(
                '# %s %7.2f (%7.2f-%7.2f) ' % (
                    repr(order_wave[i]).rjust(3), central_wave[i], beg_wave[i], end_wave[i]
                )
            )
            text_file.write(
                '%5.3f %.5e %5.1f %8.1f %s \n' % (
                    EFF_order_central[i], N_OBJ_order_central[i], SN_pxl_order_central[i],
                    SN_bin_order_central[i], repr(int(SATURATION[i]*100.)).rjust(10)
                )
            )

    central_wave = np.array(central_wave)
    order_wave = np.array(order_wave)
    wavelengths_nm = np.array(wavelengths_nm)
    total_effs = np.array(total_effs)
    N_OBJ_arr = np.array(N_OBJ)

    SN_pxl_Y = np.nanmean(
        SN_pxl_order_central[((central_wave >= wlnmin_y) & (central_wave <= wlnmax_y))]
    )
    SN_pxl_J = np.nanmean(
        SN_pxl_order_central[((central_wave >= wlnmin_j) & (central_wave <= wlnmax_j))]
    )
    SN_pxl_H = np.nanmean(
        SN_pxl_order_central[((central_wave >= wlnmin_h) & (central_wave <= wlnmax_h))]
    )

    SN_bin_Y = np.nanmean(
        SN_bin_order_central[((central_wave >= wlnmin_y) & (central_wave <= wlnmax_y))]
    )
    SN_bin_J = np.nanmean(
        SN_bin_order_central[((central_wave >= wlnmin_j) & (central_wave <= wlnmax_j))]
    )
    SN_bin_H = np.nanmean(
        SN_bin_order_central[((central_wave >= wlnmin_h) & (central_wave <= wlnmax_h))]
    )

    mean_eff_Y = np.nanmean(
        total_effs[((wavelengths_nm >= wlnmin_y) & (wavelengths_nm <= wlnmax_y))]
    )
    mean_eff_J = np.nanmean(
        total_effs[((wavelengths_nm >= wlnmin_j) & (wavelengths_nm <= wlnmax_j))]
    )
    mean_eff_H = np.nanmean(
        total_effs[((wavelengths_nm >= wlnmin_h) & (wavelengths_nm <= wlnmax_h))]
    )

    if save_details:
        text_file.close()
    print("-----------------------------------------------------------------")
    print("\n SIGNAL TO NOISE RATIO:\n")
    print("Mean S/N of center of all orders: %5.1f (ph/pxl) | %5.1f (ph/res elem)" % (
        np.nanmean(SN_pxl_order_central[(SN_pxl_order_central > 0)]),
        np.nanmean(SN_bin_order_central[SN_bin_order_central > 0])
    ))
    print("Mean S/N of center of orders (ph/pxl): Y=%5.1f | J=%5.1f | H=%5.1f" % (
        SN_pxl_Y, SN_pxl_J, SN_pxl_H
    ))
    print("         (ph/res elem): Y=%5.1f | J=%5.1f | H=%5.1f\n\n" % (
        SN_bin_Y, SN_bin_J, SN_bin_H
    ))
    print(
        'S/N in H ('+str(hband_wave)+' nm): %5.1f (ph/pxl) | %5.1f (ph/res elem)\n' % (
            SNR_pxl_H, SNR_bin_H
        )
    )
    print("-----------------------------------------------------------------")
    print("\nMean Efficiency: %5.3f " % (np.nanmean(total_effs)))
    print("Mean Efficiencies Y=%5.3f | J=%5.3f | H=%5.3f \n" % (
        mean_eff_Y, mean_eff_J, mean_eff_H
    ))
    print("-----------------------------------------------------------------")

    N_OBJ_arr = np.array(N_OBJ)
    wavelengths_nm = np.array(wavelengths_nm)

    stypes_rv = np.array(['K3V', 'K7V', 'M0V', 'M1V', 'M2V', 'M3V', 'M4V', 'M5V', 'M6V', 'M7V', 'M8V', 'M9V'])
    if len(np.where(stypes_rv == st)[0]) == 1:
        rv_total_spirou, rv_total_eniric = nec.calc_rv_precision(
            st, obs_mode, sn_h, spirou_fit_qvalues_file,
            phoenix_Q_conversions_file, phoenix_eniric_Qfactors_file,
            N_OBJ_arr, wavelengths_nm, bandpass,
            wlnmax_y, wlnmax_j, wlnmax_h,
            wlnmin_y, wlnmin_j, wlnmin_h)
    else:
        print("# It's not possible to calculate the RV precision of this spectral type with the NIRPS ETC currently #\n")
        rv_total_spirou = "--"
        rv_total_eniric = [rv_total_spirou] * 4

    if plot:
        ################################################
        #
        #
        #        Plot SNR vs Wavelength
        #
        #
        ################################################

        fig = plt.figure(figsize=(15,  8))
        plt.plot(wavelengths_nm, S_N_pxl, label="All wavelength", color='gainsboro', zorder=1)
        plt.axvspan(wlnmin_y, wlnmax_y,  facecolor='#2ca02c',  alpha=0.1)
        plt.text(1045, 0, 'Y', fontsize=18)
        plt.axvspan(wlnmin_j, wlnmax_j,  facecolor='#2ca02c',  alpha=0.1)
        plt.text(1265, 0, 'J', fontsize=18)
        plt.axvspan(wlnmin_h, wlnmax_h,  facecolor='#2ca02c',  alpha=0.1)
        plt.text(1622.5, 0, 'H', fontsize=18)
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        cmap = plt.get_cmap('jet')
        # for (cw, snr, sat) in zip(central_wave, S_N_pxl_mean, SATURATION):
        #     plt.scatter(cw, snr, color=cmap(sat), zorder=2)
        for (cw, snr, sat) in zip(central_wave, SN_pxl_order_central, SATURATION):
            plt.scatter(cw, snr, color=cmap(sat), zorder=2)
        plt.title('S/N vs Wavelength graph',  fontsize=18)
        plt.ylabel('S/N (ph/pxl)',  fontsize=16)
        plt.xlabel('Wavelength (nm)',  fontsize=16)
        ax1 = fig.add_axes([0.91,  0.1,  0.01,  0.8])
        norm = mpl.colors.Normalize(vmin=0,  vmax=100)
        cb1 = mpl.colorbar.ColorbarBase(ax1,  cmap=cmap,
                                        norm=norm,
                                        orientation='vertical')
        cb1.set_label('Saturation (%)', fontsize=16)
        plt.xlim = (900, 1900)
        plt.plot(wavelengths_nm, S_N_pxl, label="All wavelength", color='gainsboro', zorder=1)
        if save_plot:
            plt.savefig(f"outputs/snr_wln{name_suffix}.pdf")
        if show:
            plt.show()
        else:
            plt.close()

    print("=================================================================\n\n")

    return SN_pxl_order_central, SN_pxl_Y, SN_pxl_J, SN_pxl_H, order_wave, rv_total_spirou, rv_total_eniric
