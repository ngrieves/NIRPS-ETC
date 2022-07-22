import pandas as pd
import numpy as np

from nirps_etc_driver import run_nirps_etc

input_targets_file = 'etc_targets_input.txt'
output_targets_file = 'etc_targets_output.txt'

input_targets = pd.read_csv(input_targets_file, sep=r"\s+", header=0)
targets = input_targets['target']
sts = input_targets['st']
obs_modes = input_targets['obs_mode']
seeings = input_targets['seeing']
airmasses = input_targets['airmass']
Hmags = input_targets['H']
t_exps = input_targets['t_exp']
bandpasses = input_targets['bandpass']
waveselect = 'FSR'

output_targets = open(output_targets_file, "w")
# output_targets.write('target Mean_S/N_(ph/pxl) Mean_S/N_Y(ph/pxl) Mean_S/N_J(ph/pxl) Mean_S/N_H(ph/pxl) Hmax_1619nm_S/N(ph/pxl) EchelleOrd90_S/N(ph/pxl) spRV enRV_vsini01 enRV_vsini1 enRV_vsini5 enRV_vsini10 \n')
output_targets.write(
    'target Mean_S/N_(ph/pxl) Mean_S/N_Y(ph/pxl) Mean_S/N_J(ph/pxl) Mean_S/N_H(ph/pxl) EchelleOrd90_S/N(ph/pxl) spRV enRV_vsini01 enRV_vsini1 enRV_vsini5 enRV_vsini10 \n'
)

### Telescope parameters ###
# airmass = 1.0                     # Airmass

for itarget in range(len(targets)):
    target = targets[itarget]
    obs_mode = obs_modes[itarget]
    seeing = seeings[itarget]
    airmass = airmasses[itarget]
    H = Hmags[itarget]
    t_exp = t_exps[itarget]
    st = sts[itarget]
    bandpass = bandpasses[itarget]

    (SN_pxl_order_central, SN_pxl_Y, SN_pxl_J, SN_pxl_H, order_wave,
     rv_total_spirou, rv_total_eniric) = run_nirps_etc(obs_mode, st, H, seeing,
                                                       airmass, t_exp,
                                                       bandpass,
                                                       waveselect, plot=True,
                                                       name=target, show=False,
                                                       save_details=False,
                                                       save_plot=False)

    if isinstance(rv_total_spirou, str):
        spRV = rv_total_spirou
        enRV_vsini01 = rv_total_eniric[0]
        enRV_vsini1 = rv_total_eniric[1]
        enRV_vsini5 = rv_total_eniric[2]
        enRV_vsini10 = rv_total_eniric[3]
    else:
        spRV = '{:.2f}'.format(rv_total_spirou)
        enRV_vsini01 = '{:.2f}'.format(float(rv_total_eniric[0]))
        enRV_vsini1 = '{:.2f}'.format(float(rv_total_eniric[1]))
        enRV_vsini5 = '{:.2f}'.format(float(rv_total_eniric[2]))
        enRV_vsini10 = '{:.2f}'.format(float(rv_total_eniric[3]))

    print ("=================================================================\n\n")
    output_sn = target+' %.1f %.1f %.1f %.1f %.1f ' % (
        np.nanmean(SN_pxl_order_central[(SN_pxl_order_central > 0)]),
        SN_pxl_Y, SN_pxl_J, SN_pxl_H,
        SN_pxl_order_central[order_wave == 90][0])
    output_rv = spRV+' '+enRV_vsini01+' '+enRV_vsini1+' '+enRV_vsini5+' '+enRV_vsini10
    output_line = output_sn + output_rv + '\n'
    output_targets.write(output_line)

output_targets.close()
