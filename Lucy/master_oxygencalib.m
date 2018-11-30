%Lucy's master script for oxygen calibration 
%first run HYPM_1 to load in all data and create plots without calibrations
    load_HYPM
    close all
%Next run wfp_Irminger_winklercalibration_allyears which makes initial gain
%correction
%%
    wfp_Irminger_winklercalibration_allyrs
%Next run HYPM_2 which makes plots with initial gain calculation and grids
%onto isotherms
    HYPM_O2_gaincorr
%Next run wfp_deepIsotherm_driftCorrection to correct for oxygen drift overtime
    wfp_deepIsotherm_driftCorrection
%Make plots with drift correction
    plot_driftcorr
