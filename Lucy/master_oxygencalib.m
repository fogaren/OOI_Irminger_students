%Lucy's master thesis script 
%first run HYPM_1 to load in all data and create plots without calibrations
    HYPM_1
    close all
%Next run wfp_Irminger_winklercalibration_allyears which makes initial gain
%correction
%%
    wfp_Irminger_winklercalibration_allyrs
%Next run HYPM_2 which makes plots with initial gain calculation and grids
%onto isotherms
    HYPM_2
%Next run wfp_deepIsotherm_driftCorrection to correct for oxygen drift overtime
    wfp_deepIsotherm_driftCorrection
