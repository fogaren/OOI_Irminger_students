%% Master script for Lucy Wanzer's thesis
%Irminger005_CTD loads in all the CTD data from the Irminger 005 cruise and
%makes depth plots as well as a map that shows where the CTD casts were
%taken
    Irminger005_CTD
    close all

%% Plotting the raw OOI data
%load_HYPM loads in all the OOI data (which should already be
%on your path), grids the data onto depth and temperature grids, and plots 
%the uncalibrated data
    load_HYPM
    close all
    
%% %Correcting the oxygen data
%wfp_Irminger_winklercalibration_allyears loads in the discrete oxygen data from previous cruises, 
%compares this data to the profiler mooring data, and calculates the initial gain correction and makes plots
%illustrating the initial gain correction
    wfp_Irminger_winklercalibration_allyrs
    
%HYPM_O2_gaincorr plots the gain corrected data
    HYPM_O2_gaincorr
    
%wfp_deepIsotherm_driftCorrection calculates oxygen drift overtime and
%creates plots that visualize the oxygen drift overtime
    wfp_deepIsotherm_driftCorrection
    
%plot_driftcorr takes the calculations made in the script above and applies
%them to the data, and then produces final plots with fully calibrated data
    plot_driftcorr
    
%% Caclaulate and plot respiration rates for each year
%IrmingerRespiration smoothes the oxygen data, identifies the beginning and
%end of stratification for each season, and makes plots visualizing
%respiraiton in the water column
    IrmingerRespiration
    
%O2_Chl_w_redlines puts red lines on the oxygen, chl, and backscatter time
%series to visualize the beginning and end of the stratified season
    O2_Chl_w_redlines
