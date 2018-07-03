  %% File Loading
        %Apex Near Surface Mooring 
        NearSurface = ['deployment0003_GI01SUMO-RID16-06-DOSTAD000-recovered_host-dosta_abcdjm_dcl_instrument_recovered_20160710T174508.675000-20170812T204759.773000.nc']; 
        
        %Apex Surface Mooring 
        Surface = ['deployment0003_GI01SUMO-SBD11-04-DOSTAD000-recovered_host-dosta_abcdjm_dcl_instrument_recovered_20160710T174509.583000-20170812T200026.036000 (1).nc'];
        
        %Flanking A
        FlankingA = ['deployment0003_GI03FLMA-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20160712T160001-20160830T060001.nc'];
        
        %Flanking B
        FlankingB = ['deployment0003_GI03FLMB-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20160713T160001-20170112T181502.nc']; 

%% Variable Naming 
    %Time
    ApexS.Time =convertTime(ncread(Surface,'time'));
    NearS.Time = convertTime(ncread(NearSurface,'time'));
    FlankingA_strt.Time = convertTime(ncread(FlankingA, 'time'));
    FlankingB_strt.time = convertTime(ncread(FlankingB, 'time'));
    %%
    %Oxygen
   % Check if "dissolved oxygen" is right
   %'Dosta' appears to be a non-adjusted reading, it consistently reads
   %extremely high. 
    ApexS.O2_dosta = ncread(Surface,'dosta_abcdjm_cspp_tc_oxygen');
    NearS.O2_dosta = ncread(NearSurface,'dosta_abcdjm_cspp_tc_oxygen');
    FlankingA_strt.O2_dosta = ncread(FlankingA, 'dosta_abcdjm_cspp_tc_oxygen');
    FlankingB_strt.O2_dosta = ncread(FlankingB, 'dosta_abcdjm_cspp_tc_oxygen');
    %Here I am naming O_2, I think that O_2 is correct 
    ApexS.O2 = ncread(Surface,'dissolved_oxygen');
    NearS.O2 = ncread(NearSurface,'dissolved_oxygen');
    FlankingA_strt.O2 = ncread(FlankingA, 'dissolved_oxygen');
    FlankingB_strt.O2 = ncread(FlankingB, 'dissolved_oxygen');
   %% 
    %Temperature (C °)
    %temp_ApexS = ncread(Surface,'temp_compensated_phase');
    NearS.temp = ncread(NearSurface,'temp');
    FlankingA_strt.temp= ncread(FlankingA, 'ctdmo_seawater_temperature');
    FlankingB_strt.temp = ncread(FlankingB, 'ctdmo_seawater_temperature');
   %% 
    %salinity
    ApexS.sal= ncread(Surface,'met_salsurf');
    NearS.sal = ncread(NearSurface,'practical_salinity');
    FlankingA_strt.sal = ncread(FlankingA, 'practical_salinity');
    FlankingB_strt.sal = ncread(FlankingB, 'practical_salinity');
    %%
    %pressure 
    ApexS.depth= ncread(Surface,'ct_depth');
    NearS.p = ncread(NearSurface,'pressure');
    FlankingA_strt.p = ncread(FlankingA, 'ctdmo_seawater_pressure');
    FlankingB_strt.p = ncread(FlankingB, 'ctdmo_seawater_pressure');
    %%
    
    %Sea Surface temp
    ApexS.surfacetemp= ncread(Surface,'sea_surface_temperature');
   
    %LatAndLon
    ApexS.Lat = ncread(Surface, 'lat');
    ApexS.Lon = ncread(Surface, 'lon');
    NearS.Lat = ncread(NearSurface, 'lat');
    NearS.Lon = ncread(NearSurface, 'lon');
    FlankingA_strt.Lat = ncread(FlankingA, 'lat');
    FlankingA_strt.Lon = ncread(FlankingA, 'lon');
    FlankingB_strt.Lat = ncread(FlankingB, 'lat');
    FlankingB_strt.Lon = ncread(FlankingB, 'lon');
    
    %%
    %Expected O2
    ApexS.O2_expected = gsw_O2sol_SP_pt(ApexS.sal , ApexS.surfacetemp)
    NearS.O2_expected = gsw_O2sol_SP_pt(NearS.sal , NearS.temp)
    FlankingA_strt.O2_expected = gsw_O2sol_SP_pt(FlankingA_strt.sal , FlankingA_strt.temp)
    FlankingB_strt.O2_expected = gsw_O2sol_SP_pt(FlankingB_strt.sal , FlankingB_strt.temp)
   
    