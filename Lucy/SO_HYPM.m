%HYPM Code for Southern Ocean for LOWER PROFILER
%Extract Year 1 data

%Wire-following profiler, Year 1, DOSTA    
filename = ['deployment0001_GS02HYPM-WFP03-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20150220T020205-20151225T000845.nc']; ncdisp(filename)
    Yr1_wfp.time_dosta = ncread(filename,'time');
    Yr1_wfp.lon_dosta = ncread(filename,'lon');
    Yr1_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr1_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_dosta = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_dosta = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr1_wfp.oxygen = ncread(filename,'dissolved_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr1_wfp.optode_temperature = ncread(filename,'optode_temperature_qc_results'); %long_name = 'Optode Temperature' units = 'deg_C'
        %Note that these points fall close to, but not exactly on, a 1:1
        %line with the CTD temperature. Points with zero value of optode
        %temperature may be indicator of bad oxygen data points.
    %Convert to matlab time
    Yr1_wfp.time_dosta_mat = convertTime(Yr1_wfp.time_dosta);

%Wire-following profiler, Year 1, Fluorometer    
filename = ['deployment0001_GS02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150220T000204-20151121T114916.nc']; ncdisp(filename)
    Yr1_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
    Yr1_wfp.lon_flord = ncread(filename,'lon');
    Yr1_wfp.lat_flord = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr1_wfp.temperature_flord = ncread(filename,'raw_internal_temp'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_flord = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_flord = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr1_wfp.backscatter = ncread(filename,'optical_backscatter'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr1_wfp.scat_total = ncread(filename,'seawater_scattering_coefficient'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr1_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
        %Note that there are other scattering based results that I don't
        %really understand - look back at these if I actually want to use
        %the backscatter or scat_total results
    %Convert to matlab time
    Yr1_wfp.time_flord_mat = convertTime(Yr1_wfp.time_flord);
        
%% Extract Year 2 data
%Wire-following profiler, Year 2, DOSTA    
filename = ['deployment0002_GS02HYPM-WFP03-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20151216T040203-20161211T105326.nc']; ncdisp(filename)
    Yr2_wfp.time_dosta = ncread(filename,'time');
    Yr2_wfp.lon_dosta = ncread(filename,'lon');
    Yr2_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr2_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr2_wfp.pracsal_dosta = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr2_wfp.pressure_dosta = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr2_wfp.oxygen = ncread(filename,'dissolved_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr2_wfp.optode_temperature = ncread(filename,'optode_temperature_qc_results'); %long_name = 'Optode Temperature' units = 'deg_C'   
    %Convert to matlab time
    Yr2_wfp.time_dosta_mat = convertTime(Yr2_wfp.time_dosta);

%Wire-following profiler, Year 2, Fluorometer    
filename = ['deployment0002_GS02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20151216T000205-20161211T064635.nc']; ncdisp(filename)
    %All of these are the same as the DOSTA points
%     Yr2_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
%     Yr2_wfp.lon_flord = ncread(filename,'lon');
%     Yr2_wfp.lat_flord = ncread(filename,'lat');
%     Yr2_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
%     Yr2_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
%     Yr2_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr2_wfp.backscatter = ncread(filename,'optical_backscatter'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr2_wfp.scat_total = ncread(filename,'seawater_scattering_coefficient'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr2_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'  
    
%% Extract Year 3 Data
%Wire-following profiler, Year 3, DOSTA    
filename = ['deployment0003_GS02HYPM-WFP03-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20161202T040203-20171130T215714.nc']; ncdisp(filename)
    Yr3_wfp.time_dosta = ncread(filename,'time');
    Yr3_wfp.lon_dosta = ncread(filename,'lon');
    Yr3_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr3_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr3_wfp.pracsal_dosta = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr3_wfp.pressure_dosta = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr3_wfp.oxygen = ncread(filename,'dissolved_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr3_wfp.optode_temperature = ncread(filename,'optode_temperature_qc_results'); %long_name = 'Optode Temperature' units = 'deg_C'
        %Note that these points fall close to, but not exactly on, a 1:1
        %line with the CTD temperature. Points with zero value of optode
        %temperature may be indicator of bad oxygen data points.
    %Convert to matlab time
    Yr3_wfp.time_dosta_mat = convertTime(Yr1_wfp.time_dosta);

%Wire-following profiler, Year 3, Fluorometer    
filename = ['deployment0003_GS02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20161202T000205-20171130T184256.nc']; ncdisp(filename)
    Yr3_wfp.time_flord = ncread(filename,'time');
    Yr3_wfp.lon_flord = ncread(filename,'lon');
    Yr3_wfp.lat_flord = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr3_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr3_wfp.pracsal_flord = ncread(filename,'practical_salinity'); %standard_name = 'sea_water_practical_salinity'
    Yr3_wfp.pressure_flord = ncread(filename,'int_ctd_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr3_wfp.backscatter = ncread(filename,'optical_backscatter'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr3_wfp.scat_total = ncread(filename,'seawater_scattering_coefficient'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr3_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
        %Note that there are other scattering based results that I don't
        %really understand - look back at these if I actually want to use
        %the backscatter or scat_total results
    %Convert to matlab time
    Yr3_wfp.time_flord_mat = convertTime(Yr1_wfp.time_flord);
    

 %% Assign profile indices prior to gridding
Yr1_wfp.depth_dosta = -gsw_z_from_p(Yr1_wfp.pressure_dosta,Yr1_wfp.lat_dosta);
    [Yr1_wfp.profile_index,Yr1_wfp.updown_index] = profileIndex(Yr1_wfp.depth_dosta);

Yr2_wfp.depth_dosta = -gsw_z_from_p(Yr2_wfp.pressure_dosta,Yr2_wfp.lat_dosta);
    [Yr2_wfp.profile_index,Yr2_wfp.updown_index] = profileIndex(Yr2_wfp.depth_dosta);
    
Yr3_wfp.depth_dosta = -gsw_z_from_p(Yr3_wfp.pressure_dosta,Yr3_wfp.lat_dosta);
    [Yr3_wfp.profile_index,Yr3_wfp.updown_index] = profileIndex(Yr3_wfp.depth_dosta);

%% Calculate density in raw profiles to enable gridding on density surfaces
Yr1_wfp.pdens = gsw_p_from_z(-Yr1_wfp.depth_dosta, Yr1_wfp.lat_dosta); %potential density function
Yr2_wfp.pdens = gsw_p_from_z(-Yr2_wfp.depth_dosta, Yr2_wfp.lat_dosta);
Yr3_wfp.pdens = gsw_p_from_z(-Yr3_wfp.depth_dosta, Yr3_wfp.lat_dosta);

%% Grid data to consistent depth intervals for each profile
depth_grid = [150:5:2600];
secinday = 60*60*24;

%All profiles for year 1
scivars = [Yr1_wfp.temperature_dosta, Yr1_wfp.pracsal_dosta, Yr1_wfp.oxygen, Yr1_wfp.optode_temperature...
        Yr1_wfp.backscatter, Yr1_wfp.scat_total, Yr1_wfp.chla];
[Yr1_wfpgrid] = glider_grid(Yr1_wfp.time_dosta,Yr1_wfp.lat_dosta,Yr1_wfp.lon_dosta,Yr1_wfp.depth_dosta,Yr1_wfp.profile_index,Yr1_wfp.updown_index',scivars,depth_grid);
    Yr1_wfpgrid.depth_grid = depth_grid;
Yr1_wfpgrid.time_start = convertTime(Yr1_wfpgrid.time_start);
Yr1_wfpgrid.duration = Yr1_wfpgrid.duration/secinday;
Yr1_wfpgrid.updown = Yr1_wfpgrid.profile_direction;

%All profiles for year 2
scivars = [Yr2_wfp.temperature_dosta, Yr2_wfp.pracsal_dosta, Yr2_wfp.oxygen, Yr2_wfp.optode_temperature...
        Yr2_wfp.backscatter, Yr2_wfp.scat_total, Yr2_wfp.chla];
[Yr2_wfpgrid] = glider_grid(Yr2_wfp.time_dosta,Yr2_wfp.lat_dosta,Yr2_wfp.lon_dosta,Yr2_wfp.depth_dosta,Yr2_wfp.profile_index,Yr2_wfp.updown_index',scivars,depth_grid);
    Yr2_wfpgrid.depth_grid = depth_grid;
Yr2_wfpgrid.time_start = convertTime(Yr2_wfpgrid.time_start);
Yr2_wfpgrid.duration = Yr2_wfpgrid.duration/secinday;
Yr2_wfpgrid.updown = Yr2_wfpgrid.profile_direction;  

%All profiles for year 3
scivars = [Yr3_wfp.temperature_dosta, Yr3_wfp.pracsal_dosta, Yr3_wfp.oxygen, Yr3_wfp.optode_temperature...
        Yr3_wfp.backscatter, Yr3_wfp.scat_total, Yr3_wfp.chla];
[Yr3_wfpgrid] = glider_grid(Yr3_wfp.time_dosta,Yr3_wfp.lat_dosta,Yr3_wfp.lon_dosta,Yr3_wfp.depth_dosta,Yr3_wfp.profile_index,Yr3_wfp.updown_index',scivars,depth_grid);
    Yr3_wfpgrid.depth_grid = depth_grid;
Yr3_wfpgrid.time_start = convertTime(Yr3_wfpgrid.time_start);
Yr3_wfpgrid.duration = Yr3_wfpgrid.duration/secinday;
Yr3_wfpgrid.updown = Yr3_wfpgrid.profile_direction; 

%% Take mean of paired up and down profiles
    tol = 1; %only combine profiles where time_start is < 1 day apart
[Yr1_wfpgrid.scivars_pair,Yr1_wfpgrid.ind_pair] = profilePairMean(Yr1_wfpgrid,tol);
[Yr2_wfpgrid.scivars_pair,Yr2_wfpgrid.ind_pair] = profilePairMean(Yr2_wfpgrid,tol);
[Yr3_wfpgrid.scivars_pair,Yr3_wfpgrid.ind_pair] = profilePairMean(Yr3_wfpgrid,tol);

%% Unpack scivars in gridded form
%When using scivars, gets all profiles (both up and down)
%When using scivars_pair, takes mean of paired up and down profiles

% Year 1
Yr1_wfpgrid.T = squeeze(Yr1_wfpgrid.scivars_pair(:,1,:));
Yr1_wfpgrid.S = squeeze(Yr1_wfpgrid.scivars_pair(:,2,:));
Yr1_wfpgrid.O2conc = squeeze(Yr1_wfpgrid.scivars_pair(:,3,:));
Yr1_wfpgrid.optode_temperature = squeeze(Yr1_wfpgrid.scivars_pair(:,4,:));
Yr1_wfpgrid.backscatter = squeeze(Yr1_wfpgrid.scivars_pair(:,5,:));
Yr1_wfpgrid.scat_total = squeeze(Yr1_wfpgrid.scivars_pair(:,6,:));
Yr1_wfpgrid.chla = squeeze(Yr1_wfpgrid.scivars_pair(:,7,:));
Yr1_wfpgrid.pdens = gsw_sigma0(Yr1_wfpgrid.S,Yr1_wfpgrid.T)+1000; 
Yr1_wfpgrid.press = gsw_p_from_z(repmat(-Yr1_wfpgrid.depth_grid,length(Yr1_wfpgrid.profile_ind),1)',...
        repmat(Yr1_wfpgrid.lat,1,length(Yr1_wfpgrid.depth_grid))');

    
% Year 2
Yr2_wfpgrid.T = squeeze(Yr2_wfpgrid.scivars_pair(:,1,:));
Yr2_wfpgrid.S = squeeze(Yr2_wfpgrid.scivars_pair(:,2,:));
Yr2_wfpgrid.O2conc = squeeze(Yr2_wfpgrid.scivars_pair(:,3,:));
Yr2_wfpgrid.optode_temperature = squeeze(Yr2_wfpgrid.scivars_pair(:,4,:));
Yr2_wfpgrid.backscatter = squeeze(Yr2_wfpgrid.scivars_pair(:,5,:));
Yr2_wfpgrid.scat_total = squeeze(Yr2_wfpgrid.scivars_pair(:,6,:));
Yr2_wfpgrid.chla = squeeze(Yr2_wfpgrid.scivars_pair(:,7,:));
Yr2_wfpgrid.pdens = gsw_sigma0(Yr2_wfpgrid.S,Yr2_wfpgrid.T)+1000; 
Yr2_wfpgrid.press = gsw_p_from_z(repmat(-Yr2_wfpgrid.depth_grid,length(Yr2_wfpgrid.profile_ind),1)',...
        repmat(Yr2_wfpgrid.lat,1,length(Yr2_wfpgrid.depth_grid))');      
    
         % Year 3
Yr3_wfpgrid.T = squeeze(Yr3_wfpgrid.scivars_pair(:,1,:));
Yr3_wfpgrid.S = squeeze(Yr3_wfpgrid.scivars_pair(:,2,:));
Yr3_wfpgrid.O2conc = squeeze(Yr3_wfpgrid.scivars_pair(:,3,:));
Yr3_wfpgrid.optode_temperature = squeeze(Yr3_wfpgrid.scivars_pair(:,4,:));
Yr3_wfpgrid.backscatter = squeeze(Yr3_wfpgrid.scivars_pair(:,5,:));
Yr3_wfpgrid.scat_total = squeeze(Yr3_wfpgrid.scivars_pair(:,6,:));
Yr3_wfpgrid.chla = squeeze(Yr3_wfpgrid.scivars_pair(:,7,:));
Yr3_wfpgrid.pdens = gsw_sigma0(Yr3_wfpgrid.S,Yr3_wfpgrid.T)+1000; 
Yr3_wfpgrid.press = gsw_p_from_z(repmat(-Yr3_wfpgrid.depth_grid,length(Yr3_wfpgrid.profile_ind),1)',...
        repmat(Yr3_wfpgrid.lat,1,length(Yr3_wfpgrid.depth_grid))'); 
   
    
%Calculate O2 saturation    
    O2equil = gsw_O2sol_SP_pt(Yr1_wfpgrid.S,Yr1_wfpgrid.T);
Yr1_wfpgrid.O2sat = (Yr1_wfpgrid.O2conc./O2equil - 1)*100;   
    O2equil = gsw_O2sol_SP_pt(Yr2_wfpgrid.S,Yr2_wfpgrid.T);
Yr2_wfpgrid.O2sat = (Yr2_wfpgrid.O2conc./O2equil - 1)*100;
    O2equil = gsw_O2sol_SP_pt(Yr3_wfpgrid.S,Yr3_wfpgrid.T);
Yr3_wfpgrid.O2sat = (Yr3_wfpgrid.O2conc./O2equil - 1)*100;

%% If using separate up and down profiles, show comparison (don't do this if using profilePairMean above)
%plotUpDownProfileComparisonWFP

%% Visualize gridded data
for i = 1:3
    if i == 1
        plotting = Yr1_wfpgrid;
    elseif i == 2
        plotting = Yr2_wfpgrid;
    elseif i == 3
        plotting = Yr3_wfpgrid;
    end


figure(i); clf
set(gcf,'color','w')
x0=1;
y0=1;
width=28;
height=20;
set(gcf,'units','centimeters','position',[x0,y0,width,height]) 
    subplot(4,2,1)
imagesc(plotting.T); colorbar; title('Temperature');
    subplot(4,2,2)
imagesc(plotting.S); colorbar; title('Salinity');
    subplot(4,2,3)
imagesc(plotting.pdens); colorbar; title('Density');
    subplot(4,2,4)
imagesc(plotting.O2conc); colorbar; caxis([240 300]); title('O_2 concentration');
    subplot(4,2,5)
imagesc(plotting.O2sat); colorbar; caxis([-25 0]); title('O_2 sat');
    subplot(4,2,6)
imagesc(plotting.backscatter); colorbar; caxis([0 0.002]); title('Backscatter');
    subplot(4,2,7)
imagesc(plotting.scat_total); colorbar; title('Scat Total');
    subplot(4,2,8)
imagesc(plotting.chla); colorbar; caxis([0 0.3]); title('Chlorophyll a');

end

%% Merge all years of data into a single dataset (no corrections)
wfpmerge.time = [Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair); Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair); Yr3_wfpgrid.time_start(Yr3_wfpgrid.ind_pair)];
wfpmerge.depth_grid = Yr1_wfpgrid.depth_grid;
wfpmerge.T = [Yr1_wfpgrid.T Yr2_wfpgrid.T Yr3_wfpgrid.T];
wfpmerge.S = [Yr1_wfpgrid.S Yr2_wfpgrid.S Yr3_wfpgrid.S];
wfpmerge.pdens = [Yr1_wfpgrid.pdens Yr2_wfpgrid.pdens Yr3_wfpgrid.pdens];
wfpmerge.O2conc = [Yr1_wfpgrid.O2conc Yr2_wfpgrid.O2conc Yr3_wfpgrid.O2conc];
wfpmerge.O2sat = [Yr1_wfpgrid.O2sat Yr2_wfpgrid.O2sat Yr3_wfpgrid.O2sat];
wfpmerge.backscatter = [Yr1_wfpgrid.backscatter Yr2_wfpgrid.backscatter Yr3_wfpgrid.backscatter];
wfpmerge.chla = [Yr1_wfpgrid.chla Yr2_wfpgrid.chla Yr3_wfpgrid.chla];

%% Plot merged data
%Adjustable parameters for plotting
    mindepth = 150; maxdepth = 2600;
    cints = 60; %number of contour intervals
    C = cmocean('Dense'); %set colormap
    C2 = cmocean('Algae'); 

%Make plotting grid
[X,Y] = meshgrid(wfpmerge.time, wfpmerge.depth_grid);

figure(i + 1); clf;
    subplot(411) %Density
cmin = 27.5; cmax = 27.8; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.pdens - 1000,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('\sigma_\theta', 'Fontsize', 12)

    subplot(412) %Temperature
cmin = 2; cmax = 6; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.T,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Temperature (deg C)', 'Fontsize', 12)

    subplot(413) %Oxygen concentration
cmin = 230; cmax = 320; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.O2conc,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration (mol/L)', 'Fontsize', 12)

    subplot(414) %Oxygen saturation
cmin = -25; cmax = 0; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.O2sat,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen saturation (%)', 'Fontsize', 12)

figure (i + 2); clf;
    subplot(211) %Backscatter
cmin = 4E-4; cmax = 1.5E-3; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.backscatter,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)'); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Backscatter')

    subplot(212) %Chlorophyll
cmin = 0; cmax = 0.3; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.chla,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)'); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Chlorophyll (�g/L)', 'Fontsize', 15)

