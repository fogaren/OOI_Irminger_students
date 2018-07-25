%% Extract Year 1 data

%Wire-following profiler, Year 1, DOSTA    
filename = ['deployment0001_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20140912T050204-20150812T103930.nc']; ncdisp(filename)
    Yr1_wfp.time_dosta = ncread(filename,'time');
    Yr1_wfp.lon_dosta = ncread(filename,'lon');
    Yr1_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr1_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr1_wfp.oxygen = ncread(filename,'dosta_ln_wfp_abs_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr1_wfp.optode_temperature = ncread(filename,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'
        %Note that these points fall close to, but not exactly on, a 1:1
        %line with the CTD temperature. Points with zero value of optode
        %temperature may be indicator of bad oxygen data points.
    %Convert to matlab time
    Yr1_wfp.time_dosta_mat = convertTime(Yr1_wfp.time_dosta);

%Wire-following profiler, Year 1, Fluorometer    
filename = ['deployment0001_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20140912T050204-20150812T103930.nc']; ncdisp(filename)
    Yr1_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
    Yr1_wfp.lon_flord = ncread(filename,'lon');
    Yr1_wfp.lat_flord = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr1_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr1_wfp.backscatter = ncread(filename,'flort_kn_bback_total'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr1_wfp.scat_total = ncread(filename,'scat_seawater'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr1_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
        %Note that there are other scattering based results that I don't
        %really understand - look back at these if I actually want to use
        %the backscatter or scat_total results
    %Convert to matlab time
    Yr1_wfp.time_flord_mat = convertTime(Yr1_wfp.time_flord);
        
%% Extract Year 2 data
%Wire-following profiler, Year 2, DOSTA    
filename = ['deployment0002_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20150817T030206-20160628T060527.nc']; ncdisp(filename)
    Yr2_wfp.time_dosta = ncread(filename,'time');
    Yr2_wfp.lon_dosta = ncread(filename,'lon');
    Yr2_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr2_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr2_wfp.pracsal_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr2_wfp.pressure_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr2_wfp.oxygen = ncread(filename,'dosta_ln_wfp_abs_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr2_wfp.optode_temperature = ncread(filename,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'   
    %Convert to matlab time
    Yr2_wfp.time_dosta_mat = convertTime(Yr2_wfp.time_dosta);

%Wire-following profiler, Year 2, Fluorometer    
filename = ['deployment0002_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150817T030206-20160628T060527.nc']; ncdisp(filename)
    %All of these are the same as the DOSTA points
%     Yr2_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
%     Yr2_wfp.lon_flord = ncread(filename,'lon');
%     Yr2_wfp.lat_flord = ncread(filename,'lat');
%     Yr2_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
%     Yr2_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
%     Yr2_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr2_wfp.backscatter = ncread(filename,'flort_kn_bback_total'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr2_wfp.scat_total = ncread(filename,'scat_seawater'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr2_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'  

 %% Assign profile indices prior to gridding
Yr1_wfp.depth_dosta = -gsw_z_from_p(Yr1_wfp.pressure_dosta,Yr1_wfp.lat_dosta);
    [Yr1_wfp.profile_index,Yr1_wfp.updown_index] = profileIndex(Yr1_wfp.depth_dosta);

Yr2_wfp.depth_dosta = -gsw_z_from_p(Yr2_wfp.pressure_dosta,Yr2_wfp.lat_dosta);
    [Yr2_wfp.profile_index,Yr2_wfp.updown_index] = profileIndex(Yr2_wfp.depth_dosta);

%% Calculate density in raw profiles to enable gridding on density surfaces
Yr1_wfp.pdens = gsw_p_from_z(-Yr1_wfp.depth_dosta, Yr1_wfp.lat_dosta); %potential density function
Yr2_wfp.pdens = gsw_p_from_z(-Yr2_wfp.depth_dosta, Yr2_wfp.lat_dosta);

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

%% Take mean of paired up and down profiles
    tol = 1; %only combine profiles where time_start is < 1 day apart
[Yr1_wfpgrid.scivars_pair,Yr1_wfpgrid.ind_pair] = profilePairMean(Yr1_wfpgrid,tol);
[Yr2_wfpgrid.scivars_pair,Yr2_wfpgrid.ind_pair] = profilePairMean(Yr2_wfpgrid,tol);

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
Yr1_wfpgrid.pdens = sw_dens0(Yr1_wfpgrid.S,Yr1_wfpgrid.T); 
Yr1_wfpgrid.press = sw_pres(repmat(Yr1_wfpgrid.depth_grid,length(Yr1_wfpgrid.profile_ind),1)',...
        repmat(Yr1_wfpgrid.lat,1,length(Yr1_wfpgrid.depth_grid))');

    
% Year 2
Yr2_wfpgrid.T = squeeze(Yr2_wfpgrid.scivars_pair(:,1,:));
Yr2_wfpgrid.S = squeeze(Yr2_wfpgrid.scivars_pair(:,2,:));
Yr2_wfpgrid.O2conc = squeeze(Yr2_wfpgrid.scivars_pair(:,3,:));
Yr2_wfpgrid.optode_temperature = squeeze(Yr2_wfpgrid.scivars_pair(:,4,:));
Yr2_wfpgrid.backscatter = squeeze(Yr2_wfpgrid.scivars_pair(:,5,:));
Yr2_wfpgrid.scat_total = squeeze(Yr2_wfpgrid.scivars_pair(:,6,:));
Yr2_wfpgrid.chla = squeeze(Yr2_wfpgrid.scivars_pair(:,7,:));
Yr2_wfpgrid.pdens = sw_dens0(Yr2_wfpgrid.S,Yr2_wfpgrid.T); 
Yr2_wfpgrid.press = sw_pres(repmat(Yr2_wfpgrid.depth_grid,length(Yr2_wfpgrid.profile_ind),1)',...
        repmat(Yr2_wfpgrid.lat,1,length(Yr2_wfpgrid.depth_grid))');      
    
%Calculate O2 saturation    
    O2equil = O2sol(Yr1_wfpgrid.S,Yr1_wfpgrid.T);
Yr1_wfpgrid.O2sat = (Yr1_wfpgrid.O2conc./O2equil - 1)*100;   
    O2equil = O2sol(Yr2_wfpgrid.S,Yr2_wfpgrid.T);
Yr2_wfpgrid.O2sat = (Yr2_wfpgrid.O2conc./O2equil - 1)*100;

%% If using separate up and down profiles, show comparison (don't do this if using profilePairMean above)
%plotUpDownProfileComparisonWFP

%% Visualize gridded data
for i = 1:2
    if i == 1
        plotting = Yr1_wfpgrid;
    elseif i == 2
        plotting = Yr2_wfpgrid;
    end
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