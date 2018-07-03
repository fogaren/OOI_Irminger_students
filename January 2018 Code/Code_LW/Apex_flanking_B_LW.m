%This line will display information about the netcdf file you have
%downloaded (change the filename it points to)
filename = ['deployment0003_GI01SUMO-RID16-06-DOSTAD000-recovered_host-dosta_abcdjm_dcl_instrument_recovered_20160710T174508.675000-20170812T204759.773000.nc']; ncdisp(filename)

%These are the names of the variables I chose to save from the file above.
%You will probably also want time, lat, and lon and (if available) the sea
%surface temperature and salinity from the METBK data, but the actual
%sensor data will be different. You will want to make sure the variable
%name in each ncread function is the correct one, based on reading through
%the information that was displayed about the file
    Yr3_rid.time = ncread(filename,'time');
    Yr3_rid.lat = ncread(filename,'lat');
    Yr3_rid.lon = ncread(filename,'lon');
    Yr3_rid.oxygen_apexns = ncread(filename,'dissolved_oxygen'); %From ncdisp informationL units = umol kg-1, long_name = 'DO - Temp Sal Corrected (METBK)'
    %CTD data - For some but not all sensors, these will be automatically taken from taken measurements by the CTD
    Yr3_rid.temperature_dosta_apexns = ncread(filename,'temp');
    Yr3_rid.pracsal_dosta_apexns = ncread(filename,'practical_salinity');
    Yr3_rid.pressure_dosta = ncread(filename,'pressure');
    %METBK data - For some sensors near the surface, these will be automatically taken from measurements by the surface meteorological package
    %Yr3_rid.SSS_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf');
    %Yr3_rid.SST_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature');
    %Convert to matlab time (make sure the function convertTime is in your
    %Matlab path)
    Yr3_rid.time_mat = convertTime(Yr3_rid.time);
    
    %This line will display information about the netcdf file you have
%downloaded (change the filename it points to)
filename = ['deployment0003_GI03FLMB-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20160713T160001-20170112T181502.nc']; ncdisp(filename)

%These are the names of the variables I chose to save from the file above.
%You will probably also want time, lat, and lon and (if available) the sea
%surface temperature and salinity from the METBK data, but the actual
%sensor data will be different. You will want to make sure the variable
%name in each ncread function is the correct one, based on reading through
%the information that was displayed about the file
    Yr3_rid.time_flankb = ncread(filename,'time');
    Yr3_rid.lat = ncread(filename,'lat');
    Yr3_rid.lon = ncread(filename,'lon');
    Yr3_rid.oxygen_flankb = ncread(filename,'dissolved_oxygen'); %From ncdisp informationL units = umol kg-1, long_name = 'DO - Temp Sal Corrected (METBK)'
    %CTD data - For some but not all sensors, these will be automatically taken from taken measurements by the CTD
    Yr3_rid.temperature_dosta_flankb = ncread(filename,'ctdmo_seawater_temperature');
    Yr3_rid.pracsal_dosta_flankb = ncread(filename,'practical_salinity');
    Yr3_rid.pressure_dosta = ncread(filename,'ctdmo_seawater_pressure');
    %METBK data - For some sensors near the surface, these will be automatically taken from measurements by the surface meteorological package
    %Yr3_rid.SSS_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf');
    %Yr3_rid.SST_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature');
    %Convert to matlab time (make sure the function convertTime is in your
    %Matlab path)
    Yr3_rid.time_mat_flankb = convertTime(Yr3_rid.time_flankb);
    
    %%
    figure(11); clf
    plot(Yr3_rid.time_mat_flankb,Yr3_rid.oxygen_flankb,'k.'); hold on;
    plot(Yr3_rid.time_mat,Yr3_rid.oxygen_apexns,'g.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    ylim([250 320]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('OOI Irminger Sea Oxygen (original factory calibration) from the Flanking B and Apex instrument frames')
    legend('Flanking B', 'Apex Near Surface')
     %%
    figure(7); clf
    plot(Yr3_rid.time_mat_flankb,Yr3_rid.pracsal_dosta_flankb,'k.'); hold on;
    plot(Yr3_rid.time_mat,Yr3_rid.pracsal_dosta_apexns,'g.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    %ylim([33 36]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Temperature (Celsius)','Fontsize',10)
    title('OOI Irminger Sea Salinity (original factory calibration) from the Flanking B and Apex instrument frames')
    legend('Flanking B Salinity', 'Apex Near Surface Salinity')
    %%
    figure(8); clf
    plot(Yr3_rid.time_mat_flankb,Yr3_rid.temperature_dosta_flankb,'k.'); hold on;
    plot(Yr3_rid.time_mat,Yr3_rid.temperature_dosta_apexns,'g.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    %ylim([33 36]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Temperature (Celsius)','Fontsize',10)
    title('OOI Irminger Temperature (original factory calibration) from the Flanking B and Apex instrument frames')
    legend('Flanking B','Apex Near Surface')
    
    %%
    figure(9); clf
    plot(Yr3_rid.time_mat_flankb,Yr3_rid.oxygen_flankb,'k.'); hold on;
    plot(Yr3_rid.time_mat_flankb,equil_flankb,'g.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    ylim([250 330]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('OOI Irminger Sea Oxygen (original factory calibration) from the Flanking B Instrument vs. Expected Oxygen')
    legend ('Flanking B', 'Expected in Equilibrium')
    
    %%
    figure(10); clf
    plot(Yr3_rid.time_mat,Yr3_rid.oxygen_apexns,'k.'); hold on;
    plot(Yr3_rid.time_mat,equil_apex,'g.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    %ylim([250 330]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('OOI Irminger Sea Oxygen (original factory calibration) from the Apex Near Surface Instrument vs. Expected Oxygen')
    legend ('Apex Near Surface', 'Expected in Equilibrium')