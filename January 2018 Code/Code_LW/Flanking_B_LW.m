%This line will display information about the netcdf file you have
%downloaded (change the filename it points to)
path = ['C:\Users\emmal\Dropbox\OOI_Irminger_students\SST_LW_Winter_18\OOI Data_y3\DataFile\Sensor Data\'];
file = ['deployment0003_GI03FLMB-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20160713T160001-20170112T181502.nc'];
filename= [path file]
ncdisp(filename)

%These are the names of the variables I chose to save from the file above.
%You will probably also want time, lat, and lon and (if available) the sea
%surface temperature and salinity from the METBK data, but the actual
%sensor data will be different. You will want to make sure the variable
%name in each ncread function is the correct one, based on reading through
%the information that was displayed about the file
    Yr3_rid.time = ncread(filename,'time');
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
    Yr3_rid.time_mat = convertTime(Yr3_rid.time);

    %%
%This is an example of the code I used to plot the data (which you will
%need to modify to match your own variables)
figure(2); clf
    plot(Yr3_rid.time_mat,Yr3_rid.oxygen_flankb,'m*'); hold on;
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    ylim([260 300]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',15)
    title('OOI Irminger Sea Dissolved Oxygen (original factory calibration) from the Flanking B instrument frame')
    %%
    figure(5); clf
    plot(Yr3_rid.time_mat,Yr3_rid.temperature_dosta_flankb,'k.'); hold on;
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    %ylim([3 11]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Temperature (Celsius)','Fontsize',10)
    title('OOI Irminger Sea Temperature (original factory calibration) from the Flanking B instrument frame')
    
    %%
    figure(6); clf
    plot(Yr3_rid.time_mat,Yr3_rid.pracsal_dosta_flankb,'k.'); hold on;
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    ylim([33 36]) %This allows you to set limits for the min and max y axis values - I comment this out until I've looked at the data
    datetick('x',12,'keeplimits'); %This changes the x tick labels to be dates
    ylabel('Salinity','Fontsize',10)
    title('OOI Irminger Sea Salinity (original factory calibration) from the Flanking B instrument frame')
   
  