%% Directions to create merged file optode and suna as .xlsx
%To create merged file from the files automatically created every 30
%minutes by Roo's python underway system software:
% 1. Move all files to C:\Users\Hilary\Dropbox\Irminger5\Underway
% 2. Open a command window (type "cmd" in start menu)
% 3. Use "cd" to change directory to the correct folder (hint that you can navigate one folder a a time and can type the
% first few letters and then hit "tab" to fill in the full name)
% 4. Copy all text files in the directory into a merged single text file
% with the commands: copy suna_* irminger5_underwaySuna.txt and copy
% optode_* irminger5_underwayOptode.txt
% 5. Convert the merged .txt file to .xlsx using Excel. Open the .txt file
% selecting both tab and space delimiters, and then save as .xlsx with the
% same name.

%% Read in the merged xlsx data
%Optode headers: Timestamp	OptodeType	SerialNumber	O2 Conc (uM)	Air Sat (%)	Temp (Deg)	CalPhase (Deg)	TCPhase (Deg)	C1RPh (Deg)	C2RPh (Deg)	C1Amp (mV)	C2Amp (mV)	RawTemp (mV)
optode.merge_data = xlsread('C:\Users\emmal\Dropbox\OOI_Irminger_students\Irminger5_CruiseDocs\underway\irminger5_underwayOptode.xlsx');
optode.timein = optode.merge_data(:,1); %milliseconds since Jan 1 1970

suna.merge_data = xlsread('C:/Users/emmal/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/underway/irminger5_underwaySuna_truncated.xlsx');
suna.timein = suna.merge_data(:,1); %milliseconds since Jan 1 1970

%%% Convert the optode timestamps to datenum
    e = datenum('01-jan-1970 00:00:00');
        y_optode = datestr(e + optode.timein/1000/86400);
        y_suna = datestr(e + suna.timein/1000/86400);
    optode.time = datenum(y_optode);
    suna.time = datenum(y_suna);
    
%%% Pull out the main data from each sensor
    optode.O2raw = optode.merge_data(:,4);
    suna.NO3raw = suna.merge_data(:,5);
    
%% Filter out zeros from nitrate data
indzero = find(suna.NO3raw == 0);
suna.NO3raw(indzero) = NaN;

%% Read in ship's underway data
%%% To merge ship's GPS data: copy AR_GPS10_* irminger5_gps.csv
[gps.merge_data, gps.text_data] = xlsread('C:/Users/emmal/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/underway/irminger5_gps.csv');
gps.lat = gps.merge_data(:,1);
gps.lon = gps.merge_data(:,2);

%Pull out date and time and convert to Matlab format (for lines with real data)
    gps.date = gps.text_data(3:end,1);
    gps.timeofday = gps.text_data(3:end,2);
gps.time = NaN*ones(length(gps.merge_data),1);
indnonan = find(isnan(gps.lat) == 0);
gps.time(indnonan) = datenum(gps.date(indnonan)) + datenum(gps.timeofday(indnonan)) - datenum(2018,1,1);


%% Read in ship's underway seawater data
%%% To merge ship's underway T, S, + fluo data: copy AR_SSW* irminger5_ssw.csv
[ssw.merge_data, ssw.text_data] = xlsread('C:/Users/emmal/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/underway/irminger5_ssw.csv');
ssw.SST = ssw.merge_data(:,3); %using AML temperature for SST
ssw.SSS = ssw.merge_data(:,7); %using SBE45 TSG salinity for SSS
ssw.flr = ssw.merge_data(:,5); %fluorescence from shipboard fluorometer
ssw.flow = ssw.merge_data(:,6); %flow rate to ship's TSG and fluorometer

%Pull out date and time and convert to Matlab format (for lines with real data)
    ssw.date = ssw.text_data(3:end,1);
    ssw.timeofday = ssw.text_data(3:end,2);
ssw.time = NaN*ones(length(ssw.merge_data),1);
indnonan = find(isnan(ssw.SST) == 0);
ssw.time(indnonan) = datenum(ssw.date(indnonan)) + datenum(ssw.timeofday(indnonan)) - datenum(2018,1,1);

%% Pick time ranges for plots and analysis
mintime = datenum(2018,6,5); maxtime = datenum(2018,6,25);

%% Filter everything to an even grid to plot
    time_step = 10/24/60; %10 minutes
    begtime  = mintime; 
    endtime = maxtime;

[ssw.time_filt, ssw.data_filt] = meanTimeInterval(ssw.time, [ssw.flr ssw.SST ssw.SSS ssw.merge_data(:,3)], time_step, begtime, endtime);
    ssw.flr_filt = ssw.data_filt(:,1);
    ssw.SST_filt = ssw.data_filt(:,2);
    ssw.SSS_filt = ssw.data_filt(:,3);
    ssw.Tship_filt = ssw.data_filt(:,4);
    
[gps.time_filt, gps.data_filt] = meanTimeInterval(gps.time, [gps.lon gps.lat], time_step, begtime, endtime);
    gps.lon_filt = gps.data_filt(:,1);
    gps.lat_filt = gps.data_filt(:,2);
    
[suna.time_filt, suna.NO3_filt] = meanTimeInterval(suna.time, suna.NO3raw, time_step, begtime, endtime);
%% Isolate data from transects 1&2
maxT1=datenum(2018,6,7); minT2=datenum(2018,6,22)
    time_step = 10/24/60; %10 minutes
    begtimeT1  = mintime; 
    begtimeT2 =minT2
    endtimeT1 = maxT1;
    endtimeT2 =maxtime

    %T1
[ssw.time_filtT1, ssw.data_filtT1] = meanTimeInterval(ssw.time, [ssw.flr ssw.SST ssw.SSS ssw.merge_data(:,3)], time_step, begtimeT1, endtimeT1);
    ssw.flr_filtT1 = ssw.data_filtT1(:,1);
    ssw.SST_filtT1 = ssw.data_filtT1(:,2);
    ssw.SSS_filtT1 = ssw.data_filtT1(:,3);
    ssw.Tship_filtT1 = ssw.data_filtT1(:,4);
    
[gps.time_filtT1, gps.data_filtT1] = meanTimeInterval(gps.time, [gps.lon gps.lat], time_step, begtimeT1, endtimeT1);
    gps.lon_filtT1 = gps.data_filtT1(:,1);
    gps.lat_filtT1 = gps.data_filtT1(:,2);
    
[suna.time_filtT1, suna.NO3_filtT1] = meanTimeInterval(suna.time, suna.NO3raw, time_step, begtimeT1, endtimeT1);

    %T2
[ssw.time_filtT2, ssw.data_filtT2] = meanTimeInterval(ssw.time, [ssw.flr ssw.SST ssw.SSS ssw.merge_data(:,3)], time_step, begtimeT2, endtimeT2);
    ssw.flr_filtT2 = ssw.data_filtT2(:,1);
    ssw.SST_filtT2 = ssw.data_filtT2(:,2);
    ssw.SSS_filtT2 = ssw.data_filtT2(:,3);
    ssw.Tship_filtT2 = ssw.data_filtT2(:,4);
    
[gps.time_filtT2, gps.data_filtT2] = meanTimeInterval(gps.time, [gps.lon gps.lat], time_step, begtimeT2, endtimeT2);
    gps.lon_filtT2 = gps.data_filtT2(:,1);
    gps.lat_filtT2 = gps.data_filtT2(:,2);
    
[suna.time_filtT2, suna.NO3_filtT2] = meanTimeInterval(suna.time, suna.NO3raw, time_step, begtimeT2, endtimeT2);
%% Filter oxygen data
%Calculate derivative of O2 data
O2diff = diff(optode.O2raw);
%Remove extreme outliers
    O2diff_outlier = find(abs(O2diff) > 20);
    O2diff(O2diff_outlier) = NaN;
%Calculate spikes based on remaining outliers
O2diff_spike = find(abs(O2diff) > 2.4*nanstd(O2diff));

%Cut out a time chunk before and after each spike
cutbefore = 30; cutafter = 120; %each time point is 2 seconds
O2cut = O2diff_spike; %initialize the O2cut list with spike locations
for i = 1:length(O2diff_spike)
    cutinterval = [O2diff_spike(i) - cutbefore:O2diff_spike(i) + cutafter];
    O2cut = unique([O2cut; cutinterval']);
end
%Remove any times in O2cut after end of measurements
A = find(O2cut < length(optode.time));
O2cut = O2cut(A);

figure(1); 
    subplot(2,1,1)
plot(optode.time, optode.O2raw, 'k.'); hold on;
plot(optode.time(O2cut), optode.O2raw(O2cut), 'b.'); hold on;
plot(optode.time(O2diff_spike), optode.O2raw(O2diff_spike), 'r.'); hold on;
axis([mintime maxtime 350 500])
datetick('x', 2, 'keeplimits'); title('Oxygen concentration')

    subplot(2,1,2)
plot(optode.time(2:end), O2diff, 'k.'); hold on;
plot(optode.time(O2cut), O2diff(O2cut), 'b.'); hold on;
axis([mintime maxtime -1 3])
datetick('x', 2, 'keeplimits'); title('Oxygen diff')

%Make a new variable for O2 with spikes removed
optode.O2_nospike = optode.O2raw; optode.O2_nospike(O2cut) = NaN;

%% Apply the salinity correction
%Calculate salinity on optode times
    ind = find(isnan(ssw.time) + isnan(ssw.SSS) + isnan(ssw.SST) == 0);
optode.SSS_interp = interp1(ssw.time(ind), ssw.SSS(ind), optode.time);
optode.SST_interp = interp1(ssw.time(ind), ssw.SST(ind), optode.time);
optode.O2_nospike_salcorr = aaoptode_salpresscorr(optode.O2_nospike, optode.SST_interp, optode.SSS_interp, 0, 0);

%% Filter O2 to even grid
[optode.time_filt, optode.O2_nospike_salcorr_filt] = meanTimeInterval(optode.time, optode.O2_nospike_salcorr, time_step, begtime, endtime);
%O2 for Transects 1&2
[optode.time_filtT1, optode.O2_nospike_salcorr_filtT1] = meanTimeInterval(optode.time, optode.O2_nospike_salcorr, time_step, begtimeT1, endtimeT1);
[optode.time_filtT2, optode.O2_nospike_salcorr_filtT2] = meanTimeInterval(optode.time, optode.O2_nospike_salcorr, time_step, begtimeT2, endtimeT2);

%% Read in calibration data for O2
% [Winkler_BCP, Winkler_BCPtxt] = xlsread('C:/Users/Hilary/Dropbox/Irminger5/Irminger5_WinklerSamples.xlsx');
% Winkler_BCP = Winkler_BCP(10:end,:); Winkler_BCPtxt = Winkler_BCPtxt(10:end,:);
%     indUnderway = find(isnan(Winkler_BCP(:,7)) == 0);
%     Winkler_timeUnderway = datenum(Winkler_BCPtxt(indUnderway,6)) + Winkler_BCP(indUnderway,7);
% 
%     %To get densities for O2 data
%     time_uw = datenum(Winkler_BCPtxt(indUnderway,6)) + Winkler_BCP(indUnderway,7);
%     SST_uw = interp1(ssw.time_filt, ssw.SST_filt, time_uw);
%     SSS_uw = interp1(ssw.time_filt, ssw.SSS_filt, time_uw);
%     Tship_uw = interp1(ssw.time_filt, ssw.Tship_filt, time_uw);
%     dens_uw = sw_dens0(SSS_uw, Tship_uw);
    
%% Plot data by time
O2min = 290; O2max = 390; NO3min = -5; NO3max = 25; SSTmin = 4; SSTmax = 10; flrmin = 50; flrmax = 110;

figure(2); clf
    subplot(4,1,1)
plot(ssw.time, ssw.SST, 'k.'); hold on;
plot(ssw.time_filt, ssw.SST_filt, 'b.'); hold on;
axis([mintime maxtime SSTmin SSTmax])
datetick('x', 2, 'keeplimits'); title('SST')

    subplot(4,1,2)
plot(ssw.time, ssw.flr, 'k.'); hold on;
plot(ssw.time_filt, ssw.flr_filt, 'b.'); hold on;
axis([mintime maxtime flrmin flrmax])
datetick('x', 2, 'keeplimits'); title('fluo')

    subplot(4,1,3)
plot(optode.time, optode.O2_nospike_salcorr, 'k.'); hold on;
plot(optode.time_filt, optode.O2_nospike_salcorr_filt, 'b.'); hold on;
% plot(Winkler_timeUnderway, Winkler_BCP(indUnderway,16), 'r.','markersize',15); hold on;
axis([mintime maxtime O2min O2max])
datetick('x', 2, 'keeplimits'); title('Oxygen concentration')

    subplot(4,1,4)
plot(suna.time, suna.NO3raw, 'k.'); hold on;
plot(suna.time_filt, suna.NO3_filt, 'b.'); hold on;
axis([mintime maxtime NO3min NO3max])
datetick('x', 2, 'keeplimits'); title('Nitrate')

%% Plot ship track to date

figure(3); clf
latminplot = 59; latmaxplot = 65; lonminplot = -45; lonmaxplot = -20;
    subplot(221)
scatter(gps.lon_filt, gps.lat_filt, [], ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
axis([lonminplot lonmaxplot latminplot latmaxplot])
xlabel('Longitude'); ylabel('Latitude'); title('Time in cruise (days)')   
    subplot(222)
scatter(gps.lon_filt, gps.lat_filt, [], ssw.SST_filt, 'filled'); colorbar;
axis([lonminplot lonmaxplot latminplot latmaxplot])
xlabel('Longitude'); ylabel('Latitude'); title('Sea surface temperature (deg C)')   
    subplot(223)
scatter(gps.lon_filt, gps.lat_filt, [], ssw.flr_filt, 'filled'); colorbar;
axis([lonminplot lonmaxplot latminplot latmaxplot])
xlabel('Longitude'); ylabel('Latitude'); title('Chlorophyll fluorescence')
    subplot(224)
scatter(gps.lon_filt, gps.lat_filt, [], suna.NO3_filt); colorbar;
axis([lonminplot lonmaxplot latminplot latmaxplot])
xlabel('Longitude'); ylabel('Latitude'); title('Nitrate')
%% Emma nitrate vs. fluorometer figure TRANSECTS
figure (4); clf;
    subplot(211)
    plot(suna.NO3_filtT1, ssw.flr_filtT1, 'b.')
    xlabel('Nitrate'); ylabel('Fluo');
    axis ([-20, 20, 60, 90]); 
% Emma nitrate vs. O2 concentration figure
    subplot(212)
    plot(suna.NO3_filtT1, optode.O2_nospike_salcorr_filtT1, 'r.')
    xlabel('Nitrate'); ylabel('O2 Concentration')
    axis ([0, 20, 280, 380])
%% Emma spatial and time-based figures for fluo and nitrate LONGITUDE
figure (6); clf;
lonminplot= -43; lonmaxplot= -20; NO3minplot= -30; NO3maxplot= 20; Fminplot= 50; Fmaxplot= 110;
    subplot(211)
scatter(gps.lon_filt, suna.NO3_filt, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
axis([lonminplot lonmaxplot NO3minplot NO3maxplot])
xlabel('Longitude'); ylabel('Nitrate'); title('Nitrate by Space and Time')
    subplot(212)
scatter(gps.lon_filt, ssw.flr_filt, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
axis([lonminplot lonmaxplot Fminplot Fmaxplot])
xlabel('Longitude'); ylabel('Fluo'); title('Fluo by Space and Time')
%% Emma spatial and time-based figures for fluo and nitrate LATITUDE
figure (7); clf;
latminplot= 59; latmaxplot= 65; NO3minplot= -30; NO3maxplot= 20; Fminplot= 50; Fmaxplot= 110;
    subplot(211)
scatter(gps.lat_filt, suna.NO3_filt, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
axis([latminplot latmaxplot NO3minplot NO3maxplot])
xlabel('Latitude'); ylabel('Nitrate'); title('Nitrate by Space and Time')
    subplot(212)
scatter(gps.lat_filt, ssw.flr_filt, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
axis([latminplot latmaxplot Fminplot Fmaxplot])
xlabel('Latitude'); ylabel('Fluo'); title('Fluo by Space and Time')