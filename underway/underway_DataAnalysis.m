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
optode.merge_data = xlsread('C:\Users\Hilary\Dropbox\Wellesley\OOI_Irminger_students\Irminger5_CruiseDocs\underway\irminger5_underwayOptode.xlsx');
optode.timein = optode.merge_data(:,1); %milliseconds since Jan 1 1970

suna.merge_data = xlsread('C:/Users/Hilary/Dropbox/Wellesley/OOI_Irminger_students/Irminger5_CruiseDocs/underway/irminger5_underwaySuna_truncated.xlsx');
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
[gps.merge_data, gps.text_data] = xlsread('C:/Users/Hilary/Dropbox/Wellesley/OOI_Irminger_students/Irminger5_CruiseDocs/underway/irminger5_gps.csv');
gps.lat = gps.merge_data(:,1);
gps.lon = gps.merge_data(:,2);

%Pull out date and time and convert to Matlab format (for lines with real data)
    gps.date = gps.text_data(3:end,1);
    gps.timeofday = gps.text_data(3:end,2);
gps.time = NaN*ones(length(gps.merge_data),1);
indnonan = find(isnan(gps.lat) == 0);
gps.time(indnonan) = datenum(gps.date(indnonan)) + datenum(gps.timeofday(indnonan)) - datenum(2019,1,1);

%% Read in ship's underway seawater data
%%% To merge ship's underway T, S, + fluo data: copy AR_SSW* irminger5_ssw.csv
[ssw.merge_data, ssw.text_data] = xlsread('C:/Users/Hilary/Dropbox/Wellesley/OOI_Irminger_students/Irminger5_CruiseDocs/underway/irminger5_ssw.csv');
ssw.SST = ssw.merge_data(:,3); %using AML temperature for SST
ssw.SSS = ssw.merge_data(:,7); %using SBE45 TSG salinity for SSS
ssw.flr = ssw.merge_data(:,5); %fluorescence from shipboard fluorometer
ssw.flow = ssw.merge_data(:,6); %flow rate to ship's TSG and fluorometer

%Pull out date and time and convert to Matlab format (for lines with real data)
    ssw.date = ssw.text_data(3:end,1);
    ssw.timeofday = ssw.text_data(3:end,2);
ssw.time = NaN*ones(length(ssw.merge_data),1);
indnonan = find(isnan(ssw.SST) == 0);
ssw.time(indnonan) = datenum(ssw.date(indnonan)) + datenum(ssw.timeofday(indnonan)) - datenum(2019,1,1);

%% Pick time ranges for plots and analysis
mintime = datenum(2018,6,5); maxtime = datenum(2018,6,25);

%% Filter everything to an even grid to plot
    time_step = 10/24/60; %10 minutes
    begtime = mintime; 
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

%% Remove cleaning times from nitrate 
clean1t = [datenum(2018,6,15,16,15,00) datenum(2018,6,15,17,45,00)]; %Cleaned SUNA-only 
clean2t = [datenum(2018,6,19,12,00,00) datenum(2018,6,19,13,30,00)]; %tube cleaning

clean1 = find(suna.time_filt>clean1t(1) & suna.time_filt<clean1t(2));
clean2 = find(suna.time_filt>clean2t(1) & suna.time_filt<clean2t(2));
clean = [clean1, clean2]; %all NO3 samples taken during cleaning times

    suna.NO3_beforeclean = suna.NO3_filt; %to use samples INCLUDING  cleaning times, use .NO3beforeclean
    suna.NO3_filt(clean) = NaN; %now suna.NO3_filt omits cleaning times
    
%% Suggested discrete samples to analyze and SUNA cleaning times
discrete_suna = [datenum(2018,6,5,19,27,00) datenum(2018,6,6,6,42,00) datenum(2018,6,8,6,26,00) datenum(2018,6,9,8,25,00) ...
                datenum(2018,6,10,6,31,00) datenum(2018,6,12,6,31,00) datenum(2018,6,12,19,58,00) datenum(2018,6,13,18,27,00) ...
                datenum(2018,6,14,6,39,00) datenum(2018,6,15,17,40,00) datenum(2018,6,16,12,38,00) datenum(2018,6,17,18,52,00) ...
                datenum(2018,6,18,20,16,00) datenum(2018,6,19,12,42,00) datenum(2018,6,22,17,54,00) datenum(2018,6,23,6,59,00)];
discrete_suna_y = [16*ones(length(discrete_suna))];

figure(1); clf; %Suggested nutrient samples to analyze
plot(suna.time_filt, suna.NO3_filt, 'm.'); hold on; %NO3 from underway 
plot(discrete_suna, discrete_suna_y, 'rx'); %discrete samples
axis([mintime maxtime -15 20]); datetick('x',2, 'keeplimits');  
ylabel('NO_3', 'Fontsize', 15); title('NO_3/nutrients: suggested discrete samples');
%% Define Sections for Time & Spatial Analysis   
    %times = [begtime arriveOOI CTD010 CTD012 startWest CTD020 leaveOOI underwayoff]
times = [begtime datenum(2018,6,8,5,00,00) datenum(2018,6,12,19,51,00) datenum(2018,6,13,20,30,00) datenum(2018,6,16,15,00,00) datenum(2018,6,19,13,49,00) datenum(2018,6,21,15,00,00) datenum(2018,6,23,17,40,00)];
    T1 = find(ssw.time_filt>times(1) & ssw.time_filt<times(2)); %Transect out to array
    array1 = [find(ssw.time_filt>times(2) & ssw.time_filt<times(3)) , find(ssw.time_filt>times(4) & ssw.time_filt<times(5))];  
    OSNAPeast = find(ssw.time_filt>times(3) & ssw.time_filt<times(4));
    OSNAPwest = find(ssw.time_filt>times(5) & ssw.time_filt<times(6));
    array2 = find(ssw.time_filt>times(6) & ssw.time_filt<times(7));    
    T2 = find(ssw.time_filt>times(7) & ssw.time_filt<times(8));
transect = [T1 T2]; %all indices from transects
array = [array1 array2]; %all indices from array

figure(11); clf; %Shows gps track with sections color coded
    for i = 1:6
        sections = {T1 T2 OSNAPeast OSNAPwest array1 array2}; %Sections 
        sectionC = {nicecolor('Rw') nicecolor('rk') nicecolor('Bbw') nicecolor('GY') nicecolor('BRw') nicecolor('rBk')}; %Section colors
        scatter(gps.lon_filt(sections{i}), gps.lat_filt(sections{i}), [], sectionC{i}, 'filled', 'LineWidth', .75); hold on;
    end 
    axis([-42 -30 59.5 63]); legend({'Transect 1', 'Transect 2', 'OSNAP East', 'OSNAP West', 'OOI Array 1', 'OOI Array 2'}, 'Fontsize', 16);
    xlabel('Longitude', 'Fontsize', 15); ylabel('Latitude', 'Fontsize', 15); title('Sections of Cruise Track', 'Fontsize', 15);

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

figure(2); clf; %Plot optode data showing spikes to remove
    subplot(2,1,1)
plot(optode.time, optode.O2raw, 'k.'); hold on;
plot(optode.time(O2cut), optode.O2raw(O2cut), 'b.'); hold on;
plot(optode.time(O2diff_spike), optode.O2raw(O2diff_spike), 'r.'); hold on;
axis([mintime maxtime 320 550])
datetick('x', 2, 'keeplimits'); title('Oxygen concentration')
    subplot(2,1,2)
plot(optode.time(2:end), O2diff, 'k.'); hold on;
plot(optode.time(O2cut), O2diff(O2cut), 'b.'); hold on;
axis([mintime maxtime -2 4])
datetick('x', 2, 'keeplimits'); title('Oxygen diff')

%Make a new variable for O2 with spikes removed
optode.O2_nospike = optode.O2raw; optode.O2_nospike(O2cut) = NaN;

%% Apply the salinity correction
%Calculate salinity on optode times
    indT1 = find(isnan(ssw.time) + isnan(ssw.SSS) + isnan(ssw.SST) == 0);
optode.SSS_interp = interp1(ssw.time(indT1), ssw.SSS(indT1), optode.time);
optode.SST_interp = interp1(ssw.time(indT1), ssw.SST(indT1), optode.time);
optode.O2_nospike_salcorr = aaoptode_salpresscorr(optode.O2_nospike, optode.SST_interp, optode.SSS_interp, 0, 0);

%% Filter O2 to even grid
[optode.time_filt, optode.O2_nospike_salcorr_filt] = meanTimeInterval(optode.time, optode.O2_nospike_salcorr, time_step, begtime, endtime);

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
    
%% Calculate O2sol and AOU
optode.O2sol_filt = gsw_O2sol_SP_pt(ssw.SSS_filt, ssw.SST_filt); %calculates O2sol using PSU and SST
optode.aou_filt = (optode.O2sol_filt)-(optode.O2_nospike_salcorr_filt); %calculates AOU
    AOU = optode.aou_filt; %AOU including outliers 
    outAOU = find(optode.aou_filt>50);
    optode.aou_filt(outAOU) = NaN; %optode.aou_filt now OMITS outliers over 50
optode.O2Sat = (((optode.O2_nospike_salcorr_filt)-(optode.O2sol_filt))./(optode.O2sol_filt))*100; %Percent saturation
AOUmin= -100; AOUmax= 10; 

figure(3);clf; %Plot O2_sol and percent saturation
    yyaxis left 
        plot(optode.time_filt, optode.O2sol_filt); hold on;
        ylabel ('O_2');
    yyaxis right
        plot(optode.time_filt, optode.O2Sat); hold on;
        axis([mintime maxtime 0 30]);
        ylabel ('Percent Saturation')
    legend('O_2 solubility', 'Percent Saturation');
    datetick('x', 2, 'keeplimits'); title('Oxygen');
%% Correlations between NO3 and AOU
for i = 1:6
    x = find(isfinite(optode.aou_filt(sections{i}))==1); %finds only real numbers of sections
        if i == 2
            x = find((isfinite(optode.aou_filt(T2))==1) & suna.NO3_filt(T2)>-5); %only real numbers, omits outliers
        end
            % ^this takes out outliers>(-5)NO3
        if i == 5,6;
            x = find((isfinite(suna.NO3_filt(sections{i}))==1) & (optode.aou_filt(sections{i})) > ((i-1).*(-11)));
        end
            % ^this takes out outliers>(-44 if array1) or (-33 if array 2)
    P(i,:) = polyfit(suna.NO3_filt(sections{i}(x)), optode.aou_filt(sections{i}(x)), 1); %variable P created to polyfit NO3 by AOU  
end

sectionL = {'Transect 1', 'Transect 2', 'OSNAP East', 'OSNAP West', 'OOI Array 1', 'OOI Array 2'};
XtoPlot = [-5:25];
figure(10); clf; % Figure showing sections with NO3 vs AOU correlations as lines
for i = 1:6
    subplot(2,3,i)
    scatter(suna.NO3_filt(sections{i}), optode.aou_filt(sections{i}), 25, sectionC{i}, 'filled'); hold on;
    plot(XtoPlot, (P(i,1))*XtoPlot + (P(i,2)), 'r--', 'LineWidth', 1); 
    ylabel(sectionL{i}, 'FontSize', 16);
    title(sectionL{i}, 'FontSize', 16);
end

%% Plot data by time
O2min = 290; O2max = 390; NO3min = -10; NO3max = 25; SSTmin = 4; SSTmax = 10; flrmin = 50; flrmax = 110;

figure(4); clf
    subplot(4,1,1) %SST
        plot(ssw.time, ssw.SST, 'k.'); hold on;
        plot(ssw.time_filt, ssw.SST_filt, 'b.'); hold on;
        axis([mintime maxtime SSTmin SSTmax])
        datetick('x', 2, 'keeplimits'); title('SST')
    subplot(4,1,2) %Fluo
        plot(ssw.time, ssw.flr, 'k.'); hold on;
        plot(ssw.time_filt, ssw.flr_filt, 'b.'); hold on;
        axis([mintime maxtime flrmin flrmax])
        datetick('x', 2, 'keeplimits'); title('fluo')
    subplot(4,1,3) %O2
        plot(optode.time, optode.O2_nospike_salcorr, 'k.'); hold on;
        plot(optode.time_filt, optode.O2_nospike_salcorr_filt, 'b.'); hold on;
        % plot(Winkler_timeUnderway, Winkler_BCP(indUnderway,16), 'r.','markersize',15); hold on;
        axis([mintime maxtime O2min O2max])
        datetick('x', 2, 'keeplimits'); title('Oxygen concentration')
    subplot(4,1,4) %NO3
        plot(suna.time, suna.NO3raw, 'k.'); hold on;
        plot(suna.time_filt, suna.NO3_filt, 'b.'); hold on;
        axis([mintime maxtime NO3min NO3max])
        datetick('x', 2, 'keeplimits'); title('Nitrate')

%% Plot ship track to date
figure(5); clf
latminplot = 59; latmaxplot = 65; lonminplot = -45; lonmaxplot = -20;
    subplot(221) %days
        scatter(gps.lon_filt, gps.lat_filt, [], ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
        axis([lonminplot lonmaxplot latminplot latmaxplot])
        xlabel('Longitude'); ylabel('Latitude'); title('Time in cruise (days)')   
    subplot(222) %SST
        scatter(gps.lon_filt, gps.lat_filt, [], ssw.SST_filt, 'filled'); colorbar;
        axis([lonminplot lonmaxplot latminplot latmaxplot])
        xlabel('Longitude'); ylabel('Latitude'); title('Sea surface temperature (deg C)')   
    subplot(223) %Fluo
        scatter(gps.lon_filt, gps.lat_filt, [], ssw.flr_filt, 'filled'); colorbar;
        axis([lonminplot lonmaxplot latminplot latmaxplot])
        xlabel('Longitude'); ylabel('Latitude'); title('Chlorophyll fluorescence')
    subplot(224) %NO3
        scatter(gps.lon_filt, gps.lat_filt, [], suna.NO3_filt, 'filled'); colorbar;
        axis([lonminplot lonmaxplot latminplot latmaxplot])
        xlabel('Longitude'); ylabel('Latitude'); title('Nitrate')
%% Spatial and temporal patterns for fluo, nitrate, O2, LONGITUDE
figure (7); clf;
NO3minplot= -30; NO3maxplot= 20;
    subplot(311) %NO3
        scatter(gps.lon_filt, suna.NO3_filt, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
        axis([lonminplot lonmaxplot NO3minplot NO3maxplot])
        xlabel('Longitude'); ylabel('Nitrate'); title('Nitrate by Space and Time')
    subplot(312) %Fluo
        scatter(gps.lon_filt, ssw.flr_filt, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
        axis([lonminplot lonmaxplot flrmin flrmax])
        xlabel('Longitude'); ylabel('Fluo'); title('Fluo by Space and Time')
    subplot(313) %O2
        scatter(gps.lon_filt, optode.O2Sat, 10, ssw.time_filt - min(ssw.time_filt), 'filled'); colorbar;
        axis([lonminplot lonmaxplot 0 20])
        xlabel('Longitude'); ylabel('O2 Sat'); title('O2 Saturation by Space and Time')
%% Emma Spatial Ratio Comparisons
figure (9); clf; 
subplot(221)
    for i = [5,1,2,3,4,6] %better to plot out of order so that array1 doesn't cover the others
        scatter(suna.NO3_filt(sections{i}), ssw.flr_filt(sections{i}), 15, sectionC{i}, 'filled'); hold on;
    end
    xlabel('Nitrate'); ylabel('Fluo'); axis ([NO3min NO3max flrmin flrmax]); title('NO_3 vs. Fluo');
subplot(222)
    for i = [5,1,2,3,4,6]
        scatter(suna.NO3_filt(sections{i}), optode.aou_filt(sections{i}), 15, sectionC{i}, 'filled'); hold on;
    end
    xlabel('Nitrate'); ylabel('AOU'); axis ([NO3min NO3max AOUmin AOUmax]); title('NO_3 vs. AOU')
subplot(223)
    for i = [5,1,2,3,4,6]
        scatter(ssw.flr_filt(sections{i}), optode.aou_filt(sections{i}), 15, sectionC{i}, 'filled'); hold on;
    end
    xlabel('Fluo'); ylabel('AOU'); axis ([flrmin flrmax AOUmin AOUmax]); title('Fluo vs. AOU')
%% Average O2sats boxplot work
x = [optode.O2Sat(transect); optode.O2Sat(array); optode.O2Sat(array2); optode.O2Sat(array1); optode.O2Sat(OSNAPwest); ...
         optode.O2Sat(OSNAPeast); optode.O2Sat(T2); optode.O2Sat(T1)];
g = [zeros(length(transect),1); ones(length(array),1); 2*ones(length(array2),1); 3*ones(length(array1),1); 4*ones(length(OSNAPwest),1); 5*ones(length(OSNAPeast),1); 6*ones(length(T2),1); 7*ones(length(T1),1)];

figure (12); clf;
    boxplot(x, g,'symbol', 'y.', 'Labels', {'Transects', 'OOI Array','Array2' 'Array1', 'OSNAP West', 'OSNAP East','Transect2', 'Transect1' ,}); 
%Commented out for now; have not found a way to fill with color AND have median lines
%     h = findobj(gca,'Tag','Box'); 
%     for j=1:length(h)
%        patch(get(h(j),'XData'),get(h(j),'YData'), sectionC{j}); % <this is where you add the color at the end
%     end   
    ylabel('Oxygen Saturation', 'FontSize', 20);
%% AOU and Nitrate over time
figure(13); clf;
    fig = figure(13);
    left_color = nicecolor('RRk'); right_color = nicecolor('RBk');
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    yyaxis left %AOU
        scatter(optode.time_filt, optode.aou_filt, '.'); hold on;
        axis([mintime maxtime AOUmin AOUmax]);
        ylabel ('AOU');
        datetick('x', 2, 'keeplimits'); title('Oxygen');
    yyaxis right %NO3
        scatter(suna.time_filt, suna.NO3_filt, '.');
        axis([mintime maxtime NO3min NO3max]);
        ylabel ('NO_3')
%% NO3 vs AOU Shelf investigation
figure(14); clf;
subplot(221)
    scatter(ssw.SST_filt(OSNAPwest), ssw.SSS_filt(OSNAPwest), [], ssw.time_filt(OSNAPwest) - min(ssw.time_filt(OSNAPwest)), 'filled'); colorbar;    
    xlabel('SST'); ylabel('Salinity')
    %axis ([NO3min 10 -70 AOUmax]);
    title('SST vs. Salinity vs. Days');
subplot(222)
    scatter(ssw.SST_filt(OSNAPwest), ssw.SSS_filt(OSNAPwest), [], ssw.flr_filt(OSNAPwest), 'filled'); colorbar;    
    xlabel('SST'); ylabel('Salinity')
    %axis ([NO3min 10 -70 AOUmax]);
    title('SST vs. Salinity vs. Fluo');
subplot(223)
    scatter(ssw.SST_filt(OSNAPwest), ssw.SSS_filt(OSNAPwest), [], suna.NO3_filt(OSNAPwest), 'filled'); colorbar;    
     xlabel('SST'); ylabel('Salinity')
    %axis ([-70 AOUmax flrmin flrmax]);
    title('SST vs. Salinity vs. Nitrate');
subplot(224)
    scatter(ssw.SST_filt(OSNAPwest), ssw.SSS_filt(OSNAPwest), [], optode.aou_filt(OSNAPwest), 'filled'); colorbar;    
    xlabel('SST'); ylabel('Salinity')
    %axis ([-70 AOUmax flrmin flrmax]);
    title('SST vs. Salinity vs. AOU');
    
figure(16); clf;
subplot(221)
    scatter(suna.NO3_filt(OSNAPwest), optode.aou_filt(OSNAPwest), [], ssw.flr_filt(OSNAPwest),'filled'); colorbar;    
    xlabel('Nitrate'); ylabel('AOU')
    %axis ([NO3min 10 -70 AOUmax]);
    title('Nitrate vs. AOU vs. Fluo');
subplot(222)
    scatter(suna.NO3_filt(OSNAPwest), optode.aou_filt(OSNAPwest), [], ssw.SST_filt(OSNAPwest), 'filled'); colorbar;    
    xlabel('Nitrate'); ylabel('AOU')
    %axis ([NO3min 10 -70 AOUmax]);
    title('Nitrate vs. AOU vs. SST');
subplot(223)
    scatter(suna.NO3_filt(OSNAPwest), optode.aou_filt(OSNAPwest), [],ssw.SSS_filt(OSNAPwest), 'filled'); colorbar;    
    xlabel('Nitrate'); ylabel('AOU')
    %axis ([-70 AOUmax flrmin flrmax]);
    title('Nitrate vs. AOU vs. Salinity');
subplot(224)
    scatter(suna.NO3_filt(OSNAPwest), optode.aou_filt(OSNAPwest), [], ssw.time_filt(OSNAPwest) - min(ssw.time_filt(OSNAPwest)), 'filled'); colorbar;    
    xlabel('Nitrate'); ylabel('AOU')
    %axis ([-70 AOUmax flrmin flrmax]);
    title('Nitrate vs. AOU vs. Days');
    
%% Sensors over time COLORED BY SECTION
figure(15); clf
    subplot(4,1,1)
scatter(ssw.time_filt(T1), ssw.SST_filt(T1), 15, nicecolor('Rw'), 'filled'); hold on;
scatter(ssw.time_filt(T2), ssw.SST_filt(T2), 15,nicecolor('rk'), 'filled'); hold on;
scatter(ssw.time_filt(array1), ssw.SST_filt(array1), 15,nicecolor('BRw'), 'filled'); hold on;
scatter(ssw.time_filt(array2), ssw.SST_filt(array2), 15,nicecolor('rBk'), 'filled'); hold on;
scatter(ssw.time_filt(OSNAPwest), ssw.SST_filt(OSNAPwest), 15,nicecolor('GY'), 'filled'); hold on;
scatter(ssw.time_filt(OSNAPeast), ssw.SST_filt(OSNAPeast), 15 ,nicecolor('Bbw'), 'filled'); hold on;
axis([mintime maxtime SSTmin SSTmax])
datetick('x', 2, 'keeplimits'); title('Sea Surface Temperature', 'FontSize', 20)

    subplot(4,1,2)
scatter(ssw.time_filt(T1), ssw.flr_filt(T1), 15, nicecolor('Rw'), 'filled'); hold on;
scatter(ssw.time_filt(T2), ssw.flr_filt(T2), 15, nicecolor('rk'), 'filled'); hold on;
scatter(ssw.time_filt(array1), ssw.flr_filt(array1), 15, nicecolor('BRw'), 'filled'); hold on;
scatter(ssw.time_filt(array2), ssw.flr_filt(array2),15, nicecolor('rBk'), 'filled'); hold on;
scatter(ssw.time_filt(OSNAPwest), ssw.flr_filt(OSNAPwest), 15,nicecolor('GY'), 'filled'); hold on;
scatter(ssw.time_filt(OSNAPeast), ssw.flr_filt(OSNAPeast), 15,nicecolor('Bbw'), 'filled'); hold on;
axis([mintime maxtime flrmin flrmax])
datetick('x', 2, 'keeplimits'); title('Chlorophyll Fluorescence', 'FontSize', 20)

    subplot(4,1,3)
scatter(optode.time_filt(T1), optode.O2_nospike_salcorr_filt(T1), 15,nicecolor('Rw'), 'filled'); hold on;
scatter(optode.time_filt(T2), optode.O2_nospike_salcorr_filt(T2), 15, nicecolor('rk'), 'filled'); hold on;
scatter(optode.time_filt(array1), optode.O2_nospike_salcorr_filt(array1), 15, nicecolor('BRw'), 'filled'); hold on;
scatter(optode.time_filt(array2),optode.O2_nospike_salcorr_filt(array2), 15, nicecolor('rBk'), 'filled'); hold on;
scatter(optode.time_filt(OSNAPwest), optode.O2_nospike_salcorr_filt(OSNAPwest), 15, nicecolor('GY'), 'filled'); hold on;
scatter(optode.time_filt(OSNAPeast), optode.O2_nospike_salcorr_filt(OSNAPeast), 15, nicecolor('Bbw'), 'filled'); hold on;
axis([mintime maxtime O2min O2max])
datetick('x', 2, 'keeplimits'); title('Oxygen Concentration', 'FontSize', 20)

    subplot(4,1,4) %INCLUDES CLEANING TIMES STILL
scatter(suna.time_filt(T1), suna.NO3_filt(T1), 15, nicecolor('Rw'), 'filled'); hold on;
scatter(suna.time_filt(T2), suna.NO3_filt(T2), 15, nicecolor('rk'), 'filled'); hold on;
scatter(suna.time_filt(array1), suna.NO3_filt(array1), 15, nicecolor('BRw'), 'filled'); hold on;
scatter(suna.time_filt(array2),suna.NO3_filt(array2), 15, nicecolor('rBk'), 'filled'); hold on;
scatter(suna.time_filt(OSNAPwest), suna.NO3_filt(OSNAPwest), 15, nicecolor('GY'), 'filled'); hold on;
scatter(suna.time_filt(OSNAPeast), suna.NO3_filt(OSNAPeast), 15, nicecolor('Bbw'), 'filled'); hold on;
axis([mintime maxtime NO3min NO3max])
datetick('x', 2, 'keeplimits'); title('Nitrate', 'FontSize', 20)
    