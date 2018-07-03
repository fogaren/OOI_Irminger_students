%% Example of how to remove diurnal signal by separating out day and night indices
lat_irming = (60);
lon_irming = -39;
UTCoffset = 0;
tol = 1; %hrs before/after sunrise/sunset to count as daylight
%%
[ApexS.dayind,ApexS.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,ApexS.Time,tol);
% [Yr2_flmB.dayind,Yr2_flmB.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_flmB.time_mat,tol);
% [Yr2_sb.dayind,Yr2_sb.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_sb.time_mat,tol);
% [Yr2_rid.dayind,Yr2_rid.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_rid.time_mat,tol);

%% Example of how to plot the separated day and night data
figure(1); clf; %Surface buoy
plot(ApexS.Time(ApexS.dayind),ApexS.O2(ApexS.dayind),'b.'); 
%% hold on; %Year 2    
plot(ApexS.Time(ApezS.nightind),Yr2_sb.oxygen(Yr2_sb.nightind),'k.'); hold on; %Year 2
datetick('x'); title('Surface Buoy Oxygen (T,S,P corrected)')

%% Example of how to pull out only good data for final analysis
Yr2_sb.idgood = find(Yr2_sb.time_mat < datenum(2016,1,25)); %in this example, removed bad data based on time, but could also set other filters
Yr2_sb.O2ind = intersect(Yr2_sb.idgood,Yr2_sb.nightind);

%% Another possible step (generally done after calibration) is to smooth the data with a filter
filteringpointsO2 = 20; %remove spikes in data by taking median over every 20 points
Yr2_sb.O2cleaned = medfilt1(Yr2_sb.oxygen_calibrated,filteringpointsO2);