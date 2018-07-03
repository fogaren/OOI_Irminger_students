%% Removing diurnal signal by separating out day and night indices
lat_irming = (60);
lon_irming = -39;
UTCoffset = 0;
tol = 1; %hrs before/after sunrise/sunset to count as daylight

%% Adding day and night time indices to data structures 
[ApexS.dayind,ApexS.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,ApexS.Time,tol);
[NearS.dayind,NearS.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,NearS.Time,tol);
[FlankingA_strt.dayind,FlankingA_strt.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,FlankingA_strt.Time,tol);
[FlankingB_strt.dayind,FlankingB_strt.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,FlankingB_strt.time,tol);

%% Plotting Filtered Data
figure(1); clf; %Surface buoy
plot(ApexS.Time(ApexS.dayind),ApexS.O2(ApexS.dayind),'b.');hold on; %Year 2    
plot(ApexS.Time(ApexS.nightind),ApexS.O2(ApexS.nightind),'k.'); hold on; %Year 2
datetick('x'); title('Surface Buoy Oxygen')
legend('Daytime','Nighttime')
title('Daytime and Nighttime Data, Apex mooring') 



%%
figure(2); clf; %Surface buoy, near surface
plot(NearS.Time(NearS.dayind),NearS.O2(NearS.dayind),'b.');hold on; %Year 2    
plot(NearS.Time(NearS.nightind),NearS.O2(NearS.nightind),'k.'); hold on; %Year 2
datetick('x'); title('Near Surface Buoy Oxygen ')
legend('Daytime','Nighttime')
title('Daytime and Nighttime Data, Near Surface') 

%% 
figure(3); clf; %Flanking A 
plot(FlankingA_strt.Time(FlankingA_strt.dayind),FlankingA_strt.O2(FlankingA_strt.dayind),'b.');hold on; %Year 2    
plot(FlankingA_strt.Time(FlankingA_strt.nightind),FlankingA_strt.O2(FlankingA_strt.nightind),'k.'); hold on; %Year 2
datetick('x'); title('Flanking Mooring A Oxygen')
legend('Daytime','Nighttime')
title('Daytime and Nighttime Data, flanking mooring A') 

%%
figure(4); clf; %Flanking B
plot(FlankingB_strt.time(FlankingB_strt.dayind),FlankingB_strt.O2(FlankingB_strt.dayind),'b.');hold on; %Year 2    
plot(FlankingB_strt.time(FlankingB_strt.nightind),FlankingB_strt.O2(FlankingB_strt.nightind),'k.'); hold on; %Year 2
datetick('x'); title('Flanking Mooring B Oxygen')
legend('Daytime','Nighttime')
title('Daytime and Nighttime Data, flanking mooring B') 
  

%% Disssolved Oxygen, All data and just nighttime data
figure(5); clf; % Overlay nighttime data 
%all of data first 
ax1 = subplot(2,2,1) %This is now just the Apex Mooring 
    plot(ApexS.Time,ApexS.O2,'k.'); hold on;
    plot(NearS.Time,NearS.O2,'g.'); hold on;
    %plot(FlankingA.Time, FlankingA.O2, 'm.'); hold on;
    %plot(FlankingB.time, FlankingB.02, 'b.');
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) 
    %ylim([200 400]) 
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('All data for Apex Mooring')
    legend('Apex Surface', 'Apex Near Surface')%, 'Flanking A', 'Flanking B')
    
 ax2 = subplot(2,2,3) %This is the data from the flanking mooringss
    plot(FlankingA_strt.Time, FlankingA_strt.O2, 'm.'); hold on;
    plot(FlankingB_strt.time, FlankingB_strt.O2, 'b.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) 
    ylim([200 400]) 
    datetick('x',12, 'keeplimits');
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('All data for Flanking A & B')
    legend('Flanking A', 'Flanking B')
    
 ax3 =  subplot(2,2,2)
      plot(ApexS.Time(ApexS.nightind),ApexS.O2(ApexS.nightind),'k.'); hold on;
      plot(NearS.Time(NearS.nightind),NearS.O2(NearS.nightind),'g.'); %hold on;
      %ylim([200 400]) 
      datetick('x',12);
      ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
      title('Nighttime data for Apex Mooring')
      %plot(FlankingA_strt.Time(FlankingA_strt.nightind),FlankingA_strt.O2(FlankingA_strt.nightind),'k.'); hold on;
      %plot(FlankingB_strt.time(FlankingB_strt.nightind),FlankingB_strt.O2(ApexS.nightind),'k.'); hold on;

ax4 = subplot(2,2,4)
    plot(ApexS.Time(ApexS.nightind), ApexS.O2(ApexS.nightind),'m.'); hold on;
    plot(FlankingB_strt.time(FlankingB_strt.nightind),FlankingB_strt.O2(FlankingB_strt.nightind),'b.')
    ylim([200 400]) 
    datetick('x',12, 'keeplimits');
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Nightime data for Flanking A & B')
    legend('Flanking A', 'Flanking B')
    
    linkaxes([ax1,ax2,ax3, ax4],'xy')
    
    
%% Beginning of Cleaning the bad values 
 
  %% Removing bad data from O2

FlankingA_strt.valcleaned = BadData(FlankingA_strt.O2,266,400)  
FlankingB_strt.valclearned = BadData(FlankingB_strt.O2,266,400)


%% Removing bad data from CTD
    %When we looked the figures from the CTD data for the Y3, it showed
    %that there were certain points where there was sensor error. Future
    %work should remove these errors. 
  
      
    


