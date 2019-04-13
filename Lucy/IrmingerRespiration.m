%% Script runs in data analysis pipeline after master_oxygencalib

%% Filtering (smoothing) the data at each depth
NumProfilesToSmooth = 15; %set this as a variable so that you can adjust this value

%Updated this to smooth over time, which required adding the dimension "2"
%since you want to smooth over the 2nd dimension (time) rather than the 1st
%dimension (depth)
oxygen_driftcorr_smoothed = movmean(wfpmerge.oxygen_driftcorr, NumProfilesToSmooth, 2); %unsure at what interval I want to smooth oxygen data

%Example plot
figure(1); clf
    depth_id = 71; %Example plot at the 11th depth in depth_grid (which is 400 m), in depth_grid, 400m is 51st column
plot(wfpmerge.time, wfpmerge.oxygen_driftcorr(depth_id,:),'k.'); hold on;
plot(wfpmerge.time, oxygen_driftcorr_smoothed(depth_id,:),'m-'); hold on;
%plot(wfpmerge.time, max_O2_season2(1,71), 'b.'); hold on;
datetick('x',2)

%% Second filtering pass
%returns true for all elements more than three standard deviations from the
%mean. Along the second dimension which should calculate it along rows and
%not columns

    %TF=isoutlier(oxygen_driftcorr_smoothed,'movmedian', 5, 2);
    B=filloutliers(oxygen_driftcorr_smoothed,'previous', 2);

    %Added another smoothing pass after having removed the outliers
    oxygen_driftcorr_nooutliers_smoothed = movmean(B, NumProfilesToSmooth, 10);
    

%% Plotting at range of depths through water column
for i=10:20:200
    figure(i); 
    depth_id=i;
    plot(wfpmerge.time, wfpmerge.oxygen_driftcorr(depth_id,:),'k.'); hold on; %modified this to show original data
    plot(wfpmerge.time, oxygen_driftcorr_nooutliers_smoothed(depth_id,:),'m-'); hold on;
    datetick('x',2)
    title(['Oxygen time series at ' num2str(wfpmerge.depth_grid(i)) ' meters'])
end

%% Identifying 3 strat seasons and pulling out dates for each season
    %Time constraints to separately specify range of dates for O2 maximum
    %during winter ventilation and O2 minimum at end of stratified season
    %Note that you will want to fine-tune these date choices based on
    %plotting and looking at the results of these calculations
strat_beg_1_id = find(wfpmerge.time <= datenum(datetime(2015,6,1)) & wfpmerge.time >= datenum(datetime(2015,2,1)));
strat_end_1_id = find(wfpmerge.time <= datenum(datetime(2016,3,15)) & wfpmerge.time >= datenum(datetime(2015,11,1)));
oxygen_strat_beg_1 = oxygen_driftcorr_nooutliers_smoothed(:,strat_beg_1_id);
oxygen_strat_end_1 = oxygen_driftcorr_nooutliers_smoothed(:,strat_end_1_id);
    strat_beg_1_time = wfpmerge.time(strat_beg_1_id);
    strat_end_1_time = wfpmerge.time(strat_end_1_id);
    
strat_beg_2_id = find(wfpmerge.time <= datenum(datetime(2016,6,1)) & wfpmerge.time >= datenum(datetime(2016,2,1)));
strat_end_2_id = find(wfpmerge.time <= datenum(datetime(2017,3,10)) & wfpmerge.time >= datenum(datetime(2016,12,20)));
oxygen_strat_beg_2 = oxygen_driftcorr_nooutliers_smoothed(:,strat_beg_2_id);
oxygen_strat_end_2 = oxygen_driftcorr_nooutliers_smoothed(:,strat_end_2_id);
    strat_beg_2_time = wfpmerge.time(strat_beg_2_id);
    strat_end_2_time = wfpmerge.time(strat_end_2_id);
    
strat_beg_3_id = find(wfpmerge.time <= datenum(datetime(2017,6,1)) & wfpmerge.time >= datenum(datetime(2017,3,15)));
strat_end_3_id = find(wfpmerge.time <= datenum(datetime(2018,3,15)) & wfpmerge.time >= datenum(datetime(2017,12,1)));
oxygen_strat_beg_3 = oxygen_driftcorr_nooutliers_smoothed(:,strat_beg_3_id); %takes all O2 values within the specificed time (as the column) and all rows
oxygen_strat_end_3 = oxygen_driftcorr_nooutliers_smoothed(:,strat_end_3_id);
    strat_beg_3_time = wfpmerge.time(strat_beg_3_id);
    strat_end_3_time = wfpmerge.time(strat_end_3_id);
    
    
%% Implemented for 1st stratified season - separately find max and min O2 at beginning and end of stratified season
%Finding max and min O2 at each depth
  for j = 1: length(wfpmerge.depth_grid)
      [max_O2_season1(j), id_max_season1(j)] = max(oxygen_strat_beg_1(j,:)); %creates variables max_O2_season1 and id_max_season1 
      [min_O2_season1(j), id_min_season1(j)] = min(oxygen_strat_end_1(j,:));
  end
  
 for j = 1: length(wfpmerge.depth_grid)
      [max_O2_season2(j), id_max_season2(j)] = max(oxygen_strat_beg_2(j,:));
      [min_O2_season2(j), id_min_season2(j)] = min(oxygen_strat_end_2(j,:));
 end
 
 for j = 1: length(wfpmerge.depth_grid)
      [max_O2_season3(j), id_max_season3(j)] = max(oxygen_strat_beg_3(j,:));
      [min_O2_season3(j), id_min_season3(j)] = min(oxygen_strat_end_3(j,:));
 end
  
%Use the row indices of the maximum and minimum O2 at each depth to
%determine when those values occur
maxdate_O2_season1 = strat_beg_1_time(id_max_season1); %maxdate is a 491x1, so there is a dif maxdate for each depth 
mindate_O2_season1 = strat_end_1_time(id_min_season1);

maxdate_O2_season2 = strat_beg_2_time(id_max_season2);
mindate_O2_season2 = strat_end_2_time(id_min_season2);

maxdate_O2_season3 = strat_beg_3_time(id_max_season3);
mindate_O2_season3 = strat_end_3_time(id_min_season3);

%% Create array of linear fit values
for i = 1:length(depth_grid) %loop over all depths
    %Identify points at this depth that have non-nan O2 values
    idnonan = find(isnan(oxygen_driftcorr_nooutliers_smoothed(i,:)) == 0);
    %Stratified season 1
        %Pull out times within season
        time_id1 = find(wfpmerge.time > maxdate_O2_season1(i) & wfpmerge.time < mindate_O2_season1(i));
        %Find times within season that have non-nan O2 values
        time_id1_nonan = intersect(time_id1,idnonan);
        %Calculate slope and intercept of linear fit over smoothed oxygen
        %data within this stratified season
        [p1(i,:)] = polyfit(wfpmerge.time(time_id1_nonan), oxygen_driftcorr_nooutliers_smoothed(i,time_id1_nonan)',1);
    %Stratified season 2
        %Pull out times within season
        time_id2 = find(wfpmerge.time > maxdate_O2_season2(i) & wfpmerge.time < mindate_O2_season2(i));
        %Find times within season that have non-nan O2 values
        time_id2_nonan = intersect(time_id2,idnonan);
        %Calculate slope and intercept of linear fit over smoothed oxygen
        %data within this stratified season
        [p2(i,:)] = polyfit(wfpmerge.time(time_id2_nonan), oxygen_driftcorr_nooutliers_smoothed(i,time_id2_nonan)',1);
    %Stratified season 3
        %Pull out times within season
        time_id3 = find(wfpmerge.time > maxdate_O2_season3(i) & wfpmerge.time < mindate_O2_season3(i));
        %Find times within season that have non-nan O2 values
        time_id3_nonan = intersect(time_id3,idnonan);
        %Calculate slope and intercept of linear fit over smoothed oxygen
        %data within this stratified season
        [p3(i,:)] = polyfit(wfpmerge.time(time_id3_nonan), oxygen_driftcorr_nooutliers_smoothed(i,time_id3_nonan)',1);
end

%% Test plots to evaluate reasonableness of beginning/end dates and respiration rates for stratified seasons
    top_id = 31; bot_id = 241; %top and bottom depths to use in these test plots
    
%Show calculated dates for beginning and end of stratified season
figure(2); clf
plot(maxdate_O2_season1(top_id:bot_id), depth_grid(top_id:bot_id), 'ko','markerfacecolor','m'); hold on;
plot(mindate_O2_season1(top_id:bot_id), depth_grid(top_id:bot_id), 'm.');
plot(maxdate_O2_season2(top_id:bot_id), depth_grid(top_id:bot_id), 'ko','markerfacecolor','r');
plot(mindate_O2_season2(top_id:bot_id), depth_grid(top_id:bot_id), 'r.');
plot(maxdate_O2_season3(top_id:bot_id), depth_grid(top_id:bot_id), 'ko','markerfacecolor','g');
plot(mindate_O2_season3(top_id:bot_id), depth_grid(top_id:bot_id), 'g.');
datetick('x',2); set(gca,'YDir','reverse');
ylabel('Depth (m)'); title('Beginning (filled circles) and end (dots) dates of each stratified season')

%Show calculated respiration rates (slopes) for each stratified season
figure(3); clf
plot(p1(top_id:bot_id)*-365, depth_grid(top_id:bot_id), 'm.'); hold on;
plot(p2(top_id:bot_id)*-365, depth_grid(top_id:bot_id), 'r.');
plot(p3(top_id:bot_id)*-365, depth_grid(top_id:bot_id), 'g.');
set(gca,'YDir','reverse');
ylabel('Depth (m)')
xlabel('O_2 decrease rate (\mumol kg^{-1} yr^{-1})')
title('Respiration rates based on slope of linear fit over each stratified season')

%% %Makes figures at adjustable depths and plots max and min O2 to see if it matches with the smoothed data
for i=21:10:201 %101 is 650 meters, 41 is 350 meters
figure(i); clf
    depth_id = i; %Fixed to set depth_id to i within the loop
plot(wfpmerge.time, oxygen_driftcorr_smoothed(depth_id,:),'k.'); hold on;
plot(wfpmerge.time, oxygen_driftcorr_nooutliers_smoothed(depth_id,:),'m-'); hold on;

plot(maxdate_O2_season1(i), max_O2_season1(i), 'm.', 'MarkerSize', 25); hold on;
plot(maxdate_O2_season2(i), max_O2_season2(i), 'r.', 'MarkerSize', 25); hold on;
plot(maxdate_O2_season3(i), max_O2_season3(i), 'g.', 'MarkerSize', 25); hold on;

plot(mindate_O2_season1(i), min_O2_season1(i), 'm.', 'MarkerSize', 25); hold on;
plot(mindate_O2_season2(i), min_O2_season2(i), 'r.', 'MarkerSize', 25); hold on;
plot(mindate_O2_season3(i), min_O2_season3(i), 'g.', 'MarkerSize', 25); hold on;

time_id1 = find(wfpmerge.time > maxdate_O2_season1(i) & wfpmerge.time < mindate_O2_season1(i));
plot(wfpmerge.time(time_id1), wfpmerge.time(time_id1)*p1(i,1) + p1(i,2), 'b-','linewidth',2);
time_id2 = find(wfpmerge.time > maxdate_O2_season2(i) & wfpmerge.time < mindate_O2_season2(i));
plot(wfpmerge.time(time_id2), wfpmerge.time(time_id2)*p2(i,1) + p2(i,2), 'b-','linewidth',2);
time_id3 = find(wfpmerge.time > maxdate_O2_season3(i) & wfpmerge.time < mindate_O2_season3(i));
plot(wfpmerge.time(time_id3), wfpmerge.time(time_id3)*p3(i,1) + p3(i,2), 'b-','linewidth',2);

datetick('x',2)
%ylim([280 320])
title(['Oxygen time series at ' num2str(wfpmerge.depth_grid(i)) ' meters'])
end

%%
strat_season_length1 = (mindate_O2_season1 - maxdate_O2_season1);
strat_season_length2 = (mindate_O2_season2 - maxdate_O2_season2);
strat_season_length3 = (mindate_O2_season3 - maxdate_O2_season3);

%Calculate resp rates for each year
% resp_rate1 = (strat_season_length1/365).*p1(:,1);
% resp_rate2 = (strat_season_length2/365).*p2(:,1);
% resp_rate3 = (strat_season_length3/365).*p3(:,1);
% 
% sum_resp1 = sum(abs(resp_rate1))
% sum_resp2 = sum(abs(resp_rate2))
% sum_resp3 = sum(abs(resp_rate3))

resp_rate1 = (strat_season_length1).*-p1(:,1);
resp_rate2 = strat_season_length2.*-p2(:,1);
resp_rate3 = strat_season_length3.*-p3(:,1);
%maxdate-mindate and multiple by slopes 

%% Plot calculated respiration rates
figure(4); clf
plot(resp_rate1(top_id:bot_id), depth_grid(top_id:bot_id), 'm.'); hold on;
plot(resp_rate2(top_id:bot_id), depth_grid(top_id:bot_id), 'r.');
plot(resp_rate3(top_id:bot_id), depth_grid(top_id:bot_id), 'g.');
plot([0 0],[depth_grid(top_id) depth_grid(bot_id)], 'k--');
set(gca,'YDir','reverse');
ylabel('Depth (m)')
xlabel('O_2 decrease (\mumol kg^{-1})')
title('Total stratified season respiration based on slope of linear fit and length of stratified season')

%%
%Plot CHL at various depths 
for i=41:10:101 %101 is 650 meters, 41 is 350 meters
figure(i)
    i=depth_id; %Example plot at the 11th depth in depth_grid (which is 400 m), in depth_grid, 400m is 51st column
plot(wfpmerge.time, wfpmerge.chla(depth_id,:),'k.'); hold on;
ylim([0 0.3])
datetick('x',2)
end

%%
figure (1);
    top_id = 21;
subplot(1,3,1)
plot(max_O2_season1(:,top_id:211) - min_O2_season1(:,top_id:211), wfpmerge.depth_grid(:,top_id:211),'k.')
%datetick('x',2)
axis ij
xlabel ('O_2 decrease (\mumol kg^{-1})', 'FontSize', 14)
ylabel ('Depth (m)', 'FontSize', 14)
title('Year 1', 'FontSize', 16)

subplot(1,3,2)
plot(max_O2_season2(:,11:211) - min_O2_season2(:,11:211), wfpmerge.depth_grid(:,11:211),'k.')
%plot(max_O2_season2 - min_O2_season2, wfpmerge.depth_grid,'k.')
%datetick('x',2)
axis ij
xlabel ('O_2 decrease (\mumol kg^{-1})', 'FontSize', 14)
ylabel ('Depth (m)', 'FontSize', 14)
title('Year 2', 'FontSize', 16)


subplot(1,3,3)
plot(max_O2_season3(:,51:211) - min_O2_season3(:,51:211), wfpmerge.depth_grid(:,51:211),'k.')
%plot(max_O2_season3 - min_O2_season3, wfpmerge.depth_grid,'k.')
%datetick('x',2)
axis ij
xlabel ('O_2 decrease (\mumol kg^{-1})', 'FontSize', 14)
ylabel ('Depth (m)', 'FontSize', 14)
title('Year 3', 'FontSize', 16)

%%
%now calculate respiration rates throughout water column by integrating. From Max to Min at each depth makes a curve, then find area under curve 
%% Alternative version of figure above calculated from slopes
    top_id = 15;
figure (10);
subplot(1,3,1)
plot(strat_season_length1(top_id:bot_id).*-p1(top_id:bot_id,1), wfpmerge.depth_grid(top_id:bot_id),'k.')
%datetick('x',2)
axis ij
xlabel ('O_2 decrease (\mumol kg^{-1})', 'FontSize', 14)
ylabel ('Depth (m)', 'FontSize', 14)
title('Year 1', 'FontSize', 16)

subplot(1,3,2)
plot(strat_season_length2(top_id:bot_id).*-p2(top_id:bot_id,1), wfpmerge.depth_grid(top_id:bot_id),'k.')
%plot(max_O2_season2 - min_O2_season2, wfpmerge.depth_grid,'k.')
%datetick('x',2)
axis ij
xlabel ('O_2 decrease (\mumol kg^{-1})', 'FontSize', 14)
ylabel ('Depth (m)', 'FontSize', 14)
title('Year 2', 'FontSize', 16)


subplot(1,3,3)
plot(strat_season_length3(top_id:bot_id).*-p3(top_id:bot_id,1), wfpmerge.depth_grid(top_id:bot_id),'k.')
%plot(max_O2_season3 - min_O2_season3, wfpmerge.depth_grid,'k.')
%datetick('x',2)
axis ij
xlabel ('O_2 decrease (\mumol kg^{-1})', 'FontSize', 14)
ylabel ('Depth (m)', 'FontSize', 14)
title('Year 3', 'FontSize', 16)

%% Integrate to calculate full stratified season respiration rate
%Calculate respiration rate in each depth interval:
   %Total resp rate (in umol O2/kg) x density at each depth interval
   %(kg/m3) x interval between depths (m) x (conversion from umol to mol) =
   %rate for each depth bin in mol/m2
ThermResp1 = (max_O2_season1 - min_O2_season1)'.*nanmean(Yr1_wfpgrid.pdens,2).*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000; %respiration per depth interval (mol/m2)
ThermResp1_slope = resp_rate1.*nanmean(Yr1_wfpgrid.pdens,2).*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000; %respiration per depth interval (mol/m2)
%Specify how deep to integrate (you could use similar approach to specify where to start from at top)    
    topdepth = 240;
    id_topdepth = find(wfpmerge.depth_grid == topdepth);
    intdepth = 1000; %choose integration depth for base of seasonal thermocline
    id_basetherm = find(wfpmerge.depth_grid == intdepth);
%Take integral of ThermResp from toptherm to basetherm using rectangle sum method
ThermResp_Int1 = nansum(ThermResp1(id_topdepth:id_basetherm)); %mol O2 m-2 respired and ventilated during winter
ThermResp_Int1_slope = nansum(ThermResp1_slope(id_topdepth:id_basetherm)); 
%Calculate cumulative integral from top therm to basetherm using trapezoid
%method (note that you need to choose topdepth to avoid NaNs with this
%approach)
ThermResp_CumInt1 = cumtrapz(ThermResp1(id_topdepth:id_basetherm));

%%
%Year 2
ThermResp2 = (max_O2_season2 - min_O2_season2)'.*nanmean(Yr2_wfpgrid.pdens,2).*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000; %respiration per depth interval (mol/m2)
ThermResp2_slope = resp_rate2.*nanmean(Yr2_wfpgrid.pdens,2).*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000; %respiration per depth interval (mol/m2)
%Specify how deep to integrate (you could use similar approach to specify where to start from at top)    
    topdepth = 240;
    id_topdepth = find(wfpmerge.depth_grid == topdepth);
    intdepth = 1100; %choose integration depth for base of seasonal thermocline
    id_basetherm = find(wfpmerge.depth_grid == intdepth);
%Take integral of ThermResp from toptherm to basetherm using rectangle sum method
ThermResp_Int2 = nansum(ThermResp2(id_topdepth:id_basetherm)); %mol O2 m-2 respired and ventilated during winter
ThermResp_Int2_slope = nansum(ThermResp2_slope(id_topdepth:id_basetherm)); 
%Calculate cumulative integral from top therm to basetherm using trapezoid
%method (note that you need to choose topdepth to avoid NaNs with this
%approach)
ThermResp_CumInt2 = cumtrapz(ThermResp2(id_topdepth:id_basetherm));

%%
%Year 3
ThermResp3 = (max_O2_season3 - min_O2_season3)'.*nanmean(Yr3_wfpgrid.pdens,2).*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000; %respiration per depth interval (mol/m2)
ThermResp3_slope = resp_rate3.*nanmean(Yr3_wfpgrid.pdens,2).*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000;
%Specify how deep to integrate (you could use similar approach to specify where to start from at top)    
    topdepth = 400;
    id_topdepth = find(wfpmerge.depth_grid == topdepth);
    intdepth = 1300; %choose integration depth for base of seasonal thermocline
    id_basetherm = find(wfpmerge.depth_grid == intdepth);
%Take integral of ThermResp from toptherm to basetherm using rectangle sum method
ThermResp_Int3 = nansum(ThermResp3(id_topdepth:id_basetherm)); %mol O2 m-2 respired and ventilated during winter
ThermResp_Int3_slope = nansum(ThermResp2_slope(id_topdepth:id_basetherm)); 
%Calculate cumulative integral from top therm to basetherm using trapezoid
%method (note that you need to choose topdepth to avoid NaNs with this
%approach)
ThermResp_CumInt3 = cumtrapz(ThermResp3(id_topdepth:id_basetherm));

%%  
% %%
%   for j = 1: length(wfpmerge.depth_grid)
%       max_O2_season2(j) = max(oxygen_strat_season_2(j,:));
%       min_O2_season2(j) = min(oxygen_strat_season_2(j,:));  
%   end
%   
% %%
%   for j = 1: length(wfpmerge.depth_grid)
%       max_O2_season3(j) = max(oxygen_strat_season_3(j,:));
%       min_O2_season3(j) = min(oxygen_strat_season_3(j,:));
%   end
  
%%
%Strat season is different for each depth in water column, so need to write something that identifies the strat season for each depth
% for j = 1: length(wfpmerge.depth_grid)) 
%     strat_season_1(j) = max (
      
%%      
% wfpmerge.time_smoothed = movmean(wfpmerge.time, 50);
% plot(wfpmerge.time_smoothed, oxygen_driftcorr_smoothed, 'm-'); %having
% trouble with getting time to show up as real dates and not numbers--
% fixed that using code below

% id_maxDO= max(oxygen_driftcorr_smoothed)'; %want to find max and min DO at everydepth over stratified season --what is in each column of DO -- this function took the max and min of each column
% id_minDO= min(oxygen_driftcorr_smoothed)';
%%
% figure(1); clf;
% [X,Y] = meshgrid(wfpmerge.time, wfpmerge.depth_grid);
% contourf(X,Y,oxygen_driftcorr_smoothed,cvec,'linecolor','none'); hold on;
% axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Oxygen concentration (mol/L)', 'Fontsize', 12)
%%
% figure(2); clf;
% [X,Y] = meshgrid(Trans_wfpmerge.time, wfpmerge.depth_grid);
% contourf(X,Y,id_maxDO-id_minDO,cvec,'linecolor','none'); hold on;
% axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Oxygen concentration (mol/L)', 'Fontsize', 12)
% Trans_wfpmerge.time = wfpmerge.time';

% figure(2);
% plot (wfpmerge.time, id_maxDO-id_minDO, 'm.'); hold on; 
% datetick('x',2,'keeplimits');
% %plot (wfpmerge.time, id_maxDO-id_minDO);
% %fitline = polyfit(wfpmerge.time, id_maxDO-id_minDO, 1);
% axis([wfpmerge.time wfpmerge.time]);


%strat_season1 = find(wfpmerge.time <= november 2016 & stationdata.Year >= june 2016) %figure out how to put in time that matlab recognizes
%then within this season, identify max and min O2 at each depth
        