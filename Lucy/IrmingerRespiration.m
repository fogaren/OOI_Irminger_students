
%plot(wfpmerge.time,wfpmerge.oxygen_driftcorr)

oxygen_driftcorr_smoothed = movmean(wfpmerge.oxygen_driftcorr, 200); %unsure at what interval I want to smooth oxygen data
% wfpmerge.time_smoothed = movmean(wfpmerge.time, 50);
% plot(wfpmerge.time_smoothed, oxygen_driftcorr_smoothed, 'm-'); %having
% trouble with getting time to show up as real dates and not numbers--
% fixed that using code below

id_maxDO= max(oxygen_driftcorr_smoothed)'; %want to find max and min DO at everydepth over stratified season --what is in each column of DO -- this function took the max and min of each column
id_minDO= min(oxygen_driftcorr_smoothed)';
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
plot (wfpmerge.time, id_maxDO-id_minDO, 'm.'); hold on; 
datetick('x',2,'keeplimits');
%plot (wfpmerge.time, id_maxDO-id_minDO);
%fitline = polyfit(wfpmerge.time, id_maxDO-id_minDO, 1);
axis([wfpmerge.time wfpmerge.time]);

%next step is to identify the strat season and then use the find function to pick out only the strat season dates (but not sure how to do this in a way that matlab recognizes)  

%strat_season1 = find(wfpmerge.time <= november 2016 & stationdata.Year >= june 2016) %figure out how to put in time that matlab recognizes
%then within this season, identify max and min O2 at each depth




        
        
    