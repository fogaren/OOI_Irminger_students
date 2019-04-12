
%Plot merged data
%Adjustable parameters for plotting
    mindepth = 150; maxdepth = 1700;
    cints = 60; %number of contour intervals
    C = cmocean('Dense'); %set colormap
    C2 = cmocean('Algae'); 

%Make plotting grid
[X,Y] = meshgrid(wfpmerge.time, wfpmerge.depth_grid);

figure(5);clf;
%subplot(311) %Oxygen_corr concentration
cmin = 260; cmax = 320; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.oxygen_corr,cvec,'linecolor','none'); hold on; %(11:271,:) for the Z variable?
%contourf(X,Y,max_O2_season1,cvec,'r.'); hold on;
%plot(strat_beg_1_time(id_max_season1(1,15:1,231)), wfpmerge.depth_grid(1,15:1,231), 'r.'); hold on; %want the y axis to be just the location of the max or min O2
plot(strat_beg_1_time(id_max_season1(23:191)), wfpmerge.depth_grid(23:191), 'r.'); hold on; %1200 m
plot(strat_end_1_time(id_min_season1(11:191)), wfpmerge.depth_grid(11:191), 'y.'); hold on; %1200 m
plot(strat_beg_2_time(id_max_season2(11:191)), wfpmerge.depth_grid(11:191), 'r.'); hold on;
plot(strat_end_2_time(id_min_season2(11:131)), wfpmerge.depth_grid(11:131), 'y.'); hold on;
plot(strat_beg_3_time(id_max_season3(11:131)), wfpmerge.depth_grid(11:131), 'r.'); hold on; %800 m
plot(strat_end_3_time(id_min_season3(51:131)), wfpmerge.depth_grid(51:131), 'y.'); hold on; %400 to 800 m
%plot(mindate_O2_season1, min_O2_season1, 'r.');
%scatter(maxdate_O2_season1, max_O2_season1,25,'r','*')
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 15); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration (mol/L)', 'Fontsize', 17)

%%

% figure(1);clf;
% %subplot(311) %Oxygen_corr concentration
% cmin = 260; cmax = 320; %manually set min and max
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
% contourf(X,Y,wfpmerge.oxygen_corr,cvec,'linecolor','none'); hold on;
% %contourf(X,Y,max_O2_season1,cvec,'r.'); hold on;
% plot(strat_end_1_time(id_min_season1), wfpmerge.depth_grid, 'r.'); hold on; 
% plot(strat_beg_1_time(id_max_season1), wfpmerge.depth_grid, 'y.'); hold on;
% plot(strat_end_2_time(id_min_season2), wfpmerge.depth_grid, 'r.'); hold on;
% plot(strat_beg_2_time(id_max_season2), wfpmerge.depth_grid, 'y.'); hold on;
% plot(strat_end_3_time(id_min_season3), wfpmerge.depth_grid, 'r.'); hold on;
% plot(strat_beg_3_time(id_max_season3), wfpmerge.depth_grid, 'y.'); hold on;
% %plot(mindate_O2_season1, min_O2_season1, 'r.');
% %scatter(maxdate_O2_season1, max_O2_season1,25,'r','*')
% axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
% colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
% datetick('x',2,'keeplimits');
% title('Oxygen concentration (mol/L)', 'Fontsize', 12)

figure(2);clf; %Oxygen saturation
cmin = -15; cmax = 0; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.O2satcorr,cvec,'linecolor','none'); hold on;
plot(strat_beg_1_time(id_max_season1), wfpmerge.depth_grid, 'r.'); hold on;
plot(strat_end_1_time(id_min_season1), wfpmerge.depth_grid, 'y.'); hold on;
plot(strat_beg_2_time(id_max_season2), wfpmerge.depth_grid, 'r.'); hold on;
plot(strat_end_2_time(id_min_season2), wfpmerge.depth_grid, 'y.'); hold on;
plot(strat_beg_3_time(id_max_season3), wfpmerge.depth_grid, 'r.'); hold on;
plot(strat_end_3_time(id_min_season3), wfpmerge.depth_grid, 'y.'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen saturation (%)', 'Fontsize', 12)

figure(3);clf; %Chlorophyll
cmin = 0; cmax = 0.3; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.chla,cvec,'linecolor','none'); hold on;
plot(strat_beg_1_time(id_max_season1), wfpmerge.depth_grid, 'r.'); hold on;
plot(strat_end_1_time(id_min_season1), wfpmerge.depth_grid, 'y.'); hold on;
plot(strat_beg_2_time(id_max_season2), wfpmerge.depth_grid, 'r.'); hold on;
plot(strat_end_2_time(id_min_season2), wfpmerge.depth_grid, 'y.'); hold on;
plot(strat_beg_3_time(id_max_season3), wfpmerge.depth_grid, 'r.'); hold on;
plot(strat_end_3_time(id_min_season3), wfpmerge.depth_grid, 'y.'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)'); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Chlorophyll (µg/L)', 'Fontsize', 15)