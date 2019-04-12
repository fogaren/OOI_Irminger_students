%% Correct O2 data years 1-4
    %makes an array of times profiles were collected, the same size as the paired profile data
    %note that this is a throwaway variable that is re-written for each year (though you could change and save it in the structure)
    timegrid = repmat(Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair),1,length(depth_grid))';
Yr1_wfpgrid.oxygen_driftcalib = (timegrid - Yr1_wfpgrid.time_start(1)) * wfp_O2drift(1,1); %calculates cumulative drift since deployment = (today's date - day 1 of deployment)* (drift rate of O2 per day)
    timegrid = repmat(Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair),1,length(depth_grid))';
Yr2_wfpgrid.oxygen_driftcalib = (timegrid - Yr2_wfpgrid.time_start(1)) * wfp_O2drift(1,2);
    timegrid = repmat(Yr3_wfpgrid.time_start(Yr3_wfpgrid.ind_pair),1,length(depth_grid))';
Yr3_wfpgrid.oxygen_driftcalib = (timegrid - Yr3_wfpgrid.time_start(1)) * wfp_O2drift(1,3);
    timegrid = repmat(Yr4_wfpgrid.time_start(Yr4_wfpgrid.ind_pair),1,length(depth_grid))';
Yr4_wfpgrid.oxygen_driftcalib = (timegrid - Yr4_wfpgrid.time_start(1)) * wfp_O2drift(1,4);

%Note the sign of the drift above is negative (sensor drifts down), so
%subtract the driftcalib value to correct the value back up to correct number
Yr1_wfpgrid.oxygen_driftcorr = Yr1_wfpgrid.oxygen_corr - Yr1_wfpgrid.oxygen_driftcalib; %basically subtracting the drift calculation
Yr2_wfpgrid.oxygen_driftcorr = Yr2_wfpgrid.oxygen_corr - Yr2_wfpgrid.oxygen_driftcalib;
Yr3_wfpgrid.oxygen_driftcorr = Yr3_wfpgrid.oxygen_corr - Yr3_wfpgrid.oxygen_driftcalib;
Yr4_wfpgrid.oxygen_driftcorr = Yr4_wfpgrid.oxygen_corr - Yr4_wfpgrid.oxygen_driftcalib;

%%
%Calculate O2 saturation with drift correction   
Yr1_wfpgrid.O2sat_driftcorr = (Yr1_wfpgrid.oxygen_driftcorr./Yr1_wfpgrid.O2equil - 1)*100;
Yr2_wfpgrid.O2sat_driftcorr = (Yr2_wfpgrid.oxygen_driftcorr./Yr2_wfpgrid.O2equil - 1)*100;
Yr3_wfpgrid.O2sat_driftcorr = (Yr3_wfpgrid.oxygen_driftcorr./Yr3_wfpgrid.O2equil - 1)*100;
Yr4_wfpgrid.O2sat_driftcorr = (Yr4_wfpgrid.oxygen_driftcorr./Yr4_wfpgrid.O2equil - 1)*100;

%%
wfpmerge.oxygen_driftcorr = [Yr1_wfpgrid.oxygen_driftcorr Yr2_wfpgrid.oxygen_driftcorr Yr3_wfpgrid.oxygen_driftcorr Yr4_wfpgrid.oxygen_driftcorr];
wfpmerge.O2sat_driftcorr = [Yr1_wfpgrid.O2sat_driftcorr Yr2_wfpgrid.O2sat_driftcorr Yr3_wfpgrid.O2sat_driftcorr Yr4_wfpgrid.O2sat_driftcorr];
%% Plot merged data
%Adjustable parameters for plotting
    mindepth = 150; maxdepth = 2600;
    cints = 60; %number of contour intervals
    C = cmocean('Dense'); %set colormap
    C2 = cmocean('Algae'); 

%Make plotting grid
[X,Y] = meshgrid(wfpmerge.time, wfpmerge.depth_grid);

figure(i + 1); clf;
    subplot(411) %Density
cmin = 27.5; cmax = 27.8; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.pdens - 1000,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('\sigma_\theta', 'Fontsize', 12)

    subplot(412) %Temperature
cmin = 2; cmax = 6; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.T,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Temperature (deg C)', 'Fontsize', 12)

    subplot(413) %Oxygen_corr concentration
cmin = 260; cmax = 320; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.oxygen_driftcorr,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration (mol/L)', 'Fontsize', 12)

    subplot(414) %Oxygen saturation
cmin = -15; cmax = 0; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.O2sat_driftcorr,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen saturation (%)', 'Fontsize', 12)

figure (i + 2); clf;
    subplot(211) %Backscatter
cmin = 4E-4; cmax = 1.5E-3; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.backscatter,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)'); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Backscatter', 'Fontsize', 15)

    subplot(212) %Chlorophyll
cmin = 0; cmax = 0.3; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.chla,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)'); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Chlorophyll (µg/L)', 'Fontsize', 15)
