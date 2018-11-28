%% Correct O2 data years 1-4
%Initial gain correction of DOSTA    
Yr1_wfp.oxygen_corr = Yr1_wfp.oxygen * gain(1,1); 
Yr2_wfp.oxygen_corr = Yr2_wfp.oxygen * gain(1,2);
Yr3_wfp.oxygen_corr = Yr3_wfp.oxygen * gain(1,3);
Yr4_wfp.oxygen_corr = Yr4_wfp.oxygen * gain(1,4);

Yr1_wfpgrid.oxygen_corr = Yr1_wfpgrid.O2conc * gain(1,1); 
Yr2_wfpgrid.oxygen_corr = Yr2_wfpgrid.O2conc * gain(1,2);
Yr3_wfpgrid.oxygen_corr = Yr3_wfpgrid.O2conc * gain(1,3);
Yr4_wfpgrid.oxygen_corr = Yr4_wfpgrid.O2conc * gain(1,4);

Yr1_wfpgrid_therm.oxygen_corr = Yr1_wfpgrid_therm.O2conc * gain(1,1); 
Yr2_wfpgrid_therm.oxygen_corr = Yr2_wfpgrid_therm.O2conc * gain(1,2);
Yr3_wfpgrid_therm.oxygen_corr = Yr3_wfpgrid_therm.O2conc * gain(1,3);
Yr4_wfpgrid_therm.oxygen_corr = Yr4_wfpgrid_therm.O2conc * gain(1,4);

%%
%Calculate O2 saturation with initial correction   
Yr1_wfpgrid.O2satcorr = (Yr1_wfpgrid.oxygen_corr./Yr1_wfpgrid.O2equil - 1)*100;
Yr2_wfpgrid.O2satcorr = (Yr2_wfpgrid.oxygen_corr./Yr2_wfpgrid.O2equil - 1)*100;
Yr3_wfpgrid.O2satcorr = (Yr3_wfpgrid.oxygen_corr./Yr3_wfpgrid.O2equil - 1)*100;
Yr4_wfpgrid.O2satcorr = (Yr4_wfpgrid.oxygen_corr./Yr4_wfpgrid.O2equil - 1)*100;

%% Visualize initial gain correction data
for i = 1:4
    if i == 1
        plotting = Yr1_wfpgrid;
    elseif i == 2
        plotting = Yr2_wfpgrid;
    elseif i == 3
        plotting = Yr3_wfpgrid;
    elseif i == 4
        plotting = Yr4_wfpgrid;
    end


figure(i); clf
set(gcf,'color','w')
x0=1;
y0=1;
width=28;
height=20;
set(gcf,'units','centimeters','position',[x0,y0,width,height]) 
    subplot(4,2,1)
imagesc(plotting.T); colorbar; title('Temperature');
    subplot(4,2,2)
imagesc(plotting.S); colorbar; title('Salinity');
    subplot(4,2,3)
imagesc(plotting.pdens); colorbar; title('Density');
    subplot(4,2,4)
imagesc(plotting.oxygen_corr); colorbar; caxis([240 300]); title('O_2 concentration');
    subplot(4,2,5)
imagesc(plotting.O2sat); colorbar; caxis([-25 0]); title('O_2 sat');
    subplot(4,2,6)
imagesc(plotting.backscatter); colorbar; caxis([0 0.002]); title('Backscatter');
    subplot(4,2,7)
imagesc(plotting.scat_total); colorbar; title('Scat Total');
    subplot(4,2,8)
imagesc(plotting.chla); colorbar; caxis([0 0.3]); title('Chlorophyll a');

end

%%
wfpmerge.oxygen_corr = [Yr1_wfpgrid.oxygen_corr Yr2_wfpgrid.oxygen_corr Yr3_wfpgrid.oxygen_corr Yr4_wfpgrid.oxygen_corr];
wfpmerge.O2satcorr = [Yr1_wfpgrid.O2satcorr Yr2_wfpgrid.O2satcorr Yr3_wfpgrid.O2satcorr Yr4_wfpgrid.O2satcorr];
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
contourf(X,Y,wfpmerge.oxygen_corr,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C); set(gca,'YDir','reverse'); ylabel('Depth (m)', 'Fontsize', 10); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Oxygen concentration (mol/L)', 'Fontsize', 12)

    subplot(414) %Oxygen saturation
cmin = -15; cmax = 0; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.O2satcorr,cvec,'linecolor','none'); hold on;
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
title('Backscatter')

    subplot(212) %Chlorophyll
cmin = 0; cmax = 0.3; %manually set min and max
    cvec = [cmin:(cmax-cmin)/cints:cmax];
contourf(X,Y,wfpmerge.chla,cvec,'linecolor','none'); hold on;
axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
colormap(C2); set(gca,'YDir','reverse'); ylabel('Depth (m)'); hcb = colorbar; set(hcb,'location','eastoutside')
datetick('x',2,'keeplimits');
title('Chlorophyll (µg/L)', 'Fontsize', 15)
