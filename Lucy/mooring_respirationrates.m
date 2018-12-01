%dissolved oxygen time series at certain depths to show respiration

depth_to_plot = 200; %set depth to plot
    [~, depth_id] = min(abs(depth_grid - depth_to_plot)); %find index in depth_grid corresponding to depth_to_plot
minval = 270; %set minimum value to include in plot (to exclude low outliers)

%[X,Y] = meshgrid(wfpmerge.time, wfpmerge.depth_grid);
Z = find(wfpmerge.oxygen_corr(depth_id,:) > minval);

figure(1); clf;
scatter(wfpmerge.time(Z),wfpmerge.oxygen_driftcorr(depth_id,Z));
datetick('x',2,'keeplimits');