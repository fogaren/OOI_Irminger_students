%dissolved oxygen time series at certain depths to show respiration

[X,Y] = meshgrid(wfpmerge.time, wfpmerge.depth_grid);
Z = find(wfpmerge.oxygen_corr(11,:)>260);

%D1 = find(wfpmerge.time == );
figure(1); clf;
scatter(wfpmerge.time(Z),wfpmerge.oxygen_driftcorr(11,Z));
datetick('x',2,'keeplimits');