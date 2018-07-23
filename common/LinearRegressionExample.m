%%% Example code for plotting a linear regression

    ind = find(ssw.SSS_filt > 30 & gps.lon_filt < -25); %remove outliers, since they will cause issues for the regression
P = polyfit(ssw.SST_filt(ind), ssw.SSS_filt(ind), 1); %calculate slope and intercept (the "1" is for 1st degree polynomial, i.e. linear fit)
    XtoPlot = [4:10]; %set range of x values to plot over
figure; clf
plot(ssw.SST_filt(ind), ssw.SSS_filt(ind), 'k.'); hold on;
plot(XtoPlot, P(1)*XtoPlot + P(2), 'r--'); hold on;

[rho,df,rho_sig95] = correlate(ssw.SST_filt(ind),ssw.SSS_filt(ind)); %calculate correlation statistics