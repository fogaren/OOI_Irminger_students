    %% Intial Plotting
    %Surface Bouy
       % Surface Bouy Plots that include the dissolved oxygen, Temperature (C °),
       % salinity and exploration into what the heck ct_depth is. 
       figure(1); clf;
       subplot(2,2,1) 
       plot(ApexS.Time, ApexS.O2, "k.")
       datetick('x',12);
       ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
       title('Apex Surface Year 3 O_2')
       
       subplot(2,2,2) 
       plot(ApexS.Time, ApexS.surfacetemp, "k.")
       datetick('x',12);
       ylabel('Temperature (C °)','Fontsize',10)
       title('Apex Surface Year 3 Temperature (C °)')
       
       subplot(2,2,3) 
       plot(ApexS.Time, ApexS.sal, "k.")
       datetick('x',12);
       ylabel('Practical Salinity At Sea Surface','Fontsize',10)
       title('Apex Surface Year 3 Salinity')
       
       subplot(2,2,4) 
       plot(ApexS.Time, ApexS.depth, "k.")
       datetick('x',12);
       title('Apex Surface Year 3 Depth of Instrument')
       
       %% Near Surface Bouy
       
       figure(2); clf;
       subplot(2,2,1) 
       plot(NearS.Time, NearS.O2, "k.")
       datetick('x',12);
       ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
       title('Near Surface Year 3 O_2')
       
       subplot(2,2,2) 
       plot(NearS.Time, NearS.temp, "k.")
       datetick('x',12);
       ylabel('Temperature (C °) (C°)','Fontsize',10)
       title('Near Surface Year 3 Temperature (C°)')
       
       subplot(2,2,3) 
       plot(NearS.Time, NearS.sal, "k.")
       datetick('x',12);
       ylabel('Practical Salinity At Sea Surface','Fontsize',10)
       title('Apex Surface Year 3 Salinity')
       
       subplot(2,2,4) 
       plot(NearS.Time, NearS.p, "k.")
       datetick('x',12);
       title('Apex Surface Year 3 Pressure')
       
       
       %% Investigating Temperature (C °) anomaly Jan in Near Surface
       
        % Plot with Flanking Mooring %%
       figure(3); clf;
       plot(NearS.Time, NearS.temp, "b."); hold on;
       plot(ApexS.Time, ApexS.surfacetemp, "k."); hold on; 
       plot(FlankingB_strt.time, FlankingB_strt.temp, "g."); hold on;
       plot(FlankingA_strt.Time, FlankingA_strt.temp, "m.");
       datetick('x',12)
       ylabel('Temperature (C °) in Y 3','Fontsize',10)
       leg = legend('Near Surface', 'Surface', 'Flanking B', 'Flanking A','location', 'eastoutside')
       set(leg,'FontSize',14)
       title('Investigating Temperature (C °) anomaly!')
       
       %% Flanking A
       figure(4); clf;
       subplot(2,2,1) 
       plot(FlankingA_strt.Time, FlankingA_strt.O2, "k.")
       datetick('x',12);
       ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
       title('Flanking A Year 3 O_2')
       
       subplot(2,2,2) 
       plot(FlankingA_strt.Time, FlankingA_strt.temp, "k.")
       datetick('x',12);
       ylabel('Temperature(C°)','Fontsize',10)
       title('FlankingA Year 3 Temperature (C °)')
       
       subplot(2,2,3) 
       plot(FlankingA_strt.Time, FlankingA_strt.sal, "k.")
       datetick('x',12);
       ylabel('Practical Salinity At Sea Surface','Fontsize',10)
       title('Flanking A Year 3 Salinity')
       
       subplot(2,2,4) 
       plot(FlankingA_strt.Time, FlankingA_strt.p, "k.")
       datetick('x',12);
       title('Flanking A Year 3 Pressure')
    
       
       %% Flanking B
       figure(5); clf;
       subplot(2,2,1) 
       plot(FlankingB_strt.time, FlankingB_strt.O2, "k.")
       datetick('x',12);
       ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
       title('Flanking A Year 3 O_2')
       
       subplot(2,2,2) 
       plot(FlankingB_strt.time, FlankingB_strt.temp, "k.")
       datetick('x',12);
       ylabel('Temperature (C °)','Fontsize',10)
       title('FlankingA Year 3 Temperature (C °)')
       
       subplot(2,2,3) 
       plot(FlankingB_strt.time, FlankingB_strt.sal, "k.")
       datetick('x',12);
       ylabel('Practical Salinity At Sea Surface','Fontsize',10)
       title('Flanking A Year 3 Salinity')
       
       subplot(2,2,4) 
       plot(FlankingB_strt.time, FlankingB_strt.p, "k.")
       datetick('x',12);
       title('Flanking A Year 3 Pressure')
    
       
       %% Some Comparison Graphs 
        % Dissolved Oxygen as measured by all sensors, code order :
        % Surface, Near Surface, Flanking A, Flanking B
       figure(6); clf;
       subplot(1,2,1) %This is now just the Apex Mooring 
    plot(ApexS.Time,ApexS.O2,'k.'); hold on;
    plot(NearS.Time,NearS.O2,'g.'); hold on;
    %plot(FlankingA.Time, FlankingA.O2, 'm.'); hold on;
    %plot(FlankingB.time, FlankingB.02, 'b.');
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) 
    ylim([200 400]) 
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Dissolved Oxygen as Measured by All Sensors over the year')
    legend('Apex Surface', 'Apex Near Surface')%, 'Flanking A', 'Flanking B')
    
       subplot(1,2,2) %This is the data from the flanking mooringss
    plot(FlankingA_strt.Time, FlankingA_strt.O2, 'm.'); hold on;
    plot(FlankingB_strt.time, FlankingB_strt.O2, 'b.')
    %xlim([datenum(2015,9,1) datenum(2016,8,1)]) 
    ylim([200 400]) 
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Dissolved Oxygen as Measured by All Sensors over the year')
    legend('Flanking A', 'Flanking B')
    
    %% Observed V. Expected in equilibrium
    % Expected values as calculated from GSW O2 function using salinity and
    % Temperature (C °) 
    
    figure (7); clf; 
    subplot(2,2,1)
    %Apex SurfaceMooring
    %This plot was also used to compare the dosta oxygen readings, this
    %part of the plot has been commented out. 
    %plot(ApexS.Time, ApexS.O2_dosta, 'r.'); hold on;
    plot(ApexS.Time, ApexS.O2, 'k.');hold on;
    plot(ApexS.Time, ApexS.O2_expected, 'b.')
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Expected O_2 as compared to observed O_2') %title('Comparing Dosta Oxygen, Adjusted Oxygen, and Expected Oxygen')
    %legend('Dosta O_2', 'DO O_2', 'Expected O_2','location','eastoutside')% 'FontSize', 13) 
    %Dosta Oxygen is not adjusted, 'dissolved oxygen' (named O2 in this
    %script) is the correct reading to use 
    
    %Near Surface Mooring 
    subplot(2,2,2)
    plot(NearS.Time, NearS.O2, 'k.');hold on;
    plot(NearS.Time, NearS.O2_expected, 'b.')
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Expected O_2 as compared to observed O_2; Near Surface')
    legend('Observed O_2', 'Expected O_2', 'location', 'eastoutside')
    
    
    subplot(2,2,3)
    plot(FlankingA_strt.Time, FlankingA_strt.O2, 'k.');hold on;
    plot(FlankingA_strt.Time, FlankingA_strt.O2_expected, 'b.')
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Expected O_2 as compared to observed O_2; Flanking A')
    %legend('Observed O_2', 'Expected O_2', 'location', 'eastoutside')
    
    subplot(2,2,4)
    plot(FlankingB_strt.time, FlankingB_strt.O2, 'k.');hold on;
    plot(FlankingB_strt.time, FlankingB_strt.O2_expected, 'b.')
    datetick('x',12);
    ylabel('Dissolved oxygen (\mumol kg^{-1})','Fontsize',10)
    title('Expected O_2 as compared to observed O_2; Flanking B')
    legend('Observed O_2', 'Expected O_2', 'location', 'eastoutside')