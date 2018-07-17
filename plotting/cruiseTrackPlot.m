%This will use the m_map toolbox to plot the Irminger-5 cruise track and locations of
%the OOI moorings

%The latitude and longitude for plotting come from having run the
%underway_DataAnalysis script

%The OOI mooring locations (from cruise report anchor surveys) are saved as
% OOImooringLocations.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load OOImooringLocations.mat

figure; clf
    M = 15;
latminplot = 58; latmaxplot = 67; lonminplot = -45; lonmaxplot = -15;
m_proj('miller','lat',[latminplot latmaxplot],'lon',[lonminplot lonmaxplot])
[CS,CH]=m_tbase('contourf',[-5000:200:0],'edgecolor','none');
h = colorbar; colormap(flipud(cmocean('deep')));
m_grid('box','fancy'); m_coast('patch',nicecolor('wwk'));
m_plot(gps.lon_filt, gps.lat_filt,'r-'); hold on;
m_plot(OOImoorings.SUMO4(2), OOImoorings.SUMO4(1),'m.','markersize',M); hold on;
m_plot(OOImoorings.HYPM4(2), OOImoorings.HYPM4(1),'m.','markersize',M); hold on;
m_plot(OOImoorings.FLMA4(2), OOImoorings.FLMA4(1),'m.','markersize',M); hold on;
m_plot(OOImoorings.FLMB4(2), OOImoorings.FLMB4(1),'m.','markersize',M); hold on;
ylabel(h,'Depth (m)')