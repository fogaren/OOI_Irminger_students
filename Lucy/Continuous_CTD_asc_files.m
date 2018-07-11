%Script for importing excel data from ASC files, which contain continous
%CTD data from each cast
path = ['/Users/lucywanzer/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/'];
file = ['ASC_Continuous_CTD_Data.xlsx'];
filename= [path file]
castnums = [1:8,10:12,14:22];
for i = 1:length(castnums)
    ASCContinuousCTDData = xlsread(filename,i);
    cast{castnums(i)}.Pres= ASCContinuousCTDData (:,1);
    cast{castnums(i)}.T1= ASCContinuousCTDData (:,2); %T090C Temperature (ITS-90, deg C)
    cast{castnums(i)}.T2= ASCContinuousCTDData (:,3); %T190C Temperature (ITS-90, deg C)
    cast{castnums(i)}.C1= ASCContinuousCTDData (:,4); %COS/m Conductivity (S/m)
    cast{castnums(i)}.C2= ASCContinuousCTDData (:,5); %C1S/m Conductivity (S/m)
    cast{castnums(i)}.OR= ASCContinuousCTDData (:,6); %sbeox0V: Oxygen raw, SBE 43 [V]
    cast{castnums(i)}.F = ASCContinuousCTDData (:,7); %flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]
    cast{castnums(i)}.Turb = ASCContinuousCTDData (:,8); %turbWETntu0: Turbidity, WET Labs ECO [NTU]
    cast{castnums(i)}.D= ASCContinuousCTDData (:,11); %Depth, meters
    cast{castnums(i)}.Sal= ASCContinuousCTDData (:,13); %sal11: Salinity, Practical [PSU]
    cast{castnums(i)}.O2 = ASCContinuousCTDData (:,14)/(.022391); %Sbeox0ML/L converted to micromoles/L
    cast{castnums(i)}.SvCM = ASCContinuousCTDData (:,15); %also not sure about this one, not defined in .ros files
end
%%
figure(1); clf
    castToPlot = 8;
plot (cast{castToPlot}.T1, cast{castToPlot}.D) 
axis ij

%%
figure (2); clf
plot (cast{3}.O2, cast{3}.D); hold on;
plot (cast{4}.O2, cast{4}.D) 
axis ij

%%
figure (3); clf
plot (cast{10}.T1, cast{10}.D, 'y'); hold on;
plot (cast{11}.T1, cast{11}.D, 'm'); hold on;
plot (cast{12}.T1, cast{12}.D, 'g'); 
axis ij
%legend ('cast10', 'cast11', 'cast12')
ylabel('Depth (meters)')
xlabel('Temperature (Celcius)')

%%
figure (4); clf
plot (cast{5}.O2, cast{5}.D, 'g'); hold on;
axis ij
legend ('cast5')
ylabel('Depth (meters)')
xlabel('Raw Oxygen')

%%
figure (5); clf
plot (cast{7}.Sal, cast{7}.D); hold on; %FLMA
plot (cast{15}.Sal, cast{15}.D); hold on;%FLMB
plot (cast{5}.Sal, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('Salinity (PSU)')
legend ('Cast 7', 'Cast 15', 'Cast 5')
axis ij

%%
figure (6); clf
plot (cast{7}.T1, cast{7}.D); hold on; %FLMA
plot (cast{15}.T1, cast{15}.D); hold on;%FLMB
plot (cast{5}.T1, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('Temperature (Celcius)')
legend ('Cast 7', 'Cast 15', 'Cast 5')
axis ij

%%
figure (7); clf
plot (cast{7}.O2, cast{7}.D); hold on; %FLMA
plot (cast{15}.O2, cast{15}.D); hold on;%FLMB
plot (cast{5}.O2, cast{5}.D); hold on;%SUMO
axis ij
ylabel('Depth (meters)')
xlabel('Oxygen')
legend ('Cast 7', 'Cast 15', 'Cast 5')
%%
figure (7); clf
plot (cast{10}.F, cast{10}.D); hold on; %Glider1
plot (cast{21}.F, cast{21}.D); hold on;%Glider2
ylabel('Depth (meters)')
xlabel('Fluorescence')
axis ij

%%
figure(8); clf
plot (cast{5}.Sal, cast{5}.D); hold on;%SUMO %how to plot 2 xaxis on one graph
plot (cast{5}.Temp, cast{5}.D); hold on;%SUMO

axis ij