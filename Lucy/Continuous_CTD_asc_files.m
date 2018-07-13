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
    cast{castnums(i)}.SA = gsw_SA_from_SP(cast{castnums(i)}.Sal,cast{castnums(i)}.Pres,60,19);
    cast{castnums(i)}.CT = gsw_CT_from_t(cast{castnums(i)}.SA,cast{castnums(i)}.T1,cast{castnums(i)}.Pres);
    cast{castnums(i)}.rho0 = gsw_rho(cast{castnums(i)}.SA,cast{castnums(i)}.CT,0);
    cast{castnums(i)}.O2sol = gsw_O2sol(cast{castnums(i)}.SA,cast{castnums(i)}.CT,0,60,19);
    cast{castnums(i)}.O2sol2 = gsw_O2sol_SP_pt(cast{castnums(i)}.Sal,cast{castnums(i)}.T1);
    cast{castnums(i)}.aou = cast{castnums(i)}.O2sol - cast{castnums(i)}.O2 ;
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
subplot(1,2,1)
plot (cast{3}.rho, cast{3}.D, 'g'); hold on;
axis ij
legend ('cast22')
ylabel('Depth (meters)')
xlabel('Density')

subplot(1,2,2)
plot (cast{3}.aou, cast{3}.D, 'g'); hold on;
axis ij
legend ('cast22')
ylabel('Depth (meters)')
xlabel('AOU')

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
subplot(1,2,1)
for i = [1:8,10:12]
plot (cast{i}.rho, cast{i}.D); hold on; %FLMA
end
% plot (cast{15}.T1, cast{15}.D); hold on;%FLMB
% plot (cast{5}.T1, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('Density')
legend ('Cast 1', 'Cast 2', 'Cast 3', 'Cast 4', 'Cast 5', 'Cast 6', 'Cast 7', 'Cast 8', 'Cast 10', 'Cast 11')
axis ij

subplot(1,2,2)
for i = [1:8,10:12]
plot (cast{i}.O2, cast{i}.D); hold on; %FLMA
end
% plot (cast{15}.T1, cast{15}.D); hold on;%FLMB
% plot (cast{5}.T1, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('O2sol')
legend ('Cast 1', 'Cast 2', 'Cast 3', 'Cast 4', 'Cast 5', 'Cast 6', 'Cast 7', 'Cast 8', 'Cast 10', 'Cast 11', 'Cast 12')
axis ij


%%
figure (6); clf
subplot(1,2,1)
for i = [14:22]
plot (cast{i}.rho, cast{i}.D); hold on; %FLMA
end
% plot (cast{15}.T1, cast{15}.D); hold on;%FLMB
% plot (cast{5}.T1, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('Density')
legend ('Cast 14', 'Cast 15', 'Cast 16', 'Cast 17', 'Cast 18', 'Cast 19', 'Cast 20', 'Cast 21', 'Cast 22')
axis ij


subplot(1,2,2)
for i = [14:22]
plot (cast{i}.aou, cast{i}.D); hold on; %FLMA
end
% plot (cast{15}.T1, cast{15}.D); hold on;%FLMB
% plot (cast{5}.T1, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('aou')
legend ('Cast 14', 'Cast 15', 'Cast 16', 'Cast 17', 'Cast 18', 'Cast 19', 'Cast 20', 'Cast 21', 'Cast 22')
axis ij

%%
figure (11); clf
for i = [1:8,10:12,14:22]
plot (cast{i}.O2, cast{i}.D); hold on; %FLMA
end
% plot (cast{15}.T1, cast{15}.D); hold on;%FLMB
% plot (cast{5}.T1, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('Oxygen (micromoles/liter)')
legend ('Cast 1', 'Cast 2', 'Cast 3', 'Cast 4', 'Cast 5', 'Cast 6', 'Cast 7', 'Cast 8', 'Cast 10', 'Cast 11', 'Cast 12', 'Cast 14', 'Cast 15', 'Cast 16', 'Cast 17', 'Cast 18', 'Cast 19', 'Cast 20', 'Cast 21', 'Cast 22', 'Location','southwest')
axis ij

%%
figure (8); clf
plot (cast{7}.CT, cast{7}.D); hold on; %FLMA
plot (cast{15}.CT, cast{15}.D); hold on;%FLMB
plot (cast{5}.CT, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('Density')
legend ('Cast 7', 'Cast 15', 'Cast 5')
axis ij

%%
figure (7); clf
plot (cast{7}.O2, cast{7}.D); hold on; 
plot (cast{8}.O2, cast{8}.D); hold on;
plot (cast{11}.O2, cast{11}.D); hold on;
plot (cast{12}.O2, cast{12}.D); hold on;
ylabel('Depth (meters)')
xlabel('Oxygen (micromoles/L)')
legend ('Cast 7', 'Cast 8', 'Cast 11', 'Cast 12')
axis ij

%%
figure (10); clf
plot (cast{7}.aou, cast{7}.D); hold on; %FLMA
plot (cast{15}.aou, cast{15}.D); hold on;%FLMB
plot (cast{5}.aou, cast{5}.D); hold on;%SUMO
axis ij
ylabel('Depth (meters)')
xlabel('AOU')
legend ('Cast 7', 'Cast 15', 'Cast 5')

%%
figure (9); clf
plot (cast{10}.F, cast{10}.D); hold on; %Glider1
plot (cast{21}.F, cast{21}.D); hold on;%Glider2
ylabel('Depth (meters)')
xlabel('Fluorescence')
axis ij



%%
%figure(8); clf
%plot (cast{5}.Sal, cast{5}.D); hold on;%SUMO %how to plot 2 xaxis on one graph
%plot (cast{5}.Temp, cast{5}.D); hold on;%SUMO
%axis ij

%%
figure (6); clf
subplot(1,2,1)
plot (cast{7}.O2sol, cast{7}.D); hold on; %FLMA
plot (cast{15}.O2sol, cast{15}.D); hold on;%FLMB
plot (cast{5}.O2sol, cast{5}.D); hold on;%SUMO
ylabel('Depth')
xlabel('O2sol')
legend ('Cast 7', 'Cast 15', 'Cast 5')
axis ij

subplot(1,2,2)
plot (cast{7}.O2sol2, cast{7}.D); hold on; %FLMA
plot (cast{15}.O2sol2, cast{15}.D); hold on;%FLMB
plot (cast{5}.O2sol2, cast{5}.D); hold on;%SUMO
ylabel('Depth (meters)')
xlabel('O2sol2')
legend ('Cast 7', 'Cast 15', 'Cast 5')
axis ij
