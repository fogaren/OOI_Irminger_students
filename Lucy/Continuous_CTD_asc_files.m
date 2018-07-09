%Script for importing excel data from ASC files, which contain continous
%CTD data from each cast
path = ['C:/Users/Hilary/Dropbox/Wellesley/OOI_Irminger_students/Irminger5_CruiseDocs/'];
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
    cast{castnums(i)}.O2 = ASCContinuousCTDData (:,14); %Sbeox0ML/L not sure what this variable is, not defined in .ros files
    cast{castnums(i)}.SvCM = ASCContinuousCTDData (:,15); %also not sure about this one, not defined in .ros files
end
%%
figure(1); clf
    castToPlot = 8;
plot (cast{castToPlot}.T2, cast{castToPlot}.D) 
axis ij


