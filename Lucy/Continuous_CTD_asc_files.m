%Script for importing excel data from ASC files, which contain continous
%CTD data from each cast
path = ['/Users/lucywanzer/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/'];
file = ['ASC_Continuous_CTD_Data.xlsx'];
filename= [path file]
ASCContinuousCTDData = xlsread(filename,1);
Pres= ASCContinuousCTDData (:,1)
T1= ASCContinuousCTDData (:,2); %T090C Temperature (ITS-90, deg C)
T2= ASCContinuousCTDData (:,3); %T190C Temperature (ITS-90, deg C)
C1= ASCContinuousCTDData (:,4); %COS/m Conductivity (S/m)
C2= ASCContinuousCTDData (:,5); %C1S/m Conductivity (S/m)
OR= ASCContinuousCTDData (:,6); %sbeox0V: Oxygen raw, SBE 43 [V]
F = ASCContinuousCTDData (:,7); %flECO-AFL: Fluorescence, WET Labs ECO-AFL/FL [mg/m^3]
Turb = ASCContinuousCTDData (:,8); %turbWETntu0: Turbidity, WET Labs ECO [NTU]
D= ASCContinuousCTDData (:,11); %Depth, meters
Sal= ASCContinuousCTDData (:,13); %sal11: Salinity, Practical [PSU]
O = ASCContinuousCTDData (:,14); %Sbeox0ML/L not sure what this variable is, not defined in .ros files
SvCM = ASCContinuousCTDData (:,15); %also not sure about this one, not defined in .ros files
plot (T2, D) 
axis ij
%salinity00=ASCContinuousCTDData(:,9)
%load CTDData001
%plot (CTDData001.Sal1, CTDData001.DepSM)
