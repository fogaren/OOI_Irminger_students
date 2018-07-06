%Script for importing excel data from ASC files, which contain continous
%CTD data from each cast
path = ['/Users/lucywanzer/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/'];
file = ['ASC_Continuous_CTD_Data.xlsx'];
filename= [path file]
%xlsread(filename)
ASCContinuousCTDData = xlsread(filename)
%salinity00=ASCContinuousCTDData(:,9)