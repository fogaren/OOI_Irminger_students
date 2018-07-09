%Script for importing excel data from ASC files, which contain continous
%CTD data from each cast
path = ['/Users/lucywanzer/Dropbox/OOI_Irminger_students/Irminger5_CruiseDocs/'];
file = ['ASC Continuous CTD Data.xlsx'];
filename= [path file]
xlsread=filename
%CTD_data=xlsread ('ASC Continous CTD Data', 'Cast 001.asc');