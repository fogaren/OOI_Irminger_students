%Run the variable naming code that will load and name the files 
run('OOI_FileLoading_VariableNaming.m')

% The Files that are loaded and the structures created include :
    % Apex Surface Mooring = ApexS
    % Near Surface Mooring = NearS
    % Flanking A = FlankingA_strt
    % Flanking B = FlankingB_strt
    
% The variables that are read out of these files are as follows : 
    %Time 
    %O2_dosta  = not finalized O2 readings
    %O2 = finalized O2 readings 
    %temp = in C°, for ApexS only Sea Surfac
    %sal = practical salinity
    %p = pressure
    %Lat = Latitude
    %Lon = Longitude
    %O2_expected = calculated O2 values using the gsw O2 calcuated from 
    
%% Run the file that will create preliminary figures 
    % Commmentary for these preliminary figures is provided in the
    % Supplementary material of the files 
    
run('OO1_Figures_Y3.m')


%% Run the file seperates and plots nighttime only data

