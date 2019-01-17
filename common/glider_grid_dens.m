%% Grid data onto even density grids (based on glider_grid.m)
function [out] = glider_grid(daten,lat,lon,dens,profile_index,profile_direction,scivars,dens_grid);

% INPUT DATA
% daten --> time in Matlab format
% lat --> use latitude_interp
% lon --> use longitude_interp
% dens --> use dens_interp
% profile_index
% scivars --> array that can include multiple variables (T, S, density, O2, etc.)
% dens_grid --> potential densities on which to interpolate (e.g. [1026.5:0.015:1027.9])

%Set tolerance for number of good points in profile to actually keep any
%data
tol = 5;

%Number of profile indices
num_profiles = floor(max(profile_index));
num_scivars = min(size(scivars)); %assumes here that there will be more timestamps than variables

scivars_grid = NaN*ones(length(dens_grid),num_scivars,num_profiles);
time_start = NaN*ones(num_profiles,1);
duration = NaN*ones(num_profiles,1);
lat_profile = NaN*ones(num_profiles,1);
lon_profile = NaN*ones(num_profiles,1);
for i = 1:num_profiles
    ind = find(profile_index == i); %find all data points in the given profile
    if length(ind) > 0 & sum(~isnan(scivars(ind,1))) > tol
        scivars_grid(:,:,i) = naninterp1(dens(ind),scivars(ind,:),dens_grid);
        time_start(i) = min(daten(ind));
        duration(i) = max(daten(ind)) - min(daten(ind));
        lat_profile(i) = nanmean(lat(ind));
        lon_profile(i) = nanmean(lon(ind));   
        updown_profile(i) = nanmean(profile_direction(ind));
    end
end

%Find profile indices with good data
ind_nan = squeeze(sum(sum(isnan(scivars_grid)),2));
    tol = 0.95; %only keep data from profiles with at least 5% of data available
ind_keep = find(ind_nan < tol*length(dens_grid)*num_scivars);

%Output gridded data for all profiles with good data
out.profile_ind = ind_keep;
out.scivars = scivars_grid(:,:,ind_keep);
out.time_start = time_start(ind_keep);
out.duration = duration(ind_keep);
out.lat = lat_profile(ind_keep);
out.lon = lon_profile(ind_keep);
out.profile_direction = updown_profile(ind_keep);