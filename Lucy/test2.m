depth = Yr1_wfp.depth_dosta;
daten = Yr1_wfp.time_dosta;
profile_index = Yr1_wfp.profile_index; %%Lucy added this
num_profiles = floor(max(profile_index));
num_scivars = min(size(scivars)); %assumes here that there will be more timestamps than variables

scivars_grid = NaN*ones(length(depth_grid),num_scivars,num_profiles);
time_start = NaN*ones(num_profiles,1);
duration = NaN*ones(num_profiles,1);
lat_profile = NaN*ones(num_profiles,1);
lon_profile = NaN*ones(num_profiles,1);
updown_profile = NaN*ones(num_profiles,1);
for i = 1:num_profiles
    ind = find(profile_index == i); %find all data points in the given profile
    if length(ind) > 0
        scivars_grid(:,:,i) = naninterp1(depth(ind),scivars(ind,:),depth_grid);
        time_start(i) = min(daten(ind));
        duration(i) = max(daten(ind)) - min(daten(ind));
        lat_profile(i) = nanmean(lat(ind));
        lon_profile(i) = nanmean(lon(ind));   
        updown_profile(i) = nanmean(profile_direction(ind));
    end
end
