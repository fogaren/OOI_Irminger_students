function [time_grid, data_grid] = meanTimeInterval(timein, datain, time_step, begtime, endtime)

%Takes mean over data at regular time intervals - particularly valuable for
%noisy data collected at uneven time intervals
% INPUTS
% timein (matlab datenum format) - vector
% datain - data to take mean over, same length as timein
% time_step - interval for time_grid (in days)
% begtime and endtime - beginning and end times for time_grid
% OUTPUTS
% time_grid (matlab datenum format)
% data_grid - mean of datain at each time in time_grid (mean over interval centered around time point)

time_grid = [begtime: time_step: endtime];
data_grid = NaN*ones(length(time_grid),min(size(datain)));
for i = 1:length(time_grid)
    ind = find(timein > (time_grid(i) - time_step/2) & timein <= (time_grid(i) + time_step/2));
    data_grid(i,:) = nanmean(datain(ind,:));
end

end
