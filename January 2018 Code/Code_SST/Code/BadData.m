%Loop to remove bad data
function val_cleaned =BadData(x,low,up)
%%%% Function to remove bad data based on values

%%%% INPUTS %%%%%%
% low = lower bound for good data
% up = upper bound for good data
% X = variable of interest

%%%% OUTPUTS %%%%%%
% val_cleaned = cleaned values within bounds


%% Loop over all points in time_in to separate day and night
val_cleaned = x;
  % return bounded value clipped between bl and bu
  idLow = find(x<low);
   val_cleaned(idLow) = NaN;
  idHigh = find (x > up);
    val_cleaned(idHigh) = NaN;
end

  