gridded_in = Yr1_wfpgrid;
indup = find(gridded_in.updown > 0);
j = 0;
for i = 1:length(indup) %for each up profile
    %if subsequent profile after is down and time between profile starts is less than given tolerance
    if gridded_in.updown(indup(i) + 1) < 0
        if diff(gridded_in.time_start(indup(i): indup(i)+1)) < tol
            j = j + 1; %advance index grid
            scivars_pair(:,:,j) = mean(gridded_in.scivars(:,:,indup(i): indup(i)+1),3);
            ind_pair(j) = indup(i);
        end
    end
end