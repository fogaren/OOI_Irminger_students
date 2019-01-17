function output_disc = loaddates(num, txt)
    date = num(2:end,2); %pull out text dates - this is in the format of a string
    output_disc = NaN*ones(length(txt)-1,1); 
    for i = 1:length(date)
        tf=isempty(date(i));
         if tf==1 
             output_disc(i) = NaN;
         else
             output_disc(i) = datenum(date(i)) - 1;
        end
    end
    