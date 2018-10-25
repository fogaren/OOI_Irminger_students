function loaddates(num, txt, output_disc)
    date = datenum(num(1:end,2)); %pull out text dates from num using datenum so it knows they're dates
    output_disc = NaN*ones(length(txt)-1,1); 
    for i = 1:length(num)
        tf=isempty((i));
         if tf==1 
             output_disc(i) = NaN;
         else output_disc(i) = date(i)-1        
        end
    end