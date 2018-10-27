function output_disc = loaddates(num, txt)
    date = txt(2:end,3); %pull out text dates - this is in the format of a string
    output_disc = NaN*ones(length(txt)-1,1); 
    for i = 1:length(num)
        tf=isempty(date(i));
         if tf==1 
             output_disc(i) = NaN;
         else output_disc(i) = datenum(date(i));        
        end
    end
    