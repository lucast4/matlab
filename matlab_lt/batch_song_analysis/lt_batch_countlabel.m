function NumCasesAll = lt_batch_countlabel(batch, strtocheck)
%% lt 3/18/18 - counts how often a string is found in batch

fid = fopen(batch);

songf = fgetl(fid);

NumCasesAll = [];
while ischar(songf)
    
    try
    tmp = load([songf '.not.mat']);
    catch err
    end
    numcases = length(strfind(tmp.labels, strtocheck));
    
    NumCasesAll = [NumCasesAll numcases];
    
   songf = fgetl(fid); 
end

