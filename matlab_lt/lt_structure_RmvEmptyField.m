function Structout = lt_structure_RmvEmptyField(Structin)
%% lt 11/6/18 - removes all fields that are empty

flist = fieldnames(Structin);
for i=1:length(flist)
    fthis = flist{i};
    
    if isempty(Structin.(fthis))
        Structin = rmfield(Structin, fthis);
    end
end

Structout = Structin;