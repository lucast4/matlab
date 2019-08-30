function lt_neural_SummaryStruct_ShowDirChans(SummaryStruct, longversion)
%% lt 3/21/19 - for all birds, plots directories and channels that are extracted

for i=1:length(SummaryStruct.birds)
   
    numneur = length(SummaryStruct.birds(i).neurons);
    alldir = {};
    allsubdir = {};
    allchans = [];
    allnuerID = [];
    for ii=1:numneur
       [a, b] = fileparts(SummaryStruct.birds(i).neurons(ii).dirname);
       alldir = [alldir; a];
       allsubdir = [allsubdir; b];
       allchans = [allchans; SummaryStruct.birds(i).neurons(ii).channel];
       allnuerID = [allnuerID; ii];
    end
    
    [alldir,~, ind2] = unique(alldir);
    if longversion == 0
    for j=1:length(alldir)
        disp([alldir{j} '  ---  ' num2str(unique(allchans(ind2==j))')]);
    end
    
    elseif longversion==1
        for j=1:length(alldir)
        disp([alldir{j}]);
        cellfun(@(x)disp(x), allsubdir(ind2==j));
%         num2str(unique(allchans(ind2==j))')]);
        disp(allnuerID(ind2==j));
        end

        cellfun(@(x)disp(x), allsubdir);
        disp('---');
    end
end
