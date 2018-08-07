%% lt 4/4/18 - stores information about each experiment.
function [DATSTRUCT, Params] = lt_generalization_database(birdname, exptname)

%% ==================== FIND PARAMS FOR SPECIFIC BIRD/EXPERIMENT

[ListOfDirs_ALL, ListOfBatch, MotifsToExtract, MotifsToRename, ...
    TargSyl, TargLearnDirs, Date_WNon, Date_SwitchTimes] = ...
    lt_generalization_exptlist(birdname, exptname);

%% ---- EXTRACT WN HITS (WILL SKIP IF DETECTS ALREADY DONE)
if (0)
lt_batchsong_extractWNhit(ListOfDirs_ALL, ListOfBatch)
end


%% --- EXTRACT DAT FOR THIS BIRD/EXPT
DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_ALL, {}, ListOfBatch, MotifsToExtract);


%% determine same type and different type

IsTarg = [];
IsSame = [];
for j=1:length(MotifsToExtract)
    motifthis = MotifsToExtract{j};
    
    if any(strcmp(motifthis, TargSyl))
        IsTarg = [IsTarg 1];
    else
        IsTarg = [IsTarg 0];
    end
    
    % --- is this motif same type to any of the targets?
    issamevec = [];
    for k=1:length(TargSyl)
        tsylthis = TargSyl{k};
        [issame] = lt_neural_QUICK_SylPairType(motifthis, tsylthis);
        issamevec = [issamevec issame];
    end
    issame = any(issamevec); % if is same type as any target, can call it same
    IsSame = [IsSame issame];
end


%% ================ combine all params
Params.TargSyl = TargSyl;
Params.TargLearnDirs = TargLearnDirs;
Params.MotifList = MotifsToExtract;
Params.IsSame = IsSame;
Params.IsTarg = IsTarg;
Params.Date_WNon =  Date_WNon;
Params.Date_SwitchTimes = Date_SwitchTimes;
Params.MotifsToRename = MotifsToRename;



