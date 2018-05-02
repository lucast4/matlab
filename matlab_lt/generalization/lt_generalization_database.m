%% lt 4/4/18 - stores information about each experiment.
function [DATSTRUCT, Params] = lt_generalization_database(birdname, exptname)

%% ==================== FIND PARAMS FOR SPECIFIC BIRD/EXPERIMENT
if strcmp(birdname, 'pu69wh78') & strcmp(exptname, 'RALMANlearn1')
    
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1'};
    ListOfBatch = {...
        'BatchAll.LABELED'};
    MotifsToExtract = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)',  'jbh(h)', '(g)'};
    MotifsToRename = {};
    TargSyl = {'a(b)'};
    % --- dates
    Date_WNon =  '01Nov2017-1402';
    Date_SwitchTimes = {'01Nov2017-2045'};
    
elseif strcmp(birdname, 'wh44wh39') & strcmp(exptname, 'RALMANlearn2')
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2'};
    ListOfBatch = {...
        'Batch1139to2028'};
    MotifsToExtract = {'c(b)', 'cb(b)', '(d)', '(k)', 'k(c)', 'kc(c)', '(n)', ...
        'n(h)', 'nh(h)'};
    MotifsToRename = {...
        {'c(b)', 'cb(b)'}, '(b)_targ', ...
        { 'k(c)', 'kc(c)'}, '(c)', ...
        {'n(h)', 'nh(h)'}, '(b)_oldh', ...
        }; % this will first rescale as deviation from mean (in hz)
    % and then combine. useful e.g. if diff syl across contexts.
    
%     TargSyl = {'c(b)', 'cb(b)'};
    TargSyl = {'(b)_targ'};
    % --- dates
    Date_WNon =  '14Mar2018-1338';
    Date_SwitchTimes = {};
    
elseif strcmp(birdname, 'wh44wh39') & strcmp(exptname, 'RALMANlearn3')
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/wh44wh39/NEURAL/032118_RALMANlearn3', ...
        '/bluejay5/lucas/birds/wh44wh39/NEURAL/032418_RALMANlearn3', ...
        '/bluejay5/lucas/birds/wh44wh39/NEURAL/040118_RALMANlearn3_audioonly'};
    ListOfBatch = {...
        'BatchAll', ...
        'BatchAllUNDIR',...
        'BatchAllUNDIR'};
    
    MotifsToExtract = {'c(b)', 'cb(b)', '(k)', 'k(c)', 'kc(c)', '(n)', ...
        'n(h)', 'nh(h)'};
    MotifsToRename = {...
        {'c(b)', 'cb(b)'}, '(h)_oldb', ...
        { 'k(c)', 'kc(c)'}, '(c)', ...
        {'n(h)', 'nh(h)'}, '(h)_targ', ...
        }; % this will first rescale as deviation from mean (in hz)
    % and then combine. useful e.g. if diff syl across contexts.
    
%     TargSyl = {'n(h)', 'nh(h)'};
    TargSyl = {'(h)_targ'};
    
    % --- dates
    Date_WNon =  '21Mar2018-1215';
%     Date_SwitchTimes = {'21Mar2018-1314', '21Mar2018-1325', '21Mar2018-2347'};
    Date_SwitchTimes = {'21Mar2018-2347'}; % removed short duration when WN on and off.
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
Params.MotifList = MotifsToExtract;
Params.IsSame = IsSame;
Params.IsTarg = IsTarg;
Params.Date_WNon =  Date_WNon;
Params.Date_SwitchTimes = Date_SwitchTimes;
Params.MotifsToRename = MotifsToRename;



