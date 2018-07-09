%% lt 4/4/18 - stores information about each experiment.
function [DATSTRUCT, Params] = lt_generalization_database(birdname, exptname)

%% ==================== FIND PARAMS FOR SPECIFIC BIRD/EXPERIMENT

if strcmp(birdname, 'wh6pk36') & strcmp(exptname, 'LMANlearn2')
    
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/wh6pk36/NEURAL/033017_LMANlearn2'};
    ListOfBatch = {...
        'BatchAll'};
MotifsToExtract = {'mk(s)d', 'mks(d)', '(v)bga', 'v(b)ga', 'vbg(a)',...
    'm(n)l', '(c)hb', 'c(h)b', 'ch(b)', 'kl(c)', 'nl(c)'};
    MotifsToRename = {};
    TargSyl = {'ch(b)'};
    TargLearnDirs = [1];
    % --- dates
    Date_WNon =  '30Mar2017-1123';
    Date_SwitchTimes = {'31Mar2017-1224'};

    
elseif strcmp(birdname, 'bu77wh13') & strcmp(exptname, 'LMANlearn1')
    
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/bu77wh13/NEURAL/020517_LMANlearn1'};
    ListOfBatch = {...
        'BatchAll'};
    MotifsToExtract = {'pj(a)b','pja(b)', ...
    'ij(b)', 'ijb(h)', ...
    'n(k)', 'o(k)', 's(k)'}; % note: removed all syls following target
    MotifsToRename = {...
        {'ijb(h)'}, '(b)_oldh'};
    TargSyl = {'pja(b)'};
    TargLearnDirs = [1];
    % --- dates
    Date_WNon =  '05Feb2017-1718';
    Date_SwitchTimes = {'07Feb2017-1509', '12Feb2017-1309'};


elseif strcmp(birdname, 'or74bk35') & strcmp(exptname, 'LMANneural2')
    
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/or74bk35/NEURAL/051217_LMANneural2'};
    ListOfBatch = {...
        'BatchAll'};
    MotifsToExtract = {'ja(b)', 'ba(b)', 'an(b)', 'h(b)', ...
                        'ja(n)', 'ba(n)', ...
                        'm(a)', 'a(a)', 'j(a)'}; % NOTE: could potentially add g (but is after WN, but long gap)
                    % COULD potnetially convert n to b (since is samish
                    % syl)
    MotifsToRename = {};
    TargSyl = {'an(b)'};
    TargLearnDirs = [1];
    % --- dates
    Date_WNon =  '12May2017-1616';
    Date_SwitchTimes = {};

elseif strcmp(birdname, 'pu69wh78') & strcmp(exptname, 'RALMANlearn1')
    
    % --- PARAMS
    ListOfDirs_ALL = {...
        '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1'};
    ListOfBatch = {...
        'BatchAll.LABELED'};
    MotifsToExtract = {'a(b)', 'j(b)', 'ab(h)', 'jb(h)',  'jbh(h)', '(g)'};
    MotifsToRename = {};
    TargSyl = {'a(b)'};
    TargLearnDirs = [1];
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
    TargLearnDirs = [1];
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
    TargLearnDirs = [1];
    
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
Params.TargLearnDirs = TargLearnDirs;
Params.MotifList = MotifsToExtract;
Params.IsSame = IsSame;
Params.IsTarg = IsTarg;
Params.Date_WNon =  Date_WNon;
Params.Date_SwitchTimes = Date_SwitchTimes;
Params.MotifsToRename = MotifsToRename;



