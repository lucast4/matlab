function [allffstim, alltraindir, allexptnum, allbirdname, alldaynums] ...
    = lt_opto_ExtrAllDat(ExptList, normeachexpt, valfield, groupbybird, ...
    onlylongepoch, StartDaySkipTime)

%% === lt 11/27/18 - extract all experiments

% normeachexpt=1; % NORMALIZE EACH EXPT (by mean of means over pos and
% % neg train
% valfield = 'All_FFmeanStimMinusNostim';
% % valfield = 'All_Dprime';
% % valfield = 'All_FFmedianStimMinusNostim';
% groupbybird = 0; % if 0, then by expt.
% onlylongepoch = 0;
% StartDaySkipTime = 0;
%

%%

% ############################ [EXTRACT SUMMARY DATA]
% =============== 1) EACH EXPERIEMNT COLLECT
OUTSTRUCT = struct;

for j=1:length(ExptList)
    
    if isempty(ExptList(j).dirtoplot)
        % then tyr to extract previously saved summary data
        
        tmp = load(ExptList(j).ctxtsummarydat);
        structbird = tmp.dat;
        
    else
        structbird = lt_opto_ExtrBirdDat(ExptList(j).dirtoplot, ...
            ExptList(j).twind, ExptList(j).SwitchTimes, ...
            0, StartDaySkipTime, onlylongepoch);
    end
    
    OUTSTRUCT.exptnum(j).dat=structbird;
end


% ================ 2) [ETRACTION] COMBINE INTO COMMON MATRICES
allffstim = [];
alltraindir = [];
allexptnum = [];
allbirdname = {};
alldaynums = [];
for j=1:length(ExptList)
    
    ffdiff_stim = OUTSTRUCT.exptnum(j).dat.(valfield);
    try
    daynum = OUTSTRUCT.exptnum(j).dat.All_Daynum;
    catch err
        daynum = OUTSTRUCT.exptnum(j).dat.All_daynums;
    end
    traindir = OUTSTRUCT.exptnum(j).dat.All_traindir;
    birdname = ExptList(j).birdname;
    

    % --------------------- NORMALIZE EACH EXPT (by mean of means over pos and
    % neg train
    if (0) %NOTE: moved to below grouping, so that the noarmlization is done after
        % grouping, so that nroamlziation does not affect significance of
        % the results.
        
        if normeachexpt==1
            y = grpstats(ffdiff_stim(traindir~=0), traindir(traindir~=0), {'mean'});
            assert(length(y)==2, 'then doesnt have both upa nd dn train ...');
            normval = mean(y);
            ffdiff_stim = ffdiff_stim-normval;
        end
    end
    
    allffstim = [allffstim; ffdiff_stim];
    alltraindir = [alltraindir; traindir];
    allexptnum = [allexptnum; ones(size(traindir,1),1)*j];
    tmp = cell(size(traindir));
    tmp(:) = {birdname};
    allbirdname = [allbirdname; tmp];
    alldaynums = [alldaynums; daynum];
end



% ===================== GROUP BY BIRD OR EXPT?
if groupbybird==1
    [allexptnum] = lt_tools_grp2idx({allbirdname});
end


% ====================== [NORMALIZATION?] if norm each experiment, then
% subtract mean effect of stim.
if normeachexpt==1
    for i=1:max(allexptnum)
        indsthis = allexptnum==i;
        
        y = grpstats(allffstim(alltraindir~=0 & indsthis), ...
            alltraindir(alltraindir~=0 & indsthis), {'mean'});
        assert(length(y)==2, 'then doesnt have both upa nd dn train ...');
        normval = mean(y);
        allffstim(indsthis) = allffstim(indsthis) - normval;
    end
end
