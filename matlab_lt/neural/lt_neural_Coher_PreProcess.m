%% lt 1/3/19 - essential preprocessing after load previously saved data

%% ==== MAKE SURE SUMMARY STRCUT IS THE CORRECT ONE
SummaryStruct = MOTIFSTATS_Compiled.SummaryStruct;

%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);

%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;

%% ====== update params to match COHSTRUCT
settmp = 1;
if isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif)
    settmp=2;
    assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif));
end
PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).t_relons;
PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(settmp).motif(1).ffbins;
PARAMS.ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands

%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);

    %% ======== FOR LEARNING, GET SWITCHES
SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);

%% ====== GET DATA FOR LEARNING (RELATIVE TO LEARNING SWITCHES)

    % --- doesn't extract cohmat, just saves path to cohmat.
 baseuseallinds =0 ;
pairtoget = 'LMAN-RA';
% pairtoget = 'LMANoutside-RA';
pairtoget = {'LMANoutside-RA', 'LMAN-RAoutside', 'LMANoutside-RAoutside'};
SwitchCohStruct = lt_neural_Coher_LearnExtr2(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS, baseuseallinds);    


%% ======== EXTRACT SCALARS

twind = [-0.08 -0.03];
fwind = [25 35];
twind = [-0.08 -0.03];
fwind = [22 32];
% twind = [-0.05 -0.0];
% fwind = [25 40];
% twind = [-0.09 -0.02];
% fwind = [15 40];
twind = [-0.07 -0.03]; % all combined
fwind = [25 35];

[SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct,...
    PARAMS, twind, fwind);


%% ====== EXTRACT LFP ACROSS LEARNING
close all;
collectAllProcess =1; % then colelcts not just Cohmat, but all phi and spectra.
% if 1, then will not collect each trial but just means (as takes too much
% meory). (mean pre and post and also diff)
plotON = 0;
averagechanpairs= 0; % for each motif, average over all chan pairs [NOTE: this is not up to date]
onlyfirstswitch = 0;
removeBadSyls = 1; % LEAVE AT 1.
zscoreLFP = 3; % default 1, z-scores each t,ff bin separately.
% if 2, then doesn't zscore, instead normalizes as power proprotion
% (separately for each time slice and trial).
% if 3, then each ff window z-scored separately (good to see moudlation).
% if 4, then first 1) zscores within f, and 2) normalizes to all f (i.e.
% proportion
collectDiffMats = 0; % if 1, then collects differences (WN minus base). redundant, so leave at 0.

if (0)
[OUTSTRUCT, OUTSTRUCT_CohMatOnly] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);
end


[OUTSTRUCT] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);

%% for each case get wn minus base, mean coh scalar
cohscal_diff = [];
for i=1:length(OUTSTRUCT.bnum)
    
    indsbase = OUTSTRUCT.indsbase_epoch{i};
    indswn = OUTSTRUCT.indsWN_epoch{i};
    
    cohscal = OUTSTRUCT.cohscal{i};
    
    cohdiff = mean(cohscal(indswn)) - mean(cohscal(indsbase));
    
    cohscal_diff = [cohscal_diff; cohdiff];
end

OUTSTRUCT.cohscal_diff = cohscal_diff;

%% ====== COLLECT LEARN DIR AT TARGET

OUTSTRUCT = lt_neural_LFP_GetLearnDir(OUTSTRUCT, SwitchStruct);
    

