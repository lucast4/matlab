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
PARAMS.ffbinsedges = [20 35 80 130]; % edges, to plot timecourse in frequency bands

%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);

%% ======== FOR LEARNING, GET SWITCHES
SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct, 'multidaylearn');

%% ====== GET DATA FOR LEARNING (RELATIVE TO LEARNING SWITCHES)

% --- doesn't extract cohmat, just saves path to cohmat.
baseuseallinds =0 ;
pairtoget = 'LMAN-RA';
% pairtoget = 'LMANoutside-RAoutside';
% pairtoget = {'LMANoutside-RA', 'LMAN-RAoutside', 'LMANoutside-RAoutside'};
removeBadTrials = 0; % will only affect epoch inds.
SwitchCohStruct = lt_neural_Coher_LearnExtr2(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS, baseuseallinds, removeBadTrials);


%% ======== EXTRACT SCALARS

twind = [-0.07 -0.03]; % all combined
% twind = [-0.02 0]; % all combined
fwind = [22 32];

% ============= EXTRACT RELATIVE TO WN ONSET?
useWNtiming=0;
WNprctile = 2.5;
prewind_relWN = [-0.1 -0.05]; % rel the percentiel you want to use.

% === linearly interpolate coherence?
interpol =0; % across time  [NOT YET DONE].

[SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct,...
    SwitchStruct, PARAMS, twind, fwind, useWNtiming, WNprctile, prewind_relWN, ...
    interpol);


%% ====== EXTRACT LFP ACROSS LEARNING
close all;
collectAllProcess =1; % then colelcts not just Cohmat, but all phi and spectra.
% if 1, then will not collect each trial but just means (as takes too much
% meory). (mean pre and post and also diff)
plotON = 0;
averagechanpairs= 0; % for each motif, average over all chan pairs [NOTE: this is not up to date]
onlyfirstswitch = 0;
zscoreLFP = 3; % default 1, z-scores each t,ff bin separately.
% if 2, then doesn't zscore, instead normalizes as power proprotion
% (separately for each time slice and trial).
% if 3, then each ff window z-scored separately (good to see moudlation).
% if 4, then first 1) zscores within f, and 2) normalizes to all f (i.e.
% proportion
collectDiffMats = 0; % if 1, then collects differences (WN minus base). redundant, so leave at 0.
removeBadChans = 0; % default: 1
removeBadSyls = 1; % LEAVE AT 1.
typesToRemove = {'wn'}; % only remove syls that are bad becuase preceded by WN
% typesToRemove = {'wn', 'noise'};

if (0)
[OUTSTRUCT, OUTSTRUCT_CohMatOnly] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats, removeBadChans, typesToRemove);
end


[OUTSTRUCT] = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats, removeBadChans, typesToRemove);

