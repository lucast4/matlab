%% ANALYSIS SHOULD DO IN THIS ORDER

%% ====== [OPTIONAL] LOAD PREVIOUSLY RECALCULATED COHERENCE
% E.G. This must be done if want to subtract trial-shifted coherence
% BEFORE THIS: must first run lt_neural_LFP_RecalcCoh

cd(['/bluejay0/bluejay2/lucas/analyses/neural/COHERENCE/RecalcCoh/' PARAMS.savemarker])

tmp = load('savestruct_cohere_PlusShuff.mat');

OUTSTRUCT.CohMean_Base = tmp.savestruct.dat.CohMean_Base;
OUTSTRUCT.CohMean_WN = tmp.savestruct.dat.CohMean_WN;
OUTSTRUCT_CohMatOnly = tmp.savestruct.dat.CohAlltrials;
OUTSTRUCT_CohMatOnly_shift = tmp.savestruct.dat.CohAlltrials_shuff;

PARAMS.tbins_old = PARAMS.tbins;
PARAMS.ffbins_old = PARAMS.ffbins;
PARAMS.tbins = tmp.savestruct.dat.t;
PARAMS.ffbins = tmp.savestruct.dat.f;

clear tmp;


%% ==================== [RECALCULATE COHERENCE MATRIX and SCLALAR
% BASED ON NEW SET OF TRIALS]
% This allows me to:
% 1) decide which trials to use
% 2) whether to normalize relative to shuffle
% 3) Whether to recalculate scalar relative to timing of WN 
% NOTE: This must be run BEFORE Realign coherence matrices...

% ===== params for recalc of sclar
% twind = [-0.07 -0.03]; % all combined
% fwind = [22 32];

twind = [-0.07 -0.03]; % all combined
fwind = [22 36];
% fwind = [25 35];

% - EXTRACT RELATIVE TO WN ONSET?
useWNtiming=1;
WNprctile = 2.5;
prewind_relWN = [-0.1 -0.05]; % rel the percentiel you want to use.

% wntouse = 'half';
wntouse = 'quarter';
% wntouse = 'third';
% wntouse = 'firsthalf';
% wntouse = 'half';

RemoveIfTooFewTrials =0;

% ========= remove bad syl
removebadsyl=1;

% --- norm to shift?
normtoshuff =1;
normtype = 'minus';
% normtype = 'zscore';
if normtoshuff==1
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(SwitchStruct, OUTSTRUCT, OUTSTRUCT_CohMatOnly, ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT, RemoveIfTooFewTrials, removebadsyl, normtoshuff, normtype, OUTSTRUCT_CohMatOnly_shift);
else    
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RecalcMat(SwitchStruct, OUTSTRUCT, OUTSTRUCT_CohMatOnly, ...
    SwitchCohStruct, PARAMS, twind, fwind, wntouse, useWNtiming, ...
    prewind_relWN, COHSTRUCT, RemoveIfTooFewTrials, removebadsyl, normtoshuff);
end


%% ========= [COHERENCE MATRIX] REALIGN BY WN ONSET
% Realign all cohernece matrices by WN time of target.
useallwn = 1; % default is 0 (just epoch) but for some experiments not enough data in those epoch?
[OUTSTRUCT, PARAMS] = lt_neural_Coher_RealignbyWN(OUTSTRUCT, SwitchCohStruct, ...
    SwitchStruct, PARAMS, useallwn);


%% ###################################### PLOTTING FUNCTIONS
%% ##########################################################

% SUMMARY PLOTS (compare diff syl types ...)
close all;
sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs
plotAllSwitchRaw = 0;
clim = [-0.1 0.1];

% =============
fieldtoplot = 'coher';
% 'coher'
% 'spec'
birdstoplot = [];
expttoplot = [];
% swtoplot = [1 7 9 11];
swtoplot = [];
useAbsVal = 0; 

% --- to get specific switch types. ... [is done in addition to above
% fitlers]
% swtoget = {}; % passes if matches ANY of these
swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0], [1 -1], [-1 1]}; % passes if matches ANY of these
% firstswitchfortarget_withinday = 1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday=1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday);

% NOTE: to get raw cohgram (before subtract) cd /blucurrently need to do
% breakpoint using spec and evaluate cohgram version instead. Should
% modify to plot cohgram.
timewindowtoplot = [-0.08 0]; % for spectra.
ffbinsedges = [20 35 80 130]; % edges, to plot timecourse in frequency bands
lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, timewindowtoplot, ...
    zscoreLFP, expttoplot, swtoplot, ffbinsedges, indtoget_b_e_s, useAbsVal);



%% ################# [MOTIFPLOT - MAIN] SUMMARIZE SCALAR CHANGE IN LEARNING
% I.E. EACH syl in order...
% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS

close all;

swtoget = {[0 1], [0 -1]}; % passes if matches ANY of these
% swtoget = {[1 0], [-1 0]}; % passes if matches ANY of these
firstswitchfortarget_withinday =1; % if 1, then onlky keeps if all targets 
% for a given switch did not have a previous switch on the same day
firstswitchofday =1;
indtoget_b_e_s = lt_neural_LEARN_FilterSwitches(SwitchStruct, swtoget, ...
    firstswitchfortarget_withinday, firstswitchofday);
% indtoget_b_e_s = [];

lt_neural_Coher_SumPlotMotifs(OUTSTRUCT, SwitchStruct, MOTIFSTATS_Compiled, PARAMS, ...
    indtoget_b_e_s);

%% ################## [GOOD] [SCALAR & PITCH LERANING- COMPARE TO LEARNING RATE]
close all;
nsegs = 4; % quartiles, then 4... (for both coh and learning)
learnConvertToRate = 0; % then learning (z) is z dividided by mean time of bin (hours).
lt_neural_Coher_PitchLearnCoh(OUTSTRUCT, PARAMS, SwitchCohStruct, ...
    SwitchStruct, nsegs, learnConvertToRate);


