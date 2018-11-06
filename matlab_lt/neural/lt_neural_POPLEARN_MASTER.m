%% ==== 1) EXTRACT MOTIFSTATS POP
%% AD HOC CHANGES TO THE SETS OF NEURONS (E.G. TO MAXIMIZE DATASET SIZE)
if (0)
% ===== 1) pu69, combining sets so that can look at learning.
i=1; 
ii=1;
assert(strcmp(SummaryStruct.birds(i).birdname, 'pu69wh78'));
assert(strcmp(SummaryStruct.birds(i).exptnum_pop(ii).exptname, 'RALMANOvernightLearn1'));
assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{4} == [18 20]));
assert(all(SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons{5} == [18 19 20]));
% --- get new
newfiles = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles{4:5}];
newneurons = [18 20];
% -- remove old
SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons(4:5) = [];
SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles(4:5) = [];
% --- add
SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_neurons ...
    newneurons];
SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles = [SummaryStruct.birds(i).exptnum_pop(ii).Sets_songfiles ...
    {newfiles}];
end

%% ======================== EXTRACT SEGMENTS FOR POPULATIONS

close all; clear MOTIFSTATS_Compiled;
collectWNhit=0;
onlyCollectTargSyl=0;
LearnKeepOnlyBase = 0;
saveOn = 1;
OrganizeByExpt =0;
collectFF=1;
Params_regexp.motif_predur = [];
Params_regexp.motif_postdur = [];
Params_regexp.preAndPostDurRelSameTimept = 1;
Params_regexp.RemoveIfTooLongGapDur = [];
Params_regexp.extractDirSong = 1;

MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl, OrganizeByExpt,...
    collectFF, [], Params_regexp);


%% ==== MAKE SURE SUMMARY STRCUT IS THE CORRECT ONE
SummaryStruct = MOTIFSTATS_Compiled.SummaryStruct;


%% ==== ALIGN ALL MOTIFS TO COMMON MOTIFS FOR PLOTTING ACROSS EXPERIMENTS

% ==== 1) quickly list all motifs
disp('========================');
numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:numbirds
    bname = MOTIFSTATS_Compiled.birds(i).birdname;
    MotiflistAll = {};
    numneurons = length(MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons);
    for ii=1:numneurons
        
        motiflist = [MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(ii).motif_regexpr_str];
        MotiflistAll = [MotiflistAll motiflist];
        
    end
    MotiflistAll = unique(MotiflistAll);
    disp(['BIRD: ' birdname]);
    disp(MotiflistAll);
end


% ========================== ANNOTATE IN MOTIFSTATS COMPILED

%% ==== GET PARAMS
clear PARAMS;

PARAMS.motif_predur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_predur;
PARAMS.motif_postdur = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.params.motif_postdur;
PARAMS.alignbyonset = MOTIFSTATS_Compiled.birds(1).MOTIFSTATS.neurons(1).motif(1).Params.REGEXP.alignByOnset;
assert(~isempty(PARAMS.alignbyonset));
% PARAMS.savemarker = '14Oct2018_2147';
PARAMS.savemarker = input('what is save marker? (e.g. 14Oct2018_2147)? ', 's');

%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);

%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;



%% #######################################################################
%% ############################# OLD STUFF
%% ================ PLOT [CORRELATION WITH FF]
close all;
xcov_dattotake = [-0.01 0.05];
xcov_dattotake = [-0.08 0.04];
xcov_dattotake = [-0.075 0.025];
xcovwindmax = 0.04;
binsize_spk = 0.0025;

MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk);


%% #######################################################################
%% ############################# LEARNING STUFF;
%% ==== 2) EXTRACT LEARNING SWITCH STRUCT

SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);


%% =========== SUMMARIZE LEARNING TRAJECTORY (PLUS NEURON SETS)
close all;

% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn1'};
% motiftoplot = 'c(b)';
% BirdExptPairsToPlot = {'pu69wh78', 'RALMANOvernightLearn1'};
% motiftoplot = 'aa(b)';
BirdExptPairsToPlot = {};
motiftoplot = '';

lt_neural_POPLEARN_PlotLearnTraj(MOTIFSTATS_pop ,SwitchStruct, ...
    SummaryStruct, BirdExptPairsToPlot, motiftoplot);


%% #######################################################################
%% ############################# CROSS CORR CHANGE DURING LAERNING
%% ================ PLOT CROSS CORR WRT TO LEARNING
close all;
BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn1'};
% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn2'};
SwitchToPlot = [2];
BregionWantedList = {{'LMAN', 'RA'}};
onlyPlotIfBothPrePostTrials = 0;
lt_neural_POPLEARN_Plot(MOTIFSTATS_pop, SwitchStruct, BirdExptPairsToPlot, ...
    SwitchToPlot, BregionWantedList, onlyPlotIfBothPrePostTrials);


%% ================ SUMMARIZE CROSS CORRELATION OVER COURSE OF EXPERIMENT
% over multiple switches
close all;
exptnum = [1];
birdnum = [2];
BregionWantedList = {{'LMAN', 'RA'}};
[OUTSTRUCT, birdnum] = lt_neural_POPLEARN_Summary(MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, exptnum, BregionWantedList);


%% ======= PLOT SUMMARY OF XCOV traces
close all;
lt_neural_POPLEARN_SummaryPlot1(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, exptnum)

%% ======= PLOT MEAN XCOV
close all;
windowmean = [-0.05 0.02]; % window, in ms, relative to lag = 0;
lt_neural_POPLEARN_SummaryPlot2(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, windowmean)



%% ======== SUMMARIZE OVER ALL EXPERIMENTS
windowmean = [-0.05 0.02]; % window, in s, relative to lag = 0;
SkipIfTargsDiffSyls = 1; % skips switches where targets are different syl types.

lt_neural_POPLEARN_SummaryPlot3(OUTSTRUCT, MOTIFSTATS_pop, SwitchStruct, ...
    birdnum, windowmean, SkipIfTargsDiffSyls)




%% [GOOD] ################ FOR EACH BIRD, SUMMARIZE ACROSS ALL EXPTS
% DONE: removed dir inds
% TO DO: 1) baseline, 2) use early or late period in epoch

% INDICATE BY HAND WHICH NEURON SETS ARE APPROPRIATE
lt_neural_POPLEARN_SumTraj_Input; % go in here and select which dataset

% --- RUN
close all;
bregionwanted = {'LMAN', 'RA'};
lt_neural_POPLEARN_SumTraj(MOTIFSTATS_pop, SwitchStruct, ...
    metadatstruct, bregionwanted);

%% ================ PLOT PAIRED RASTERS WRT TO LEARNING
% [IN PROGRESS!!!]
lt_neural_POPLEARN_PairRast




%% ###################################################################
%% ############################################ COHERENCE
% OVERALL: for each motif, look at coherence of raw data, aligned to syl
% onset. Does that change during learning?

% NOTE: to see some old progress, see:
lt_neural_MasterScript_Pop;


%% ################# OLD VERSION, EXTRACTING COHERENCE FIRST.

lt_neural_Coher_Script1;


if (0)
% ===== convert from segextract to a specific song and time in song
lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, i, neurthis);
end


%% ################## NEW VERSION, FIRST EXTRACTING LFP, THEN DOING STUFF WTIH LFP
%% ==== EXTRACT LFP AND SPEC STRUCT, WITH ORGANIZATION MATCHED TO COHSTRUCT
close all;
skipifOnlyOneChan = 1; % i.e. if a given dataset/motif not paired, then skip. 
BirdsToPlot = {'pu69wh78', 'wh44wh39'};
% SetsToSkip = {'1-2-2'};
SetsToSkip = {};

LFPSTRUCT = lt_neural_LFP_ExtractStruct(MOTIFSTATS_pop, SummaryStruct, ...
    MOTIFSTATS_Compiled, skipifOnlyOneChan);

% ================ save
if (0)
    marker = '14Oct2018_2147';
    fname = ['/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_' marker '.mat'];
    save(fname, 'LFPSTRUCT');
    
    save(['/bluejay5/lucas/analyses/neural/LFP/PARAMS_' PARAMS.savemarker '.mat'], 'PARAMS');
end

% ==== load 
if (0)
    load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_14Oct2018_2147.mat');
    load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_14Oct2018_2147.mat');
    load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_14Oct2018_2147.mat');
    load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/14Oct2018_2147/COHSTRUCT.mat');
end

%% ====== CALCULATE COHERENCE USING LFP
% GOAL TO GET COHSTRUCT
close all;
savemarker = '14Oct2018_2147';
COHSTRUCT = lt_neural_LFP_GetCohStruct(LFPSTRUCT, PARAMS, SummaryStruct);

%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);

%% ====== GET DATA FOR LEARNING (RELATIVE TO LEARNING SWITCHES)

if (0)
pairtoget = 'LMAN-RA';
SwitchCohStruct = lt_neural_Coher_LearnExtr(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS);
else
    % --- doesn't extract cohmat, just saves path to cohmat.
pairtoget = 'LMAN-RA';
SwitchCohStruct = lt_neural_Coher_LearnExtr2(COHSTRUCT, MOTIFSTATS_pop, ...
    SwitchStruct, pairtoget, LFPSTRUCT, PARAMS);    
end

%% ====== EXTRACT LFP ACROSS LEARNING
close all;
collectAllProcess =1; % then colelcts not just Cohmat, but all phi and spectra.
% if 1, then will not collect each trial but just means (as takes too much
% meory). (mean pre and post and also diff)
plotON = 0;
averagechanpairs= 0; % for each motif, average over all chan pairs [NOTE: this is not up to date]
onlyfirstswitch = 0;
removeBadSyls = 1; % LEAVE AT 1.
zscoreLFP = 1; % default 1, z-scores each t,ff bin separately.
collectDiffMats = 0; % if 1, then collects differences (WN minus base). redundant, so leave at 0.
OUTSTRUCT = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);
clear SwitchCohStruct;

    
    
%% ====== SUMMARY PLOT OF COHERENCE LARNING

% SUMMARY PLOTS (compare diff syl types ...)
close all;
sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs
plotAllSwitchRaw = 0;
clim = [-0.05 0.05];
fieldtoplot = 'Spec1Mean_WNminusBase';
% Spec1Mean_WNminusBase
% Spec2Mean_WNminusBase
% CohMean_WNminusBase
% === v1 - % NOTE: this is obsolete. see below. TO RUN THIS with coher,
% make sure in extraction of OUTSTRUCT included collectDiffMats.
lt_neural_Coher_Learn_PlotSum(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot);

% ============== THIS IS BETTER - subsumes the above, more compacta ndf
% flexible.
fieldtoplot = 'coher';
% 'coher'
% 'lfp'
lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot);



%% ===== PLOT PHASE DISTRIBUTIONS

figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


assert(all(strcmp(OUTSTRUCT.bregionpair, 'LMAN-RA')), 'assumes chan1 is LAMN, chang 2 is RA...');
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

% ---- params
ffbintoplot = 4;
tbintoplot = 10;
nrand = 10; % random trials to plot
nanglebins = 16;
ntoplot = 10; % random channels to plot

% === targets
indslist = find(OUTSTRUCT.istarg==1)';
indslist = indslist(randperm(length(indslist), ntoplot));
for j=1:length(indslist)
    indthis = indslist(j);
    
    trialsbase = OUTSTRUCT.indsbase_epoch{indthis};
    trialsWN = OUTSTRUCT.indsWN_epoch{indthis};
    phimat = OUTSTRUCT.PhiMat{indthis};    
    maxtrials = min([length(trialsbase) length(trialsWN)]);
    
    % ================== BASELINE
    phithis = squeeze(phimat(:, ffbintoplot, trialsbase));
    basewn_label = 'BASELINE';
    
    % -- 1) plot a few trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('random trials');
    ylabel(basewn_label);
    indrand = randperm(size(phithis,2), nrand);
    plot(tbins, phithis(:, indrand), '-x');
    axis tight;
    
    % -- 2) plot all trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    phitmp = phithis(:, randperm(size(phithis,2), maxtrials));
    title(['all trials, N=' num2str(size(phitmp,2))]);
    ylabel('mean, std');
    plot(tbins, phitmp, 'kx');
    plot(tbins, phitmp+2*pi, 'kx');
%     phimean = mean(phithis,2);
%     phistd = std(phithis, [], 2);
%     lt_plot(tbins, phimean, {'Errors', phistd, 'Color', 'r'});
    axis tight
    ylim([-pi 2*pi])
    
    % -- 3) plot circular distribution of angles
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel(['at tbin:' num2str(tbintoplot)]);
    y = phithis(tbintoplot, :);
%     y(y<0)=y(y<0)+2*pi;
    rose(y);
    
    % -- 4) estimate phase locking value
    

end



%% ===== PLOT LFP FOR MULTIPLE TRIALS OF ALL MOTIFS

close all;



%% ===== COMPUTE SPECTROGRAMS USING EXTRACTED

