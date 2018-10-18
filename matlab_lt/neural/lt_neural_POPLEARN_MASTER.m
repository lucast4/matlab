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



%% ================ EXTRACT COHSTRUCT

COHSTRUCT = lt_neural_Coher_ExtrCohStruct(MOTIFSTATS_pop, SummaryStruct, ...
    MOTIFSTATS_Compiled);


% ============== SAVE COHERENCE
if (0)
    marker = '14Oct2018_2147';
    fname = ['/bluejay5/lucas/analyses/neural/COHERENCE/COHSTRUCT_' marker '.mat'];
    save(fname, 'COHSTRUCT');
end



%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);


%% ================= FIX PROBLEM THAT COHMATS HAVE DIFFERENT SIZES IN COHSTRUCT.
tlengthdesired = 20; % number of time bins desired - ad hoc to make sure all same dimension
fflengthdesired = 19;
COHSTRUCT = lt_neural_Coher_FixCohMatSize(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct, ...
    tlengthdesired, fflengthdesired);

%% ############### SUMMARY PLOT OF COHERANCE
% PLOT COHGRAM AND FF BANDS FOR EVERY CHANNEL PAIR AND MOTIF

close all;
lt_neural_Coher_PlotRawSum(COHSTRUCT, MOTIFSTATS_pop);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUMMARY ANALYSES
%% ================ 1) EXTRACTION (EXTRACT MEAN COHGRAMS FOR ALL)
tlengthdesired = 20; % number of time bins desired - ad hoc to make sure all same dimension
fflengthdesired = 19;

%  Series of summary analyses - here just does extraction
[All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, ...
    All_chanpair, All_bregionpair, All_bregionpair_alphaorder, ...
    All_tbins, All_ffbins] = lt_neural_Coher_Summ_Extr(COHSTRUCT, ...
    MOTIFSTATS_pop, SummaryStruct, tlengthdesired, fflengthdesired);

%% ======= PLOT SUMMARY
close all;

pairthis = 'LMAN-RA';
ffbinsedges = [15 30 80 150]; % edges, to plot timecourse in frequency bands
lt_neural_Coher_Summ_Plot1(pairthis, All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, ...
    All_chanpair, All_bregionpair, All_bregionpair_alphaorder, ffbinsedges)

%% =================== FOR EACH MOTIF AND BRAIN REGION, ONE PLOT
close all;

MotifLists = {...
    {'pu69wh78', {'(a)ab', 'a(a)b', 'aa(b)', 'aab(h)', 'aabh(h)'}}, ...
    {'pu69wh78', {'(j)jb', 'j(j)b', 'jj(b)', 'jjb(h)', 'jjbh(h)', 'jjbhh(g)'}}, ...
    {'pu69wh78', {'(j)jbhh', 'j(j)bhh', 'jj(b)hh', 'jjb(h)h'}}, ...
    {'pu69wh78', {'h(g)'}}, ...
    {'wh44wh39', {'(j)n', '(n)hh', 'n(h)h', 'nh(h)'}}, ...
    {'wh44wh39', {'(m)d', '(d)kcc', 'd(k)cc', 'dk(c)c', 'dkc(c)', 'c(b)', 'cb(b)'}}, ...
    {'wh44wh39', {'(n)h', 'n(h)'}}, ...
    {'wh44wh39', {'(d)kc', 'd(k)c', 'dk(c)'}}, ...
    };
BregionPair = 'LMAN-RA';

tbins = All_tbins{1};
ffbins = All_ffbins{1};
PARAMS.tbins = tbins;
PARAMS.ffbins = ffbins;
lt_neural_Coher_Summ_Plot2(MotifLists, BregionPair, MOTIFSTATS_Compiled, ...
    All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, All_chanpair, ...
    All_bregionpair, All_bregionpair_alphaorder, tbins, ffbins,ffbinsedges);

%% ##################################################### 
%% ################################## LEARNING [COLLECT, FIR EACH SIWTCH]

pairtoget = 'LMAN-RA';
SwitchCohStruct = lt_neural_Coher_LearnExtr(COHSTRUCT, MOTIFSTATS_pop, SwitchStruct, pairtoget);

%% ################### LEARNING [PLOT EACH EXPERIMENT]
% =================== PLOT, ONE FIGURE FOR EACH CHANNEL PAIR
% ======= ie plots all cases.
close all;
lt_neural_Coher_Learn_PlotAll;

%% ======= 1) EXTRACT
close all;
plotON=0; % ONE PLOT FOR EACH MOTIF/CHANPAIR [NOTE: use PlotChanMotif below instead]
averagechanpairs= 0; % for each motif, average over all chan pairs
onlyfirstswitch = 0;
removeBadSyls = 1; % LEAVE AT 1.
OUTSTRUCT = lt_neural_Coher_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls);


%% ======= 2) GET SCALAR MEASURE OF CHANGE IN COHERENCE IN CERTAIN WINDOW

premotor_wind = [-0.07 0.02];
premotor_wind = [-0.07 0];
ff_wind = [12 40];

% ###################### MEAN OF POST MINUS PRE
tinds = tbins>premotor_wind(1) & tbins<premotor_wind(2);
finds = ffbins>ff_wind(1) & ffbins<ff_wind(2);
CohMean_WNminusBase_scalar = nan(length(OUTSTRUCT.CohMean_WNminusBase),1);
for i=1:length(OUTSTRUCT.CohMean_WNminusBase)
    tmp = OUTSTRUCT.CohMean_WNminusBase{i}(tinds, finds); % get t,f, window desired
    CohMean_WNminusBase_scalar(i) = nanmean(tmp(:));
end
OUTSTRUCT.CohMean_WNminusBase_scalar = CohMean_WNminusBase_scalar;
PARAMS.scalar_premotor_wind = premotor_wind;
PARAMS.scalar_ff_wind = ff_wind;

% ####################### ALL TRIALS
OUTSTRUCT.CohMat_scalar = cell(length(OUTSTRUCT.CohMean_WNminusBase),1);
for i=1:length(OUTSTRUCT.CohMean_WNminusBase)
    
    cohmat = OUTSTRUCT.CohMat{i};
    
    tmp = squeeze(nanmean(nanmean(cohmat(tinds, finds, :),1),2));
    
    OUTSTRUCT.CohMat_scalar{i} = tmp;
end

%% ========== RAW PLOTS
close all;
% ---- what to plot?
birdtoplot = 'wh44wh39';
expttoplot = 'RALMANlearn4';
swnum = 1;


% ##################################### 1) ======= PLOT RAW (EACH MOTIF AND CHAN PAIR)
% ----
averagechanpairs= 0; % for each motif, average over all chan pairs
removeBadSyls = 1; % LEAVE AT 1.
lt_neural_Coher_PlotChanMotif(SwitchStruct, SwitchCohStruct, OUTSTRUCT, ...
    birdtoplot, expttoplot, swnum, averagechanpairs, PARAMS, removeBadSyls)


% ##################################### 1) ======= PLOT RAW (EACH MOTIF, AVERAGE CHAN PAIR)
% ----
averagechanpairs= 1; % for each motif, average over all chan pairs
removeBadSyls = 1; % LEAVE AT 1.
lt_neural_Coher_PlotChanMotif(SwitchStruct, SwitchCohStruct, OUTSTRUCT, ...
    birdtoplot, expttoplot, swnum, averagechanpairs, PARAMS, removeBadSyls)



% 3) ################### ONE PLOT FOR EACH CHANPAIR (AVERAGES OVER MOTIFS)
plotrawcohdiff = 1; % if 0, then takes average over motifs. if 1, then plots each raw
swtoplot = swnum;

% assert(averagechanpairs==0, 'following needs each pair broken out');
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;
clim = [-0.15 0.15];

% ======= for each channel pair, get mean for targ, nontarg
% -- gets motifs for each channel pair
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
indsgrp_unique = unique(indsgrp);

figcount=1;
subplotrows=3;
subplotcols=10;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for j=indsgrp_unique'
   
    % ----- info...
    bname = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).birdname;
    ename = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).exptnum(unique(OUTSTRUCT.enum(indsgrp==j))).exptname;
    swnum = unique(OUTSTRUCT.switch(indsgrp==j));
    chpair = OUTSTRUCT.chanpair(indsgrp==j, :);
    chpair = chpair(1,:);
    
    if ~(strcmp(birdtoplot, bname) & strcmp(ename, expttoplot) & swnum==swtoplot)
        continue
    end
    
    % ---------- target
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    if ~any(indsthis)
        continue
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    
    if isempty(cohmat)
        continue
    end
    
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title('TARG');
    ylabel({[bname '-' ename '-sw' num2str(ss)], ['ch:' num2str(chpair)]});
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotrawcohdiff);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    % ---------- same-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    if ~any(indsthis)
        continue
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('SAME');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotrawcohdiff);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    % -------- diff type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    if ~any(indsthis)
        continue
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('DIFF');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotrawcohdiff);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
    
    % ------------ TARGET MINUS SAMETYPE
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    cohmat1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    cohmat2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    cohmat = nanmean(cohmat1,3)-nanmean(cohmat2,3);
    
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG-SAME');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
 
    
     % ------------ TARGET MINUS DIFFTYPE
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    cohmat1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    cohmat2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    cohmat = nanmean(cohmat1,3)-nanmean(cohmat2,3);
    
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG-DIFF');
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
end



% ##################################### 2) PLOT TIMECOURSE OF COHERENCE SCALAR CHANGES
swnum = swtoplot;

% =================== RUN
bnumthis = find(strcmp(birdtoplot, {MOTIFSTATS_pop.birds.birdname}));
enumthis = find(strcmp(expttoplot, {MOTIFSTATS_pop.birds(bnumthis).exptnum.exptname}));

indsthis = OUTSTRUCT.bnum==bnumthis & OUTSTRUCT.enum==enumthis ...
    & OUTSTRUCT.switch==swnum;

% PLOT EACH MOTIF SEPARATELY
motifnum_unique = unique(OUTSTRUCT.motifnum(indsthis));
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for mm=motifnum_unique'
    indsthis = OUTSTRUCT.bnum==bnumthis & OUTSTRUCT.enum==enumthis ...
        & OUTSTRUCT.switch==swnum & OUTSTRUCT.motifnum==mm;
    
    motifname = SwitchCohStruct.bird(bnumthis).exptnum(enumthis).switchlist(swnum).motifnum(mm).motifname;
    istarg = unique(OUTSTRUCT.istarg(indsthis));
    issame = unique(OUTSTRUCT.issame(indsthis));
    if istarg==1
        pcoltit = 'r';
    elseif issame ==1
        pcoltit = 'b';
    else
        pcoltit = 'k';
    end
    
    t = OUTSTRUCT.tvals(indsthis); assert(length(unique(cellfun(@mean, t)))==1, 'should all be same /./');
    t = t{1};
    ff = SwitchCohStruct.bird(bnumthis).exptnum(enumthis).switchlist(swnum).motifnum(mm).ffvals;
    cohscalall = OUTSTRUCT.CohMat_scalar(indsthis);
    cohscalall = cell2mat(cellfun(@transpose, cohscalall, 'UniformOutput', 0));
    
    % ------ only get inds during base and WN
    indsbase = OUTSTRUCT.indsbase(indsthis); assert(length(unique(cellfun(@sum, indsbase)))==1, 'should all be same /./');
    indsWN = OUTSTRUCT.indsWN(indsthis); assert(length(unique(cellfun(@sum, indsWN)))==1, 'should all be same /./');
    indsbase = indsbase{1};
    indsWN = indsWN{1};
    
    t = t(indsbase | indsWN);
    ff = ff(indsbase | indsWN);
    cohscalall = cohscalall(:, indsbase | indsWN);
    
    % ================ PLOT T
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdtoplot '-' expttoplot '-sw' num2str(swnum)]);
    ylabel([motifname], 'Color', pcoltit);
    xlabel('t');
    plot(t, ff, 'xk');
    axis tight
    % -- line for WN on
    tswitch = mean(t(sum(indsbase):sum(indsbase)+1));
    line([tswitch tswitch], ylim, 'Color', 'r');
    
    % ================ PLOT COHERENCE, EACH CHANNEL
    for j=1:size(cohscalall,1)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['chanpair: ']);
        ylabel('coherence (in t/f wind)');
        xlabel('t');
        cohscal = cohscalall(j,:);
        plot(t, cohscal, 'xk');
        axis tight;
        ylim([0.1 0.9]);
    % -- line for WN on
    tswitch = mean(t(sum(indsbase):sum(indsbase)+1));
    line([tswitch tswitch], ylim, 'Color', 'r');
    end
    
    
    
end

%% ======= SUMMARIZE SCALAR RESULT, OVER ALL SWITCHES, MOTIFS, AND CHANNELS
close all;
clim = [-0.2 0.2];

indsgrp_switch = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
indsgrp_chanpair = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, ...
    OUTSTRUCT.chanpair});

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indsgrp_switch_unique = unique(indsgrp_switch);
indsgrp_chanpair_unique = unique(indsgrp_chanpair);
for i=1:length(indsgrp_switch_unique)
   swgrpthis = indsgrp_switch_unique(i);
   bnum = unique(OUTSTRUCT.bnum(indsgrp_switch==swgrpthis));
   enum = unique(OUTSTRUCT.enum(indsgrp_switch==swgrpthis));
   swnum = unique(OUTSTRUCT.switch(indsgrp_switch==swgrpthis));
   bname = MOTIFSTATS_Compiled.birds(bnum).birdname;
   ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
   
   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
   title([bname '-' ename '-sw' num2str(swnum)]);
    ylabel(['WN - base, coh']);
    
   % ====== plot each channel pair its own line
   for chanpair = indsgrp_chanpair_unique'
       indsthis = indsgrp_switch==swgrpthis & indsgrp_chanpair==chanpair;
       if ~any(indsthis)
           continue
       end
       motifID = lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)); % ---- get positions within global motif
       cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
       assert(length(unique(motifID)) == length(motifID), 'each motif shoudl have at most 1 datapoint');
       % ---- sort in order of motifs
       [~, indsort] = sort(motifID);
       motifID = motifID(indsort);
       cohscal = cohscal(indsort);
       plot(motifID, cohscal, 'o-k');
   end
   lt_plot_zeroline;
   
   % ====== overall
   indsthis = indsgrp_switch==swgrpthis;
   
   istarg = OUTSTRUCT.istarg(indsthis);
   motifs = OUTSTRUCT.motifname(indsthis);
   motifID = lt_neural_QUICK_MotifID(bname, motifs); % ---- get positions within global motif
   cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
   
   [ymean, ysem] = grpstats(cohscal, motifID, {'mean', 'sem'});
   x = unique(motifID);
   lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
   
   % -------- NOTE DOWN POSITION OF TARGET SYSL
   indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==1;
   xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
   plot(xtarg, clim(1)+0.02, '^r');
   
   % ------- NOTE POSITION OF SAME_TYPES
   indsthis = indsgrp_switch==swgrpthis & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
   xtarg = unique(lt_neural_QUICK_MotifID(bname, OUTSTRUCT.motifname(indsthis)));
   if ~isempty(xtarg)
   plot(xtarg, clim(1)+0.02, '^b');
   end
   
   % ----- labels
   [~, motiflabels] = lt_neural_QUICK_MotifID(bname);
    set(gca, 'XTick', 1:length(motiflabels));
    set(gca, 'XTickLabel', motiflabels);
    rotateXLabels(gca, 90);
    ylim(clim);
    
end



%% ========= TRIAL BY TRIAL CORRELATION BETWEEN COHERENCE AND PITCH?

% ===== 1) for all data, calculate cross correlation



%% SUMMARY PLOTS (compare diff syl types ...)
close all;
sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs
plotAllSwitchRaw = 0;
clim = [-0.1 0.1];
lt_neural_Coher_Learn_PlotSum(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim);



%%  ######################################### GRAND MEAN FIGURE:
% incorporates functions above with appropriate params.
% see notes within script.
close all;
lt_neural_Coher_Learn_PlotGrandMean;



%% ##################
%% ################### LEARNING, FOR EACH BIRD, SUMMARIZE ACROSS ALL EXPTS
% DONE: removed dir inds
% TO DO: 1) baseline, 2) use early or late period in epoch

% INDICATE BY HAND WHICH NEURON SETS ARE APPROPRIATE
lt_neural_POPLEARN_SumTraj_Input; % go in here and select which dataset

% --- RUN
close all;
bregionwanted = {'LMAN', 'RA'};
lt_neural_POPLEARN_SumTraj(MOTIFSTATS_pop, SwitchStruct, ...
    metadatstruct, bregionwanted);



%% plot single cases
if (0)
    lt_figure; hold on;
    
    % === 1) plot multiple trials
    
    
    % === 2) plot average over all trials
    lt_subplot(3,2,2); hold on;
    
    chanpair = 1;
    dim1 = size(CohAllTrials{1},1);
    dim2 = size(CohAllTrials{1},2);
    cohmat = CohAllTrials(chanpair, :);
    
    % -- convert to 3d mat
    cohmat = reshape(cell2mat(cohmat), dim1, dim2, []);
    cohmean = mean(cohmat, 3);
    imagesc(t, ff, 10*log10(cohmean'));
    
    % ==== 3) plot specific ff
    lt_subplot(3,2,3); hold on;
    plot(cohmean(:, 3), 'ok')
    
    
    % ############################### PLOT EACH FREQUENCY BIN SEPARATELY
    figcount=1;
    subplotrows=5;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    for i=1:length(ff)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        fthis = ff(i);
        title(['ffbin:' num2str(fthis)]);
        hsplots = [hsplots hsplot];
        
        cohthis = squeeze(cohmat(:, i, :));
        ymean = mean(cohthis,2);
        ysem = lt_sem(cohthis');
        x = 1:length(ymean);
        lt_plot(x, ymean, {'Errors', ysem, 'Color', 'k'});
    end
    linkaxes(hsplots, 'xy');
    
    % -- PLOT EACH FF AS A LINE WITH ERROR BARS
end

%% ############### SUMMARY PLOT - FIRST EXTRACT MEAN COHEROGRAM FOR EACH DATAPT

%% ===== convert from segextract to a specific song and time in song

lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, i, neurthis);


