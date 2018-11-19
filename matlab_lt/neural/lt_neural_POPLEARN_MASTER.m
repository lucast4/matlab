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
    disp(['BIRD: ' bname]);
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




%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ COHERENCE
% OVERALL: for each motif, look at coherence of raw data, aligned to syl
% onset. Does that change during learning?

% NOTE: to see some old progress, see:
lt_neural_MasterScript_Pop;


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% OLD VERSION, EXTRACTING COHERENCE FIRST.

lt_neural_Coher_Script1;


if (0)
% ===== convert from segextract to a specific song and time in song
lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, i, neurthis);
end


%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ################## NEW VERSION, FIRST EXTRACTING LFP, THEN DOING STUFF WTIH LFP
% NOTE: [processing steps summary]
% 1) Extract LFP for each raw song file (save adjacent)
% 2) Extract those segmented data into LFPSTRUCT
% 3) Calculate (and save) coherence using data in LFPSTRUCT

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


%% +++++++++++++++++++++++++++++++++ PLOT RAW NEURAL, LFP, LFP(filtred)
if (0)
    % individual trials.
    % 2 VERSIONS:
    % 1) choose which motif to plot.
    % 2) go thru all switches and plot all target sylalbles (if multiple syls, does each separately
    
    %% ############################ VERSION 1- plot chosed motif
    
    % ======== 1) EXTRACT RAW DATA (NOT JUST LFP)
    close all;
    
    birdplot = 'pu69wh78';
    exptplot = 'RALMANOvernightLearn1';
    swplot = 1;
    motifplot = []; % [string] leave blank for target
    extrapad = 0.05; % seconds, pre and post...
    [DatAll, t_onoff, fs, bregionlist, chanlist_toget, i, ii, mm] = ...
        lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, ...
        SwitchCohStruct, birdplot, exptplot, swplot, motifplot, ...
        extrapad);
    
    
    % ============ 2) FOR BASE AND WN, PLOT N TRIALS OF LFP, NEURAL
    close all;
    ntoplot = 1; % trials pre and post
    filt_low = 25;
    filt_hi = 35;
    % filt_fs = fs;
    
    savedir = ['/bluejay5/lucas/analyses/neural/LFP/FIGS_PlotEgRaw/' PARAMS.savemarker];
    saveON = 1;
    
    lt_neural_LFP_PlotEgRaw(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
        i, ii, swplot, mm, ntoplot, filt_low, filt_hi, SwitchCohStruct, ...
        SwitchStruct, savedir, saveON);
    
    %% ############################ VERSION 2 - GO thru all swithches
    close all;
    for i=1:length(SwitchCohStruct.bird)
        for ii=1:length(SwitchCohStruct.bird(i).exptnum)
            
            for ss=1:length(SwitchCohStruct.bird(i).exptnum(ii).switchlist)
                
                if isempty(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum)
                    continue
                end
                
                % ========== which motif to plot?
                % get target
                motifplot_list = ...
                    SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end);
                
                % ================ GO THRU ALL MOTIFS
                for mm=1:length(motifplot_list)
                    motifplot = motifplot_list{mm};
                    
                    
                    % ================ RUN
                    birdplot = SwitchStruct.bird(i).birdname;
                    exptplot = SwitchStruct.bird(i).exptnum(ii).exptname;
                    extrapad = 0.05; % seconds, pre and post...
                    [DatAll, t_onoff, fs, bregionlist, chanlist_toget, ~, ~, motifnum] = ...
                        lt_neural_LFP_PlotEgRaw_Extract(PARAMS, SwitchStruct, MOTIFSTATS_pop, SummaryStruct, ...
                        SwitchCohStruct, birdplot, exptplot, ss, motifplot, ...
                        extrapad);
                    
                    % ================== PLTO AND SAVE EXAMPLES
                    close all;
                    ntoplot = 15; % trials pre and post
                    filt_low = 25;
                    filt_hi = 35;
                    % filt_fs = fs;
                    
                    savedir = ['/bluejay5/lucas/analyses/neural/LFP/FIGS_PlotEgRaw/' PARAMS.savemarker];
                    saveON = 1;
                    
                    lt_neural_LFP_PlotEgRaw(DatAll, t_onoff, fs, bregionlist, chanlist_toget, ...
                        i, ii, ss, motifnum, ntoplot, filt_low, filt_hi, SwitchCohStruct, ...
                        SwitchStruct, savedir, saveON);
                    
                end
                
            end
        end
    end
    
    
    %% ================ KEEP LIST OF CHANNELS TO EXCLUDE
    
    {...
        {'pu69wh78', 'RALMANlearn1', 1, 14}, ... % not bad. noise some trials.
        {'pu69wh78', 'RALMANOvernightLearn1', 1, 21}, ... % noise some trials.
        {'wh44wh39', 'RALMANlearn1', 1, 15}, ... % huge slow noise sometimes (rare), leads to power increase in beta
        }; % {bird, expt, swnum, chan}
    
    % ========= GOOD ONES:
    % pu69, RALMANlearn1, 2: rare trials at base see noise in all
    % channels... seems fine.
    
    % pu69, RALMANlearn2, 1: OK
    
    % pu69, RALMANlearn2, 2: OK
    
    % pu69, RALMANlearn2, 3: OK. rare trials see high hz noise -50ms, all
    % chans.

    % pu69, RALMANlearn2, 4: OK

    % pu69, RALMANlearn2, 5: OK

    % wh44wh39-RALMANlearn1-sw2 - OK
    
    % wh44wh39-RALMANlearn1-sw7-motif8 - OK
    
    % wh44wh39-RALMANlearn1-sw7-motif9 - OK

    % wh44wh39-RALMANlearn1-sw9-motif8 - OK

    % wh44wh39-RALMANlearn1-sw9-motif9 - OK
    
    % wh44wh39-RALMANlearn2-sw1-motif8 - OK
    
    % wh44wh39-RALMANlearn3-sw1-motif2 - OK
    
    % wh44wh39-RALMANlearn3-sw1-motif3 - OK
    
    % wh44wh39-RALMANlearn4-sw1-motif7 - OK
    
    % wh44wh39-RALMANlearn4-sw2-motif6 - OK
    
    % SUMMARY: Only noise that seems to affect beta band (25-35hz) is in
    % wh44 rl1, sw1, which is the large deviations. Also sometimes in wh44
    % ral3[? I can't remember] saw some large deviations that were
    % suspiociously synchronized between LMAN and RA, but not obvious
    % qwithin eacjh region that is noise...
    
end


%% #######################################################################
%% ====== CALCULATE COHERENCE USING LFP
% GOAL TO GET COHSTRUCT
% NOTE: will save cohstruct also.
close all;
savemarker = '14Oct2018_2147';
COHSTRUCT = lt_neural_LFP_GetCohStruct(LFPSTRUCT, PARAMS, SummaryStruct);

%% ################ NOTE DOWN BRAIN REGION PAIRS

COHSTRUCT = lt_neural_Coher_GetBrRegPairs(COHSTRUCT, MOTIFSTATS_pop, SummaryStruct);

%% ############ SAVE OVERALL PARAMS (inclyding time bins)
if exist('All_tbins', 'var')
    tbins = All_tbins{1};
    ffbins = All_ffbins{1};
    PARAMS.tbins = tbins;
    PARAMS.ffbins = ffbins;
    PARAMS.ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands
else
    assert(~isempty(COHSTRUCT.bird(1).experiment(1).setnum(2).motif));
    PARAMS.tbins = COHSTRUCT.bird(1).experiment(1).setnum(2).motif(1).t_relons;
    PARAMS.ffbins = COHSTRUCT.bird(1).experiment(1).setnum(2).motif(1).ffbins;
    PARAMS.ffbinsedges = [10 25 32 80]; % edges, to plot timecourse in frequency bands
end

%% ################# ANALYSES START FROM HERE: LOAD DATA
% ==== load 

load('/bluejay5/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/LFPSTRUCT_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PARAMS_14Oct2018_2147.mat');
load('/bluejay5/lucas/analyses/neural/LFP/PROCESSED/14Oct2018_2147/COHSTRUCT.mat');

    %% ======== FOR LEARNING, GET SWITCHES
SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);

    
%% ======= PLOT RAW LFP, PHI, SPECTROGRAM, COHEROGRAM

close all;
% birdtoplot = 'pu69wh78';
% expttoplot = 'RALMANOvernightLearn1';
% swnum = 1;
birdtoplot = 'wh44wh39';
expttoplot = 'RALMANlearn2';
swnum = 1;
motiftoplot = 'jjb(h)'; % ------ WILL plot target syl, unless say otherwise.
chanpairtoplot = [];
bregiontoget = 'LMAN-RA'; % WILL take random channel (for bregion pair) if no chanpairtoplot specified.
fs = 1500;
Nplot = 5; % n trials base and dur WN.
plotonlyLFP = 1; % then doesn't plot cohgram, spectra etc.
lt_neural_Coher_LFP_plottrials(SwitchCohStruct, [], SwitchStruct, ...
    PARAMS, birdtoplot, expttoplot, swnum, motiftoplot, chanpairtoplot, fs, ...
    Nplot, bregiontoget, plotonlyLFP);



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

PARAMS.bregionpair_toget = pairtoget;
%% ########################## COH CORRELATE WITH PITCH?
%% ======== EXTRACT SCALARS

twind = [-0.08 -0.03];
fwind = [22 32];

[SwitchCohStruct, PARAMS] = lt_neural_LFP_PitchCorr(COHSTRUCT, SwitchCohStruct,...
    PARAMS, twind, fwind);

%% ======= PLOT ALL DATA [COH SCALAR vs FF]
close all;
switchestoplot = [1];
birdstoplot = {};

numbirds = length(SwitchCohStruct.bird);
for i=1:numbirds
    bname = SummaryStruct.birds(i).birdname;
    numexpt = length(SwitchCohStruct.bird(i).exptnum);
    for ii=1:numexpt
        exptname = SummaryStruct.birds(i).exptnum_pop(ii).exptname;
        numswitch = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist);
        for ss=1:numswitch
           
            if ~isempty(switchestoplot)
               if ~any(ismember(switchestoplot, ss))
                   continue
               end
            end
            
            figcount=1;
            subplotrows=6;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];

            nummotif = length(SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum);
           for mm=1:nummotif
               disp([num2str(i) '-' num2str(ii) '-' num2str(ss)]);
               datthis = SwitchCohStruct.bird(i).exptnum(ii).switchlist(ss).motifnum(mm);
               if isempty(datthis.chanpair)
                   continue
               end
               
               motifname = datthis.motifname;
               tmp = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs(1:2:end);
               if any(strcmp(tmp, motifname))
                   istarg=1;
                   learndir = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningDirs{find(strcmp(tmp, motifname))+1};
               else
                   istarg=0;                   
               learndir = [];
               end
               
               % ==================== extract coherene scalars for this dat
               % ------------ COLLECT SCALAR
               t = datthis.tvals;
               ff = datthis.ffvals;
               cohscal = datthis.cohscalar;
               assert(size(t,2) == size(cohscal,1));
               
               % ------------ WN switch times
               swtime = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum;
               
               % ================ ONE PLOT FOR EACH CHANNEL PAIR
                nchanpairs = size(cohscal,2);
                hsplots = [];
                for cc = 1:nchanpairs
                    chpairthis = datthis.chanpair(cc,:);
                   
                    X = {};
                    Y = {};

                   % -- baseline
                   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                   hsplots = [hsplots hsplot];
                   title({[bname '-' exptname '-s' num2str(ss)], ...
                       [motifname '-ch' num2str(chpairthis) '[BASE]']});
                   indtrial = t<swtime;
                   pcol = 'k';
                   
                   xlabel('coh scalar');
                   ylabel('ff');
                   x = cohscal(indtrial,cc);
                   y = ff(indtrial);
                   plot(x, y, 'x', 'Color', pcol);
                   lt_regress(y, x, 0, 0, 1, 1, 'm', 1);
                   
                   X{1} = x;
                   Y{1} = y;
                   
                   % -- WN
                   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                   hsplots = [hsplots hsplot];
                   title({[bname '-' exptname '-s' num2str(ss)], ...
                       [motifname '-ch' num2str(chpairthis) '[WN]']});
                   indtrial = t>swtime;
                   pcol = 'm';
                   
                   xlabel('coh scalar');
                   ylabel('ff');
                   x = cohscal(indtrial,cc);
                   y = ff(indtrial);
                   plot(x, y, 'x', 'Color', pcol);
                   lt_regress(y, x, 0, 0, 1, 1, 'm', 1);
                   
                   X{2} = x;
                   Y{2} = y;

                   % ---- if target, note down learn dir
                   if istarg==1
                       lt_plot_text(min(x), max(y), ['targ, dir' num2str(learndir)], 'b');
                   end
                   
                   % -------------- OVERLAY DISTRIBUTIONS FOR BASE AND WN
                   [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                   title({[bname '-' exptname '-s' num2str(ss)], ...
                       [motifname '-ch' num2str(chpairthis)]});
                   xlabel('coh (mean, sem)');
                   ylabel('ff (mean, sem]');
                   xmean = cellfun(@mean, X);
                   ymean = cellfun(@mean, Y);
                   xerr = cellfun(@lt_sem, X);
                   yerr = cellfun(@lt_sem, Y);
                   lt_plot(xmean(1), ymean(1), {'Errors', yerr(1), 'Xerrors', xerr(1), 'Color', 'k'});
                   lt_plot(xmean(2), ymean(2), {'Errors', yerr(2), 'Xerrors', xerr(2), 'Color', 'r'});
                   xlim([0.2 0.8]);
                end
                linkaxes(hsplots, 'xy');
                              
           end
           pause;
           close all;
        end
    end
end


%% =========== EXTRACT COH SCALAR DATA



%% ################################################# LEARNING CHANGES
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
OUTSTRUCT = lt_neural_LFP_Learn_Extr(SwitchStruct, SwitchCohStruct, ...
    plotON, averagechanpairs, PARAMS, onlyfirstswitch, removeBadSyls, ...
    collectAllProcess, zscoreLFP, collectDiffMats);
clear SwitchCohStruct;

    
%% ====== COLLECT LEARN DIR AT TARGET

OUTSTRUCT = lt_neural_LFP_GetLearnDir(OUTSTRUCT, SwitchStruct);
    
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
fieldtoplot = 'spec';
% 'coher'
% 'spec'
birdstoplot = [1 2];
% NOTE: to get raw cohgram (before subtract) currently need to do
% breakpoint using spec and evaluate cohgram version instead. Should
% modify to plot cohgram.
timewindowtoplot = [-0.08 0]; % for spectra.
lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, sumplottype, ...
    plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, timewindowtoplot, ...
    zscoreLFP);

%% ###################################################################
%% ################ COHERENCE SCALAR PLOTS
% - Can predict change in coh based on ffvs coh correlation?
% - Is there signifincat correlation between coh and FF?
% - Shuffle stats for change in coherence across experiments.

%% ====== EXTRACT BASELINE CORRELATIONS.
useonlybaseepoch = 0; % if 1, then just limited epoch. if 0, then entie baseline data
% default: 0;
corrtype = 'spearman'; % spearman. or pearson or Kendall
cohdiff_usedprime = 0; % if 0, then just mean diff. if 1, then dprime
nboot = 10; % to get bootstrap SE

OUTSTRUCT = lt_neural_Coher_CohScalExtract(OUTSTRUCT, useonlybaseepoch, ...
    corrtype, cohdiff_usedprime, nboot);


%% ===== [PHI] CALCULATE THINGS

% tbintoget = [11 12];
% fbintoget = [3 4];
tbintoget = [11]; % NOTE: currently only support single index. 
fbintoget = [4];

% ============= ORDER MATTERS FOR PHI;
assert(all(strcmp(OUTSTRUCT.bregionpairs_originalorder, 'LMAN-RA')), 'assumes that is LMAN-->RA');


% ========== 1) PLV AND MEAN PHI
OUTSTRUCT.Phi_mean_BaseWN = nan(length(OUTSTRUCT.PhiMat), 2);
OUTSTRUCT.Phi_PLV_BaseWN = nan(length(OUTSTRUCT.PhiMat), 2);
for i=1:length(OUTSTRUCT.PhiMat)
   disp(i);
    phimat = OUTSTRUCT.PhiMat{i};
    inds_base = OUTSTRUCT.indsbase_epoch{i};
    inds_WN = OUTSTRUCT.indsWN_epoch{i};
    
    % ========== 1) PLV
    % - BASE
    indsthis = inds_base;
    
    phivec= squeeze(phimat(tbintoget, fbintoget, indsthis));
    [plv, mean_cart, mean_pol, bootstats] = lt_neural_QUICK_PhaseLockVal(phivec, 1);
    
    OUTSTRUCT.Phi_mean_BaseWN(i,1) = mean_pol;
    OUTSTRUCT.Phi_PLV_BaseWN(i,1) = plv;

    % - BASE
    indsthis = inds_WN;
    
    phivec= squeeze(phimat(tbintoget, fbintoget, indsthis));
    [plv, mean_cart, mean_pol, bootstats] = lt_neural_QUICK_PhaseLockVal(phivec, 1);
    
    OUTSTRUCT.Phi_mean_BaseWN(i,2) = mean_pol;
    OUTSTRUCT.Phi_PLV_BaseWN(i,2) = plv;
    
end


%% ==== [QUICK] plots PLV across target syls. [dirty, since not controlled for sampel size]
indsthis = OUTSTRUCT.istarg==1;
figure; hold on; 
lt_plot_histogram(OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1), [], 1, 1, [], 1, 'k');
lt_plot_histogram(OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2), [], 1, 1, [], 1, 'r');
figure; hold on;
plot([1 2], [OUTSTRUCT.Phi_PLV_BaseWN(indsthis,1) OUTSTRUCT.Phi_PLV_BaseWN(indsthis,2)], '-ok');
xlim([0 3]);


%% === [QUICK] Plots random data, single trial phi matrix.

lt_figure; hold on;
colormap('parula')
for j=1:8
    
i = randi(200);
ii = randi(750);
try
phimat = OUTSTRUCT.PhiMat{ii}(:,:,i);
catch err
i = randi(200);
ii = randi(750);
end

lt_subplot(4,2,j); hold on;

imagesc(PARAMS.tbins, PARAMS.ffbins, phimat', [-pi pi]);
colorbar
axis tight;
end

%% ===== [SUMMARY PLOT] LEARNING CHANGE IN COHSCALAR AS FUNCTION OF BASE CORR? 
% i.e. can I predict the change in coherence in a given expt (switch) based
% on relationship between baseline coh/ff and direction of training
% NOTE: can test variation within and between experiments (centerdata)
% NOTE: can test both normative learn dir and actual ff change (useFF ...);
close all;
cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel)
rho_norm_to_global = 0;
centerdata = 0; % then for each switch, centers data (i.e. across channels)
% NOTE: this is useful if want to ask about, within a tgiven expt, is there
% correlation. 
onlyfirstswitch = 0; % if 0, then all siwtches.
useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir

lt_neural_Coher_CohScalPlot(OUTSTRUCT, SwitchStruct, cohdiff_norm_to_global, ...
    rho_norm_to_global, centerdata, onlyfirstswitch, useFFchangeInsteadOfLearnDir);


%% ====== [SUMMARY PLOT] PHI DURING LEARNING
% 1) Changes in Phi consistency across trials
% 2) Changes in mean Phi (Wn vs. base)
% NOTE: Assume that order of data is LMAN-->RA (i.e. negative phi means LMAN leads). Code
% will make assertions so that this must be true.

close all;
cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel)
rho_norm_to_global = 0;
centerdata = 0; % then for each switch, centers data (i.e. across channels)
% NOTE: this is useful if want to ask about, within a tgiven expt, is there
% correlation. 
onlyfirstswitch = 1; % if 0, then all siwtches.
useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir

lt_neural_Coher_PhiSumPlot(OUTSTRUCT, SwitchStruct, cohdiff_norm_to_global, ...
    rho_norm_to_global, centerdata, onlyfirstswitch, useFFchangeInsteadOfLearnDir);



%% ====== [SUMMARY PLOT] SHUFFLE test of coherence change at target

% === 1) plot distribution of coh change and overlay tareget [each expt]



%% ===== PLOT PHASE DISTRIBUTIONS
close all;
figcount=1;
subplotrows=4;
subplotcols=5;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


assert(all(strcmp(OUTSTRUCT.bregionpair, 'LMAN-RA')), 'assumes chan1 is LAMN, chang 2 is RA...');
tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

% ---- params
ffbintoplot = 4;
tbintoplot = 11;
nrand = 10; % random trials to plot
nanglebins = 16;
ntoplot = 10; % random channels to plot
equalizeN = 1; % equalizes sample size by random sample from that with more samples.
plotraw = 1; % if 1, then takes random subset ot plot. if 0, then takes all.

% === targets
indslist = find(OUTSTRUCT.istarg==1)';
if plotraw==1
indslist = indslist(randperm(length(indslist), ntoplot));
else
end

% ===== to collect
PLV_all = nan(length(indslist), 2); % [base, wn]

for j=1:length(indslist)
    indthis = indslist(j);
    
    trialsbase = OUTSTRUCT.indsbase_epoch{indthis};
    trialsWN = OUTSTRUCT.indsWN_epoch{indthis};
    phimat = OUTSTRUCT.PhiMat{indthis};    
    maxtrials = min([length(trialsbase) length(trialsWN)]);
    
    
    % ================== BASELINE
    phithis = squeeze(phimat(:, ffbintoplot, trialsbase));
    if equalizeN==1
        phithis = phithis(:, randperm(size(phithis,2), maxtrials));
    end
    basewn_label = 'BASELINE';
    
    if plotraw==1
    % -- 1) plot a few trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('random trials');
    ylabel(basewn_label);
    indrand = randperm(size(phithis,2), nrand);
    plot(tbins, phithis(:, indrand), '-x');
    axis tight;
    
    % -- 2) plot all trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel('mean, std');
    plot(tbins, phithis, 'kx');
    plot(tbins, phithis+2*pi, 'kx');
%     phimean = mean(phithis,2);
%     phistd = std(phithis, [], 2);
%     lt_plot(tbins, phimean, {'Errors', phistd, 'Color', 'r'});
    axis tight
    ylim([-pi 2*pi])
    
    % -- 3) plot circular distribution of angles
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel(['at tbin:' num2str(tbintoplot)]);
    theta = phithis(tbintoplot, :);
%     theta(theta<0) = theta(theta<0)+2*pi; % so that all are positive
    rose(theta);
   
    % ---- 4) show vector sum
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    
    [x,y] = pol2cart(theta, ones(size(theta)));
    line([0 0], [mean(x) mean(y)], 'LineWidth', 3);
    xlim([-1 1]), ylim([-1 1]);
    grid on;
    
    
    % ---- get shuffle distribution of phase lock values.
    % DAT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('plv (shuff, dat)');
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    % SHUFF
    plv_shuff_all = [];
    nshuff = 100;
    ntrials = length(theta);
    for k=1:nshuff
        % -- sample uniformly from circle.
        theta_shuff = 2*pi*rand(size(theta))-pi;
        plv_shuff = lt_neural_QUICK_PhaseLockVal(theta_shuff);
        plv_shuff_all = [plv_shuff_all; plv_shuff];
    end
    lt_plot_histogram(plv_shuff_all);
    line([plv_dat plv_dat], ylim, 'Color', 'r');
    xlim([-0.05 0.7]);
    end
    
    % ========
    theta = phithis(tbintoplot, :);
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    PLV_all(j, 1) = plv_dat;
    
    
    % ================== WN
    phithis = squeeze(phimat(:, ffbintoplot, trialsWN));
    if equalizeN==1
        phithis = phithis(:, randperm(size(phithis,2), maxtrials));
    end
    basewn_label = 'WN';
    
    if plotraw==1
    % -- 1) plot a few trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('random trials');
    ylabel(basewn_label);
    indrand = randperm(size(phithis,2), nrand);
    plot(tbins, phithis(:, indrand), '-x');
    axis tight;
    
    % -- 2) plot all trials
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phitmp,2))]);
    ylabel('mean, std');
    plot(tbins, phithis, 'kx');
    plot(tbins, phithis+2*pi, 'kx');
%     phimean = mean(phithis,2);
%     phistd = std(phithis, [], 2);
%     lt_plot(tbins, phimean, {'Errors', phistd, 'Color', 'r'});
    axis tight
    ylim([-pi 2*pi])
    
    % -- 3) plot circular distribution of angles
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    ylabel(['at tbin:' num2str(tbintoplot)]);
    theta = phithis(tbintoplot, :);
%     theta(theta<0) = theta(theta<0)+2*pi; % so that all are positive
    rose(theta);
   
    % ---- 4) show vector sum
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all trials, N=' num2str(size(phithis,2))]);
    
    [x,y] = pol2cart(theta, ones(size(theta)));
    line([0 0], [mean(x) mean(y)], 'LineWidth', 3);
    xlim([-1 1]), ylim([-1 1]);
    grid on;
    
    
    % ---- get shuffle distribution of phase lock values.
    % DAT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('plv (shuff, dat)');
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    % SHUFF
    plv_shuff_all = [];
    nshuff = 100;
    ntrials = length(theta);
    for k=1:nshuff
        % -- sample uniformly from circle.
        theta_shuff = 2*pi*rand(size(theta))-pi;
        plv_shuff = lt_neural_QUICK_PhaseLockVal(theta_shuff);
        plv_shuff_all = [plv_shuff_all; plv_shuff];
    end
    lt_plot_histogram(plv_shuff_all);
    line([plv_dat plv_dat], ylim, 'Color', 'r');
    xlim([-0.05 0.7]);
    end
    
    theta = phithis(tbintoplot, :);
    plv_dat = lt_neural_QUICK_PhaseLockVal(theta);
    PLV_all(j, 2) = plv_dat;

end

lt_figure; hold on;
x = [1 2];
plot(x, PLV_all', '-ok');
lt_plot(x+0.2, mean(PLV_all,1), {'Errors', lt_sem(PLV_all), 'Color', 'r'});
xlim([0 3]);
xlabel('base - WN');
ylabel('PLV (equal N)');
signrank(PLV_all(:,1), PLV_all(:,2))

%% ===== PLOT LFP FOR MULTIPLE TRIALS OF ALL MOTIFS

close all;



%% ===== COMPUTE SPECTROGRAMS USING EXTRACTED

