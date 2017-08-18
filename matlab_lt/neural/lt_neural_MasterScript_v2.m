%% ACROSS NEURONS, FOR ANALYSES
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% EXTRACT 
clear all; close all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'LMAN', 'X'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
% BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'X'};
% ExptToKeep = {};
% RecordingDepth = [];
% LearningOnly = 1;
% BatchesDesired = {};
% ChannelsDesired = [];
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end


%% remove song dat from metadat

lt_neural_v2_PRE_RemvSongDat(SummaryStruct)

%% refinalize neurons
if (0)
lt_neural_v2_PRE_RefinalizeNeur
end

%% EXTRACT FF AND SAVE

lt_neural_v2_PRE_ExtrFF;

%% ===== EXTRACT WN HITS

% IN PROGRESS ! see inside
lt_neural_v2_EXTRACT_WNhit(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl);


%% ======== SAVE META INFO FOR LEARNING EXPT HERE
% edit this by hand

% 1) --- LOAD
LearningMetaDat = lt_neural_v2_LoadLearnMetadat;


% 2) --- EDIT
LearningMetaDat; % OPEN AND EDIT BY HAND. 
% Note: each expt/targ syl has one column:
% row 3 and larger gives time of switch and nature of switch (4 types of
% switches possible) (escape dir)
% Al = 100% WN
% Of = off
% Up = escape up
% Dn = escape dn


% 3) --- SAVE
currdir = pwd;
cd('/bluejay5/lucas/analyses/neural/');
save('LearningMetaDat', 'LearningMetaDat');
cd(currdir)


%% ==== LIST OF MOTIFS FOR EACH BIRD/EXPERIMENT
if (0) % RUNS AUTOMATICALLY WHEN EXTRACT
SummaryStruct = lt_neural_v2_PostInfo(SummaryStruct);
end



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% &&&&&&&&&&&&&&&&&& DIAGNOSTIC STUFF &&&&&&&&&&&&&&&&

%% ======= SONG MOD METRIC - FOR EACH NEURON




%% ============ Check pitch contours

close all;
useDiffColors = 1; % 1, plots each pc diff color; if 0, shades all gray
lt_neural_v2_DIAGN_pcontours(SummaryStruct, useDiffColors);
% TO DO: overlay WN (first extract WN and save in diff code)
% Get metric of how good PC is (eg fluctuation)


%% ============= SNR of smoothed firing rate 
% do this for each neuron, for each syl (separated by context)
% get distributions across birds, neurons, syls 
% assign metric to each neuron/syl
% compare by eye to FR. also to FR corr

% make code work for learning as well

% --- for example, pick anything
i=1;
ii=1;
nn=1;
mm =1;

segextract = ...
    MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;

% -- get FR matrix (time x trials)
clustnum = 1;
segextract = lt_neural_SmoothFR(segextract, clustnum);
FRmat = [segextract.FRsmooth_rate_CommonTrialDur];

t1 = 10;
t2 = 110;
trials = 1:20;
FRmat = FRmat(t1:t2, trials);

% --- RUN
[SNR, SignalPower, NoisePower] = lt_neural_v2_SNR(FRmat);



%% &&&&&&&&&&&&&&&&&&&&&& CONTEXT &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
clear CLASSES

% &&&&&&&&&&&&& 1) ARBITRARY CONTEXTS
strtype = 'xaa'; % a is fixed, x variable, across contexts
[CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype);

% &&&&&&&&&&&&& 2) EXTRACT REGEXP STRUCT 
prms.alignWhichSyl = 2; % which syl (in order) to align to
prms.alignOnset = 1; % if 1, then onset, if 0, then offset
prms.motifpredur = 0.2;
prms.motifpostdur = 0.3;
prms.preAndPostDurRelSameTimept = 1; % 1, then pre and post both aligned at same time. if 0, then post is aligned to motif ofset.
CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms);

% &&&&&&&&&&&&&& OPTIONAL - COLLECT POSITIVE CONTROL DATA
CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype);


% &&&&&&&&&&&&& 2) PLOT MEAN FR ACROSS CONTEXTS FOR EACH BRANCH 
close all;
plotPosControl = 1; % will do if exists.
lt_neural_v2_CTXT_FRanyclass(CLASSES, SummaryStruct, prms, plotPosControl);


% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& CLASSIFIER
% =========================== VERSION 1 (SINGLE ITERATION)
prms.ClassGeneral.frtimewindow =[-0.075 0.025]; % on and off, relative to syl onset
prms.ClassGeneral.frbinsize = 0.01; % in s.
prms.ClassGeneral.Nmin = 7; % in s.

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case

prms.ClassGeneral.GetPosControl =1;
CLASSES = lt_neural_v2_CTXT_ClassGeneral(CLASSES, SummaryStruct, prms);

% ============================ VERSION 2 - VARY PARAMETERS, MULTIPLE
% ITERATIONS
TimeWindowDur = 0.04;
TimeWindowSlide = 0.01;
ListOfTimeWindows = [-prms.motifpredur:TimeWindowSlide:prms.motifpostdur-TimeWindowDur; ...
    -prms.motifpredur+TimeWindowDur:TimeWindowSlide:prms.motifpostdur]'; % N x 2 (pre and post onset)
ListOfFrBinSizes = [0.005 0.01 0.02];
savenotes = 'allXLman';

prms.ClassGeneral.GetNegControl = 1; % 1 = yes. (i.e. shuffle dat-context link).
prms.ClassGeneral.GetNegControl_N = 20; % number iterations for each case

prms.ClassGeneral.GetPosControl =1;
[savedir] = lt_neural_v2_CTXT_ClassGeneral_M(CLASSES, SummaryStruct, prms, ListOfTimeWindows, ...
    ListOfFrBinSizes, savenotes);

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& PLOTS
% ========== VERSION 1 - PLOT SINGLE ITERATION
close all;
lt_neural_v2_CTXT_PlotGeneral(CLASSES, SummaryStruct, prms);

% ========== VERSION 2- PLOT ALL
strtype = 'axa';
algnsyl = 2;
algnonset = 1;
CLASSEScompiled = lt_neural_v2_CTXT_PlotGeneral_M(strtype, algnsyl, algnonset);

plotstat = 'pihatmean';
lt_neural_v2_CTXT_PlotGeneral_M2(CLASSEScompiled, plotstat);


% ========== VERSION 3- COMBINE ALL PLOTS FOR A GIVEN MOTIF
% combine across 1) syl aligned to 2) aligned to onset or offset, 3) fr
% window size 4) fr bin size
strtype = 'axa';



%% ============ 1) CLASSIFY CONTEXT USING NEURAL

% ========== 1) EXTRACT DATA 
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
LearnKeepOnlyBase = 1; 
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase);


% ====== 2) PLOT FR TRACES
lt_neural_v2_ContextFR(MOTIFSTATS_Compiled);



% ========== 2) CLASSIFY CONTEXT
nmin = 5; % for each context. skips overwise
CLASSIFIEROUT = lt_neural_v2_CTXT_Class(MOTIFSTATS_Compiled, SummaryStruct, nmin);


% ========== 3) PLOT CLASSIFIER OUTPUT
lt_neural_v2_CTXT_Plot(CLASSIFIEROUT);



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% &&&&&&&&&&&&&&&&&& VOCLATIONS --> predict neural &&&&&&&&&&&&&&&&&&

close all;
plotRaw = 0; % individaul neuronf igs (note, hand entering neuron ind currently, in code)
% Binparams.Pretime = 0.2; % to start getting binned data (rel to onset)
% Binparams.Posttime = 0.2; % to stop getting dat (rel to onset);
% Binparams.Binsize = 0.025; % for getting spike counts
Binparams.Pretime = 0.7; % to start getting binned data (rel to onset)
Binparams.Posttime = 0.7; % to stop getting dat (rel to onset);
Binparams.Binsize = 0.025; % for getting spike counts

ConvOrDiv = 'conv';
saveOn = 0;

VOCALSTRUCTall = lt_neural_v2_ANALY_VocModel(SummaryStruct, Binparams, plotRaw, ConvOrDiv, saveOn);


% ===== 1) FILTER VOCAL STRUCT


% ===== 2) PLOT
close all;
birdnum = 1;
neuronnum=1;
timebin = 35;
plotSmoothed = 0; 
lt_neural_v2_ANALY_VocModel_plot(VOCALSTRUCTall, birdnum, neuronnum, timebin, plotSmoothed)


% ===== 3) ANOVA (done for each time bin)
% currerntly good. apply filter preceding to extract,e.g.. only one branch
% point, and so on. 
lt_neural_v2_ANALY_VocModel_anova(VOCALSTRUCTall)


% ====== 4) 
lt_neural_v2_ANALY_VocModel_glm(VOCALSTRUCTall)



%% ===== EFFECT OF MOTIF POSITION (WITHIN SONG BOUT) ON SYL FEATURES + NEURAL [JOANNE]

% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
close all;
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);




%% ===== COMPARE RAW NEURAL (FILTERED) FOR ALL RENDS FOR SEGMENT

% --- IN PROGRESS
lt_neural_v2_ANALY_PlotRawNeural(SegmentsExtract);




%% &&&&&&&&&&&&&&&&&&&&& SINGLE EXPERIMENTS &&&&&&&&&&&&&&&&&&&&&&&&&
%% ============== FOR EACH LEARNING EXPERIMENT, PLOT ALL NEURONS AND ALL MOTIFS
% ONE BIRD/LEARNING EXPT AT A TIME [BUT CAN COMBINE NEURONS]

% ==== 1) EXTRACT DATA FOR EACH NEURON, EACH MOTIF
[MOTIFSTATS, SummaryStruct_filt] = lt_neural_v2_ANALY_LearningExtractMotif(SummaryStruct);
exptname = SummaryStruct_filt.birds(1).neurons(1).exptID;
birdname = SummaryStruct_filt.birds(1).birdname;
% ---- EXTRACT TARG SYL
tmp = lt_neural_v2_LoadLearnMetadat;
indbird = strcmp({tmp.bird.birdname}, birdname);
indexpt = strcmp(tmp.bird(indbird).info(1,:), exptname);
TargSyls = tmp.bird(indbird).info(2,indexpt);
MOTIFSTATS.params.TargSyls = TargSyls;

close all;
% motifs for this bird
% learning expt id
plottype = 'bysyl'; %
% 'byneuron' - each neuron one fig will all motifs [DEFAULT]
% 'dotprod' - for each bin of trials get dot prod from IN PROGRESS
% 'bysyl' - each plot one syl, all neurons.
DivisorBaseSongs = 1;
lt_neural_v2_ANALY_Learning(SummaryStruct_filt, MOTIFSTATS, plottype, DivisorBaseSongs);




%% ========= WN RESPONSE?



%% &&&&&&&&&&&&&&&&&& SUMMARY STATISTICS ACROSS ALL EXPERIMENTS &&&&&&&&&&&&&&&&
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ================ META SUMMARY OF LEARNING
% NOTE: MODIFY TO HAVE GROSS SMOOTHED FR PLOTTED AS WELL
% TO DO - indicate for each neuron its "quality" (e.g. song mod)

close all;
lt_neural_v2_LEARNING_MetaSummary

%% =================== LEARNING 
% CAN DO MULTIPLE BIRDS
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit);


% =========== PICK OUT SAME TYPE/DIFF [LEANRING SPECIFIC]
numbirds = length(MOTIFSTATS_Compiled.birds);
for i=1:numbirds
    birdname = MOTIFSTATS_Compiled.birds(i).birdname;
    numexpts = length(MOTIFSTATS_Compiled.birds(i).exptnum);
   
    for ii=1:numexpts
       exptname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).exptname;
        % ---- EXTRACT TARG SYL
        tmp = lt_neural_v2_LoadLearnMetadat;
        indbird = strcmp({tmp.bird.birdname}, birdname);
        indexpt = strcmp(tmp.bird(indbird).info(1,:), exptname);
        TargSyls = tmp.bird(indbird).info(2,indexpt);

        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.TargSyls = TargSyls;
        
        % ----
%         TargSyls = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.TargSyls;
%         MotifsActual = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.MotifsActual;
        motif_regexpr_str = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.motif_regexpr_str;
        
        [SameTypeSyls, DiffTypeSyls, motif_regexpr_str, SingleSyls] = ...
        lt_neural_v2_extractSameType(motif_regexpr_str, TargSyls);

    % --- OUTPUT
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.SameTypeSyls = SameTypeSyls;
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.DiffTypeSyls = DiffTypeSyls;
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.SingleSyls_inorder = SingleSyls;

    end
end



% ==== FOR EACH LEARNING EXPERIMENT, PLOT TIMELINE OF NEURONS AND LEARNING
close all;
NumBirds = length(SummaryStruct.birds);
for i=1:NumBirds
   ListOfExpts = unique({SummaryStruct.birds(i).neurons.exptID});
   
   for ll = 1:length(ListOfExpts)
        
       MOTIFSTATS = MOTIFSTATS_Compiled.birds(i).exptnum(ll).MOTIFSTATS;
       SummaryStruct_tmp = MOTIFSTATS_Compiled.birds(i).exptnum(ll).SummaryStruct;
       
       birdname = SummaryStruct.birds(i).birdname;
       exptname = ListOfExpts{ll};
       
       if ~isfield(MOTIFSTATS.params, 'TargSyls')
           disp('asdasdf');
       end
       if isempty(MOTIFSTATS.params.TargSyls);
           sdafasdf
       end
       
       
      % === PLOT
      lt_neural_v2_ANALY_LearningPlot1(SummaryStruct_tmp, MOTIFSTATS);
       title([birdname '-' exptname]);
   end
end
   

%% ==== FOR EACH LEARNING EXPERIMENT, PLOT SUMMARY STATISTICS
close all;

% ***************************************************** VERSEION 1 (SORT
% BY NERUON)
% ===================== EXTRACTS LEARNING STATS
% IMPORTANT - USE INTERNAL SUMMARYSTRUCT (INTERNAL TO MOTIFSTATS_Compiled)
% and not actual full SummaryStruct (neuron inds don't match)
[MOTIFSTATS_Compiled] = lt_neural_v2_ANALY_LearningStats(MOTIFSTATS_Compiled);


% =========== PLOT TIMECOURSES
birdtoplot = 'wh6pk36'; % leave as '' to get all
expttoplot = 'LMANlearn2';
neuralmetric_toplot = 'NEURrelbase_smthFrCorr';
%         NEURrelbase_smthFrCorr
%         NEURrelbase_EuclDistance
%         NEURrelbase_Norm1Distance
%         NEURrelbase_MeanFRDiff
%           NEURrelbase_EuclDistFRcentered
%           NEUR_CVsmthFR
%           NEUR_SDsmthFR
%           NEUR_MeansmthFR
plotOnlyTargSyl = 1; %

lt_neural_v2_ANALY_LearningPLOTtcour(MOTIFSTATS_Compiled, birdtoplot, ...
    expttoplot, neuralmetric_toplot, plotOnlyTargSyl)


% =========== PLOT SUMMARY OF STATS ACROSS EXPERIMENTS [ITERATING OVER
% NEURONS
close all;
neuralmetric_toplot = 'NEURrelbase_smthFrCorr';
convertToZscore = 0;
PlotNeuralFFDeviationCorrs = 1; % zscore of neural similarity and FF, correalted fore aach syl
plotOnlyTargSyl = 1; % only matters if PlotNeuralFFDeviationCorrs=1
birdtoplot = 'bk7'; % leave as '' to get all
expttoplot = 'LearnLMAN1';
lt_neural_v2_ANALY_LearningStatsPLOT(MOTIFSTATS_Compiled, convertToZscore, ...
    neuralmetric_toplot, PlotNeuralFFDeviationCorrs, plotOnlyTargSyl, ...
    birdtoplot, expttoplot)




% ******************************************************* VERSION 2 - SWITCH AS DATAPOINT
% NOTE: neuron inds in switchstruct match those in motifs_compiled

% === PULL OUT RAW FR FOR ALL NEURONS/TRIALS
RemoveTrialsZeroFR = 1;
premotorWind = [-0.075 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset
[MOTIFSTATS_Compiled] = lt_neural_v2_ANALY_GetAllFR(MOTIFSTATS_Compiled, ...
    RemoveTrialsZeroFR, premotorWind);


% PULL OUT "SWITCHES" IN CONTINGENCY ACROSS ALLEXPERIMENTS [alsop extracts
% FF values at switch]
[MOTIFSTATS_Compiled, SwitchStruct] = lt_neural_v2_ANALY_GetSwitches(MOTIFSTATS_Compiled);


% ============ EXTRACT NEURAL FOR SWITCHES
interpolateCorrNan = 0;
RemoveTrialsZeroFR = 0; % this takes precedence over interpolateCorrNan
[MOTIFSTATS_Compiled, SwitchStruct] = ...
    lt_neural_v2_ANALY_Swtch_Extract(MOTIFSTATS_Compiled, SwitchStruct, ...
    interpolateCorrNan, RemoveTrialsZeroFR, premotorWind);


% =========== PLOT SUMMARY OF LEARNING [ITERATING SWITCHES]
close all;
lt_neural_v2_ANALY_LrnSwtchPLOT(MOTIFSTATS_Compiled, SwitchStruct);


% ============ TIMECOURSES FOR NEURAL FOR SWITCHES
close all;
birdname_get = 'br92br54'; % keep empty if want all.
exptname_get = 'LMANlearn7';
switchnum_get = [2];
plotneurzscore=0;
lt_neural_v2_ANALY_Swtch_Tcourse(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore)



% ============== SUMMARIZE FOR EACH EXPT (I.E. NEURAL AND FF CHANGE FOR
% TARG AND OTHER SYLS)
close all;
RemoveLowNumtrials = 1; % min number for both base and train(double)
MinTrials = 10; % for removing
skipMultiDir = 0;
usePeakLearn = 1; % assumes that care about learning for targ 1. (median time of max learning)

% ---- what metric to plot for main figs
fieldname_baseneur = 'AllNeurSplitCorrBase';
fieldname_trainneur = 'AllNeurSplitCorrTrain';
% fieldname_baseneur = 'AllNeurCorrShuffMean';
% fieldname_trainneur = 'AllNeurCorrDat';
% fieldname_baseneur = 'AllNeurSimBase';
% fieldname_trainneur = 'AllNeurSimTrain';

% -- following only matter if use fieldname_baseneur = 'AllNeurSimBase' & fieldname_trainneur = 'AllNeurSimTrain';
UseZscoreNeural = 0;
neuralmetricname = 'NEURvsbase_FRcorr';
% neuralmetricname = 'NEUR_meanFR';
% neuralmetricname = 'NEUR_cvFR';

% ---- filtering data by learning at target
% note: if multiple targets then will filter on just first target...
plotLearnStatsOn = 0; % 
learnsigalpha = 0.15; % for deciding that an experiment showed "no learning"
OnlyKeepSigLearn = 1; % either regression or end training has be significant.

% ---- ONLY KEEP SWITCHES STARTING FROM WN OFF
OnlyKeepWNonset =0; % if 1, then yes, if 2, then only keeps those with WN transition (not onset); if 0, then takes all

[DATSylMot, ~] = lt_neural_v2_ANALY_Swtch_Summary(MOTIFSTATS_Compiled, SwitchStruct, RemoveLowNumtrials, ...
    MinTrials, UseZscoreNeural, neuralmetricname, fieldname_baseneur, fieldname_trainneur, ...
    skipMultiDir, usePeakLearn, plotLearnStatsOn, learnsigalpha, OnlyKeepSigLearn, ...
    OnlyKeepWNonset);




% ================ IN LINEAR MODEL IS THERE EFFECT OF SYLLABLE TYPE AFTER
% CONTROLLING FOR VARIOUS THINGS?
lt_neural_v2_ANALY_Swtch_LME(DATSylMot)





