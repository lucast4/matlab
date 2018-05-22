%% ACROSS NEURONS, FOR ANALYSES
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% EXTRACT
clear all; close all; fclose all;
BirdsToKeep = {}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
BrainArea = {'RA', 'LMAN', 'X'}; % if want Sam/Mel data, must include "RAmel"
% ExptToKeep = {'RAlearn1', 'RALMANlearn1', 'LMANsearch'};
ExptToKeep = {};
RecordingDepth = [];
LearningOnly = 0;
BatchesDesired = {};
ChannelsDesired = [];
extractpreDatenums = 1;
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired, ...
    extractpreDatenums);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end

%% ############################################### DISP (stuff)
%% ############################################################
%% get bird and neuron number for desired neruon

dirname = '/bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2/Chan14amp-Batch1058to1106';
birdname = 'wh44wh39';

% =========== RUN
birdnum = find(strcmp({SummaryStruct.birds.birdname}, 'wh44wh39'))
neurid = find(strcmp({SummaryStruct.birds(birdnum).neurons.dirname}, dirname))




%% ========== plot all units info

lt_neural_DISP_AllUnits(SummaryStruct);

lt_neural_DISP_AllPopUnits(SummaryStruct);

%% check fs for all

numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        
        disp([num2str(i) '-' num2str(ii)]);
        cd(SummaryStruct.birds(i).neurons(ii).dirname)
        tmp = load('times_data.mat');
        tmp2 =load('MetaDat.mat');
        
        assert(unique([tmp2.metaDat.fs]) == tmp.par.sr, 'problem, fs of sopng and neural not equal');
        
    end
    
end

%% take snapshot of current raw data - i.e. label files, extracted FF

[outdir] = lt_neural_SnapshotCurrDat(SummaryStruct);

%% remove song dat from metadat

lt_neural_v2_PRE_RemvSongDat(SummaryStruct)

%% refinalize neurons
if (0)
    lt_neural_v2_PRE_RefinalizeNeur
end

%% EXTRACT FF AND SAVE
close all;
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
% Al = 100% WN [NOTE: is currently coded same as Of (for switch struct extraction)
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
%% ====== CONFIRM THAT ALL CLUSTERING HAVE "FORCED"
% i.e. nonclustered spikes have been forced to a spike class (or to noise
% class)

NotForced = []; % n x [bird neur] 
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        disp(SummaryStruct.birds(i).neurons(ii).dirname);
        cd(SummaryStruct.birds(i).neurons(ii).dirname);
        tmp = load('times_data.mat');
        if ~any(tmp.forced==1)
            % note this down
            NotForced = [NotForced; [i ii]];
        end
        
    end
    
end

if isempty(NotForced)
    disp('Great - all forced!');
else
    disp('Not good - some not forced: ');
disp(NotForced);
end

%% ====== check if SU or MU


%% ======= CHECK IF GOOD FR MODULATION



%% ============= CHECK WHETHER ANY UNITS HAVE OVERLAPPING DATA - IF SO, REMOVE ONE OF THEM

[SummaryStruct, NeurToRemove] =lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct);


%% ============ DISPLAY SONGS AND LABELS
stoponbird = 1; % if 1, then pauses and closes; if 0 then not;
lt_neural_v2_DIAGN_DispLabels(SummaryStruct, stoponbird);


%% ===== CHECK CLUSTERING QUALITY - compare filtered neural with cluster
% also plot other channels if want to compare
close all;
displaymode = 'rand';
skipnum = 5; % only matters if is skip mode
% Birdname = 
% Exptname
lt_neural_v2_DIAGN_Rawdat(SummaryStruct, displaymode, skipnum)


%% ====== PLOT RAW NEURAL FOR ANY GIVEN BIRD, EXPT, MOTIF

close all; 
BirdToPlot = 'O55Pu53';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = [6]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
% motiflist = {'a(b)', 'jbh(h)g'};
% motiflist = {'(d)kcc', 'dk(c)c', '(n)hh', 'c(b)'};
motiflist = {'a(a)b', 'i(a)b', 'c(a)b'};

% motifpredur = 0.15;
% motifpostdur = 0.15;
motifpredur = 0.15;
motifpostdur = 0.1;
preAndPostDurRelSameTimept = 1;

% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

saveON = 0; 

lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
    motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON);


%% ======= [MODIFICATION] plot raw neural dat, all main motifs, and overlaying
close all; fclose all;

% == params, entries aligned
BirdsToPlotList = {'pu69wh78', 'wh44wh39'};
NeurToPlotList = {[], []}; % leave empty if plot all
MotifsToPlot = {{'(j)jjbhhg', '(a)abhhg'}, {'(n)hh', '(d)kccbb'}};

% ======== general params
motifpredur = 0.05;
motifpostdur = 0.05;
preAndPostDurRelSameTimept = 0;
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both
saveON = 1;


for i=1:length(BirdsToPlotList)
    BirdToPlot = BirdsToPlotList{i};
    NeurToPlot = NeurToPlotList{i};
    motiflist = MotifsToPlot{i};
    
    lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
        motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON);
end

%% ======= SONG MOD METRIC - FOR EACH NEURON




%% ============ Check pitch contours

close all;
useDiffColors = 1; % 1, plots each pc diff color; if 0, shades all gray
plotbysyl = 0; % if 1, then each plot one syl. if 0, then go by neuron. [IN PROGRESS]
dolinkaxes = 0;
lt_neural_v2_DIAGN_pcontours(SummaryStruct, useDiffColors, plotbysyl, dolinkaxes);
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


%% ================================== PLOT RASTER AND SMOOTHED FR FOR ANY MOTIF
close all
BirdToPlot = 'pu26y2';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = [7]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
% motiflist = {'a(b)', 'jbh(h)g'};
motiflist = {'i(a)c', 'a(a)c', 'e(a)c'};
plotbytime = 0; % links rasters for all motifs by time of song.

motifpredur = 0.1;
motifpostdur = 0.1;

plotIndivRaster = 0; % one raster for each neuron/motif
plotCombRast = 1; % one figure, all rasters
plotSmFR = 1; % all smoothed FR.

% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong )


if (0)
% --- 2) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong )
end


%% ===== COMPARE RAW NEURAL (FILTERED) FOR ALL RENDS FOR SEGMENT

% --- IN PROGRESS
lt_neural_v2_ANALY_PlotRawNeural(SegmentsExtract);



%% ################# DIR vs UNDIR ###
%% ##################################################################

lt_neural_v2_DirUndir_Master;

%% ################# LEARNING (SINGLE CHANNEL)  ###
%% ##################################################################

lt_neural_v2_LEARN_Master;


%% ################# ENCODING OF SONG PARAMS (OTHER THAN CONTEXT) ###
%% ##################################################################

%% ======= CORRELATION BETWEEN SINGLE UNITS AND VARIOUS PARAMS




%% ===== EFFECT OF MOTIF POSITION (WITHIN SONG BOUT) ON SYL FEATURES + NEURAL [JOANNE]

% === IMPORTANT!! - ASSUMES THAT ALL NEURONS FOR A GIVEN BIRD HAVE SAME
% MOTIFS IN SAME ORDER (I.E. IN lt_neural_v2_PostInfo)
close all;
PlotRaw = 0;
lt_neural_v2_ANALY_BoutPositionJC(SummaryStruct,PlotRaw);






%% &&&&&&&&&&&&&&&&&&&&&& NGRAMS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% like "CONTEXT" but takes all ngrams for a neuron to get distrubitions for
% each neruon

lt_neural_NGRAMS_MasterScript;



%% &&&&&&&&&&&&&&&&&&&&&& CONTEXT &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

lt_neural_v2_CTXT_MASTER;


%% ============ 1) CLASSIFY CONTEXT USING NEURAL
close all;
% ========== 1) EXTRACT DATA
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
LearnKeepOnlyBase = 1;
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase);


% ====== 2) PLOT FR TRACES
close all;
plotSTD =0;
doplot_bymotif = 0;
doplot_bysinglesyl = 0;
doplot_LMANRA = 1;
lt_neural_v2_ContextFR(MOTIFSTATS_Compiled, plotSTD, doplot_bymotif, doplot_bysinglesyl, ...
    doplot_LMANRA);



% ========== 2) CLASSIFY CONTEXT
nmin = 5; % for each context. skips overwise
CLASSIFIEROUT = lt_neural_v2_CTXT_Class(MOTIFSTATS_Compiled, SummaryStruct, nmin);


% ========== 3) PLOT CLASSIFIER OUTPUT
lt_neural_v2_CTXT_Plot(CLASSIFIEROUT);


%% %%%%%%%%%%%%%%%%%%%%%%% CORRELATIONS BETWEEN SYLLABLES

% ========== 1) EXTRACT DATA
collectWNhit=0; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
LearnKeepOnlyBase = 1;
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase);


%% get pairwise stats [between syllables]

% TO DO; need to define motifs and note down whether is on same motif, etc.

i = 1;
ii=1;
nn=1;

nummotifs = length(MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif);

for m = 1:nummotifs
    
    motifdat1 = ...
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(m);
    
    for mm = m+1:nummotifs
        
        motifdat2 = ...
            MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif(mm);
        
        
        % =============== COMPUTE THINGS FOR THIS PAIR OF SYLLABLES
        
        % ---------------- SAME TYPE?
        
        
        % ---------------- SAME MOTIF? DISTANCE ON MOTIF?
        
        
        % ---------------- PITCH (mean and corr)
        
        
        % ---------------- SPIKE COUNT, CORR ETC
        [motifdat1.SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]
        [motifdat2.SegmentsExtract.global_tokenind_DatAlignedToOnsetOfThis]
        
        
        
        
        
    end
end

%%  get pairwise stats [between neurons]



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% &&&&&&&&& VOCLATIONS --> predict neural using context/FF, etc &&&&&&&&&&&&&&&&&&

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





%% &&&&&&&&&&&&&&&&&&&&&& POPULATION LEARNING &&&&&&&&&&&&&&&&&&&

%% ==== 1) EXTRACT MOTIFSTATS POP

% ======================== EXTRACT SEGMENTS FOR POPULATIONS

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


%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);

%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;


%% ================ PLOT [CORRELATION WITH FF]
close all;
xcov_dattotake = [-0.01 0.05];
xcov_dattotake = [-0.08 0.04];
xcov_dattotake = [-0.075 0.025];
xcovwindmax = 0.04;
binsize_spk = 0.0025;

MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk);

 %% ==== 2) EXTRACT LEARNING SWITCH STRUCT

SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);

%% ================ PLOT CROSS CORR WRT TO LEARNING
close all; 
BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn4'};
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


%% =========== SUMMARIZE LEARNING TRAJECTORY (PLUS NEURON SETS)
close all;

% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn1'};
% motiftoplot = 'c(b)';
BirdExptPairsToPlot = {'pu69wh78', 'RALMANOvernightLearn1'};
motiftoplot = 'aa(b)';

lt_neural_POPLEARN_PlotLearnTraj(MOTIFSTATS_pop ,SwitchStruct, ...
    SummaryStruct, BirdExptPairsToPlot, motiftoplot);


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
BirdExptPairsToPlot = {};
SwitchToPlot = [2];
TypeOfPairToPlot = {'LMAN-RA'}; % e.g. 'LMAN-RA' (in alphabetical order)

numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    
    numexpts = length(SwitchStruct.bird(i).exptnum);
    birdname = SwitchStruct.bird(i).birdname;
    
    for ii=1:numexpts
        
        numswitches = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        % ----------------- ONLY PLOT SPECIFIC BIRD?
        if ~isempty(BirdExptPairsToPlot)
           
            ind1 = find(strcmp(BirdExptPairsToPlot, birdname));
            ind2 = find(strcmp(BirdExptPairsToPlot, exptname));
            
            if ~any(ind1+1 == ind2)
                disp(['SKIPPED ' birdname '-' exptname]);
                continue
            end
            
        end
        
        % ----------------- GO THRU ALL SWITCHES
        for iii=1:numswitches
            
            if ~isempty(SwitchToPlot)
               if ~any(SwitchToPlot == iii)
                   continue
               end
            end
            
            % ---- for this switch, figure out which populations have data
            % overlapping the onset (i.e. has data both pre and post swictch)
            swthis = SwitchStruct.bird(i).exptnum(ii).switchlist(iii);
            numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);
            
            for ss = 1:numsets
                songfiles = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_songfiles{ss};
                songtimes = datenum(songfiles, 'yymmdd_HHMMSS');
                
                inds_pre = find(songtimes>swthis.switchdnum_previous & songtimes<swthis.switchdnum);
                inds_post = find(songtimes>swthis.switchdnum & songtimes<swthis.switchdnum_next);
                
                if isempty(inds_pre) | isempty(inds_post)
                    continue
                else
                    disp(['analyzing: ' birdname '-' exptname '-sw' num2str(iii) '-neurset' num2str(ss)]);
                end
                
                
                
                % ############################################### ANALYSIS/PLOTS
                DAT = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss);
                motiflist = {DAT.motif.regexpstr};
                neurlist = MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons{ss};
                
                % ============ for each pair of neurons, plot paired
                % rasters
                % -- go thru all pairs of neurons, only plot if is desired
                % type of pair
                for j=1:length(neurlist)
                    for jj=j+1:length(neurlist)
                   
                        n1 = neurlist(j);
                        n2 = neurlist(jj);
                        
                        % ----- check what pair of brain region
                        
                        
                    end
                end
                
                
            end
        end
    end
end


%% ###################################### DIR VS. UNDIR (LEARNING)

%% ==== 1) EXTRACT MOTIFSTATS POP

% ======================== EXTRACT SEGMENTS FOR POPULATIONS
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



%% ==== 2) EXTRACT LEARNING SWITCH STRUCT

SwitchStruct = lt_neural_LEARN_getswitch(SummaryStruct);


%% ==== 3) PLOT (for each switch, plot all neurons and syls)
close all;
% BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn3'};
BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn3'};
SwitchToPlot = [2];
plotIndTrial = 0;

lt_neural_v2_DirUndir_LearnPlot(MOTIFSTATS_Compiled, SwitchStruct, BirdExptPairsToPlot, ...
    SwitchToPlot, plotIndTrial, SummaryStruct);

%% ==== 4) PLOT (for each motif, plot across all switches for a 
% given experiment
close all;
BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn3'};
% SwitchToPlot = [1];
plotIndTrial = 0;
xwindplot = [-0.1 0.05];
onlyPlotIfBothPrePostTrials = 0;
motiftoplot = {'dk(c)c'};

lt_neural_v2_DirUndir_LearnPlot2(MOTIFSTATS_Compiled, SwitchStruct, BirdExptPairsToPlot, ...
    plotIndTrial, SummaryStruct, motiftoplot, xwindplot, onlyPlotIfBothPrePostTrials)

%% =====  SUMMARY PLOT OF TIMECOURSE (all syllables and units)
close all;
BirdExptPairsToPlot = {'wh44wh39', 'RALMANlearn1'};
xwindplot = [-0.055 0.035];
matchDirUndirTrials = 1;
usemedian = 1; % instead of mean (of FR across trials, before do corr)

lt_neural_v2_DirUndir_LearnSummary1(MOTIFSTATS_Compiled, SwitchStruct, BirdExptPairsToPlot, ...
    SummaryStruct, xwindplot, matchDirUndirTrials, usemedian);


%% ###############################################################
%% ############################################  POPULATION

lt_neural_v2_POP_Master;


%% ################################################################
%% ########################################### WN RESPONSE

% ========= for eac bird/expt, plot all trials separated by WN/noWN
close all; clear MOTIFSTATS_Compiled;
collectWNhit=1; % NOTE!!! temporary - need to change so don't need to extract audio each time (i.e. do and save)
onlyCollectTargSyl=1;
LearnKeepOnlyBase = 0;
saveOn = 0;
MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, LearnKeepOnlyBase, saveOn, onlyCollectTargSyl);


% 
% lt_neural_v2_PseudoPop_Master



%% #################################################################
%% #################################################################
%% PSEUDO POPULATION

lt_neural_PseudoPop_Master;


%% #############################################################
%% ############################ RA vs LMAN (mean firing)


lt_neural_v2_FullMotifActivity;

