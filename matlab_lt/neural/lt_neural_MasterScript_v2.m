%% ACROSS NEURONS, FOR ANALYSES
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% EXTRACT
clear all; close all; fclose all;
BirdsToKeep = {'pu69wh78', 'wh44wh39'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BirdsToKeep = {'wh72pk12'}; % {birdname , neuronstokeep} if neuronstokeep = [], then gets all;
% BrainArea = {'LMAN', 'RA', 'LMANoutside', 'RAoutside'}; % IF DOING NEGATIVE CONTROLS.
BrainArea = {'LMAN', 'RA'}; % 
% BrainArea = {}; % if want Sam/Mel data, must include "RAmel"
ExptToKeep = {'RALMANOvernightLearn1', 'RALMANlearn1', 'RALMANlearn3'};
RecordingDepth = [];
LearningOnly = 1;
BatchesDesired = {'Batch0922to1529', 'Batch0922to1300UndirDir', 'Batch1435to1629DirUndir', ...
    'Batch1435to1719DirUndir', 'Batch0957to2238UndirDir', 'Batch1016to1301',  'Batch1434to2030UndirDir', ...
    'Batch1326to1740DirUndir', 'Batch1743to2319', 'Batch1345to2242DirUndir', ...
    'Batch1033to1427DirUndir', 'Batch1550to2036DirUndir', 'Batch1126to2105', ...
    'Batch1126to1800'};
ChannelsDesired = [];
extractpreDatenums = 1;
onlySpikes = 1; % [default =1]if 1, then only keeps if have spiking data; if 0, then gets anything.
% if 2, then is for LFP analysis - i.e. for wh72 gets only if not SU or MU
% set to 2 if want to get data for cohere/xcov.
% NOTE: use 1 if want to look at spiking experiments across birds...
[NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database(BirdsToKeep, ...
    BrainArea, ExptToKeep, RecordingDepth, LearningOnly, BatchesDesired, ChannelsDesired, ...
    extractpreDatenums, onlySpikes);

% --- load all neurons
if (0)
    [NeuronDatabase, SummaryStruct] = lt_neural_v2_ConvertSummary2Database;
end


%% ===== PARAMS TO USE FOR DIFFERENT DATASETS 



%% ==== LIST OF PREPROCESSING STEPS TO DO WHEN ADD NEW NEURONS
% 1) remove song dat
% 2) extract FF (compare to birds own analysis)
% 3) PostInfo (motifs)
% 4) learning information
% 5) get LFP
% 6) get WN
% 7) plot example motif to see if stable.
% 8) syllables to skip if noisy, WN overlap etc


%% ========= [NEGATIVE CONTROLS] CHANNELS OUTSIDE OF SONG SYSTEM
% RUN THIS AFTER EXTRACTING. % USE ESPECIALYL FOR POPLEARN ANALYSIS.

SummaryStruct = lt_neural_CONTROL_Extract(SummaryStruct);

%% ############################################### DISP (stuff)
%% ############################################################
%% get bird and neuron number for desired neruon

dirname = '/bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2/Chan14amp-Batch1058to1106';
birdname = 'wh44wh39';

% =========== RUN
birdnum = find(strcmp({SummaryStruct.birds.birdname}, 'wh44wh39'))
neurid = find(strcmp({SummaryStruct.birds(birdnum).neurons.dirname}, dirname))




%% ========== plot all units info
close all;
lt_neural_DISP_AllUnits(SummaryStruct);

lt_neural_DISP_AllPopUnits(SummaryStruct);


%% ============= [SHOW METADAT] FOR EACH BIRD AND NEURON, SHOW DIRECTORY
% for each directory, displays channels to get.
longversion = 1; % then shows each extracted neuron
lt_neural_SummaryStruct_ShowDirChans(SummaryStruct, longversion);

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

%% ==================== 1) Get list of all song files and correspond chans.
[Allbird_Fnames, Allbird_chanlist, Allbird_birdnum] = lt_neural_tools_allsongs(SummaryStruct);


%% *******************************************************************
%% ***************************************** EXTRACTION CODE
%%  ==== GET COHERENCE AND SPECTROGRAMS OF LFP FOR ALL SONGS
% NOTE: this is not current - current is to first get LFP, and then do
% stuff with the LFP...
close all;
lt_neural_Coher_Extract(SummaryStruct);


%% ############ GET LFP FOR ALL SONGS

lt_neural_tools_LFPextract(SummaryStruct);

%% ############# FILTER BANK OVER ALL LFP
% === 1) DESIGN FILTER BANK
freqvals = [5:150];
jai_waveletfilt_DesignFilt(1500, freqvals)

% ---- extracts "spectrogram" output using filter bank.
freqvals = [10:5:130];
skipifdone =1;
lt_neural_tools_LFPfiltbank(SummaryStruct, freqvals, skipifdone);

%% ==== EXTRACTION - WITHIN SONG TIMING OF EACH IND
% NOTE: THIS IS OLD VERSION, new version directly extracts when does RegExp
% analysis.
numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    numneurons = length(SummaryStruct.birds(i).neurons);
    for ii=1:numneurons
        lt_neural_v2_EXTRACT_WithinSongTimings(SummaryStruct, i, ii);
    end
end

%% EXTRACT FF AND SAVE
close all;
lt_neural_v2_PRE_ExtrFF;


%% ===== EXTRACT WN HITS

if (0)
% IN PROGRESS ! see inside
lt_neural_v2_EXTRACT_WNhit(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl);
end

% === new version, matches up with original song files
% see lt_generalization...
[Allbird_Fnames, Allbird_chanlist, Allbird_birdnum] = lt_neural_tools_allsongs(SummaryStruct);
lt_batchsong_extractWNhit([], [], Allbird_Fnames);




%% *******************************************************************
%% ***************************************** METADATA CODE
%% ======== SAVE META INFO FOR LEARNING EXPT HERE
% edit this by hand
% NOTE: make sure that the first switch for every syllable starts with Of
% [if it is true that there is period withouit WN that you would want to
% analyze]

% 1) --- LOAD
LearningMetaDat = lt_neural_v2_LoadLearnMetadat;


% 2) --- EDIT
LearningMetaDat; % OPEN AND EDIT BY HAND.
% Note: each expt/targ syl has one column:
% row 3 and larger gives time of switch and nature of switch (4 types of
% switches possible) (escape dir)
% Al = 100% WN [NOTE: is currently coded same as Of (for switch struct extraction)
% [note: in some cases have random WN - also call Al]
% Of = off
% Up = escape up
% Dn = escape dn


% 3) --- SAVE
currdir = pwd;
cd('/bluejay5/lucas/analyses/neural/');
save('LearningMetaDat', 'LearningMetaDat');
cd(currdir)

if (0)
currdir = pwd;
cd('/bluejay5/lucas/analyses/neural/');
save('LearningMetaDat_BACKUP120518', 'LearningMetaDat'); % this is right before combined pu69 - RALMANlearn2 switches into one (i.e. WN off periods combined into the entire learning period)
cd(currdir)
end

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
[a, b] = lt_neural_v2_DIAGN_RemoveOlap(SummaryStruct)

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
% BirdsToPlotList = {'pu69wh78', 'wh44wh39'};
% NeurToPlotList = {[], []}; % leave empty if plot all
% MotifsToPlot = {{'(j)jjbhhg', '(a)abhhg'}, {'(n)hh', '(d)kccbb'}};
BirdsToPlotList = {'gr48bu5'};
NeurToPlotList = {[16 17 18]}; % leave empty if plot all
MotifsToPlot = {{'a(b)'}};

% ======== general params
motifpredur = 0.12;
motifpostdur = 0.075;
preAndPostDurRelSameTimept = 1;
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both
saveON = 0;

Nmax = 44;

for i=1:length(BirdsToPlotList)
    BirdToPlot = BirdsToPlotList{i};
    NeurToPlot = NeurToPlotList{i};
    motiflist = MotifsToPlot{i};
    
    lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
        motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON, ...
        Nmax);
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
BirdToPlot = 'gr48bu5';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = []; % 4 % vector (e.g. [5 7]) - if [] then plots all;
motiflist = {'(a)b'};
% motiflist = {'a(r)dabh', 'r(r)dabh'};
plotbytime = 0; % links rasters for all motifs by time of song.

motifpredur = 0.1;
motifpostdur = 0.15;

plotIndivRaster = 0; % one raster for each neuron/motif
plotCombRast = 1; % one figure, all rasters
plotSmFR = 1; % all smoothed FR.

LearnKeepOnlyBase = 0;

% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase)


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



%% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%% &&&&&&&&&&&&&&&&&&& EXTRACTION OF RAW NEURAL DATA

i=1;
nn=1;
mm=1;

segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
motifpredur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_predur;
PrePostRelSameTime = MOTIFSTATS_Compiled.birds(i).Params_regexp.preAndPostDurRelSameTimept;
motifpostdur = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_postdur;
chanamp = SummaryStruct.birds(i).neurons(nn).channel;

[a, b] = fileparts(SummaryStruct.birds(i).neurons(nn).dirname);

ntrials = length(segextract);
tic
for j=1:ntrials
    
    % ======== figure out filename and time within file
    fname = segextract(j).song_filename;
    motifdur = segextract(j).global_offtime_motifInclFlank - segextract(j).global_ontime_motifInclFlank;
    ontime_token = segextract(j).WithinSong_TokenOns;
    ontime_motif = ontime_token - motifpredur;
    if PrePostRelSameTime==1
        offtime_motif = ontime_token + motifpostdur;
    else
        offtime_motif = ontime_motif + motifdur;
    end
    
    
    % ========= LOAD SONG DATA
    [amplifier_data,~,frequency_parameters, ~, ...
    ~, amplifier_channels, ~, ~] =...
    pj_readIntanNoGui([a '/' fname]);
    
    dat = amplifier_data([amplifier_channels.chip_channel] == chanamp, :);
    
    % ----------------- extract location within dat that you care about.
    t = [1:length(dat)]/frequency_parameters.amplifier_sample_rate;
    indthis = t>=ontime_motif & t<=offtime_motif;
    datthis = dat(indthis);
    
    
    % ##################### GET LFP 
    % ---- FILTER (0.5 to 400hz) (2pole bessel?)
    [datfilt,neuralFiltLow,neuralFiltHi] =lt_neural_filter(datthis, frequency_parameters, 1, ...
    5, 400);
    
    if (0)
    figure; hold on;
    t = [1:length(datthis)]/30000;
    plot(t, datthis, '-k'); plot(t, datfilt, 'r');
    end
    
    % ---- DOWNSAMPLE
    fs_new = 1500;
    factor = floor(frequency_parameters.amplifier_sample_rate/fs_new);
    datfilt = downsample(datfilt, factor);
    datfilt = single(datfilt);
    if (0)
        lt_figure; hold on;
        plot(t, datfilt, 'b');
        plot(t(1:factor:end), datfilt_dn, 'r');
    end
    
    
end
toc


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

lt_neural_POPLEARN_MASTER;


%% ###################################### DIR VS. UNDIR (LEARNING)

lt_neural_DirUndirLearn_MASTER;


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


%% ################################### SONG BOUT ANALYSES
% ============= EXTRACT DATA



% ============= ANALYSES



%% ################################## [TOOLS FOR MOTIFSTATS_Compiled]

%% ==== REMOVE DIR SONG
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);
