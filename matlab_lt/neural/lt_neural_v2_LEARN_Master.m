%% &&&&&&&&&&&&&&&&&&&&& [LEARNING] SINGLE EXPERIMENTS &&&&&&&&&&&&&&&&&&&&&&&&&
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
plottype = 'byneuron'; %
% 'byneuron' - each neuron one fig will all motifs [DEFAULT]
% 'dotprod' - for each bin of trials get dot prod from IN PROGRESS
% 'bysyl' - each plot one syl, all neurons.
DivisorBaseSongs = 1;
lt_neural_v2_ANALY_Learning(SummaryStruct_filt, MOTIFSTATS, plottype, DivisorBaseSongs);




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
% MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
%     collectWNhit);
    Params_regexp.motif_predur = 0.15;
    Params_regexp.motif_postdur = 0.1;
    Params_regexp.preAndPostDurRelSameTimept = 1;
    Params_regexp.RemoveIfTooLongGapDur = [];
    Params_regexp.extractDirSong = 0;
    MOTIFSTATS_Compiled = lt_neural_v2_ANALY_MultExtractMotif(SummaryStruct, ...
    collectWNhit, 0, 1, 0, 1, 1, [], Params_regexp);



%% =========== [PREPROCESSING REQUIRED]
% =============== PICK OUT SAME TYPE/DIFF [LEANRING SPECIFIC]
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



%% ==== FOR EACH LEARNING EXPERIMENT, PLOT TIMELINE OF NEURONS AND LEARNING
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
        if isempty(MOTIFSTATS.params.TargSyls)
            sdafasdf
        end
        
        
        % === PLOT
        lt_figure; hold on;
        newfig=0;
        lt_neural_v2_ANALY_LearningPlot1(SummaryStruct_tmp, MOTIFSTATS, newfig);
        title([birdname '-' exptname]);
%         keyboard
    end
end

%% ======================== PLOT LEARNING (hz), SIMPLE VERSION
bnum = 1;
enum = 4;

summarystruct_tmp = MOTIFSTATS_Compiled.birds(bnum).exptnum(enum).SummaryStruct;
motifstats = MOTIFSTATS_Compiled.birds(bnum).exptnum(enum).MOTIFSTATS;

exptname = MOTIFSTATS_Compiled.birds(bnum).exptnum(enum).exptname;
birdname = MOTIFSTATS_Compiled.birds(bnum).birdname;

lt_neural_v2_ANALY_LearningPlot1(summarystruct_tmp, motifstats);
title([birdname '-' exptname]);

%% ================================== PLOT LEARNING (ALL SYLS)
close all
MeanSubtract =1; % subtract baseline mean?
BirdsToPlot = {'pu69wh78'};
ExptToPlot = {'RALMANOvernightLearn1'}; % expt and bird must intersect.
lt_neural_v2_ANALY_LearnAllSylPlot(MOTIFSTATS_Compiled, ...
    MeanSubtract, BirdsToPlot,ExptToPlot);

%% ============ save learning expt struct

cd /bluejay5/lucas/analyses/neural/LEARNING/

tstamp;

%% ==== FOR EACH LEARNING EXPERIMENT, PLOT SUMMARY STATISTICS
close all;

if (0)
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
end



% ######################################################################
% ####################################### VERSION 2 - SWITCH AS DATAPOINT
% NOTE: neuron inds in switchstruct match those in motifs_compiled

% === PULL OUT RAW FR FOR ALL NEURONS/TRIALS
RemoveTrialsZeroFR = 1;
premotorWind = [-0.05 0]; % [-a b] means "a" sec before onset and "b" sec after offset
% premotorWind = [-0.03 0.025]; % [-a b] means "a" sec before onset and "b" sec after offset
% premotorWind = [-0.05 0]; % [-a b] means "a" sec before onset and "b" sec after offset
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

% ============ DISPLAY LABELS INFORMATION FOR ALL SWITCHES AND NEURONS
% [LEARNING SPECIFIC]
onlyWNonset = 0; % if 0, then all switches, if 1, then only WN on
lt_neural_v2_ANALY_Swtch_DispLabs(MOTIFSTATS_Compiled, SwitchStruct, onlyWNonset);


% #############################################################################
% ==== SEPARATION - compare separation between contexts pre and end of
% learning
% --------- 1) extract dat
onlyUseDatOnSwitchDays=1; % if 1, then restricts analyses to just dat on day of swtich
[DatstructSep, getdprime] = lt_neural_v2_ANALY_Swtch_Separation(MOTIFSTATS_Compiled, ...
    SwitchStruct, onlyUseDatOnSwitchDays);

% --------- 2) plot
close all;
% expttypewanted = 'one targ context';
% expttypewanted = 'mult targ context - samedir';
% expttypewanted = 'mult targ context - diff dir';
expttypewanted=''; % COLLECT ALL
lt_neural_v2_ANALY_Swtch_SepPlot(DatstructSep, MOTIFSTATS_Compiled, ...
    SwitchStruct, getdprime, expttypewanted);



% #############################################################################

% ===========================================
lt_figure; hold on;
lt_plot_text(0, 0.5, 'bird, expt, and neuron all match between SummaryStruct, Motifstats, and Switchstruct', 'r')

% =========== PLOT SUMMARY OF LEARNING [ITERATING SWITCHES]
% PITCH AT END.
close all;
lt_neural_v2_ANALY_LrnSwtchPLOT(MOTIFSTATS_Compiled, SwitchStruct);


% ============ TIMECOURSES FOR NEURAL FOR SWITCHES
close all;
birdname_get = 'gr48bu5'; % keep empty if want all.
exptname_get = 'RALMANLearn3';
switchnum_get = [1];
plotneurzscore=0;
onlyPlotTargNontarg=1;
lt_neural_v2_ANALY_Swtch_Tcourse(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, ...
    onlyPlotTargNontarg)


% ========================= TIMECOURSES, BINNING BY TIME, showing smoothed
% FR and rasters [GOOD]
close all;
birdname_get = 'gr48bu5'; % keep empty if want all.
exptname_get = 'RALMANLearn3';
switchnum_get = [1];
plotneurzscore=0;
FFzscore =1;
onlyPlotTargNontarg=1;
saveFigs =0;
lt_neural_v2_ANALY_Swtch_Tcourse2(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg, saveFigs)


%% ==== BINNED LEARNING
 
close all;
birdname_get = ''; % keep empty if want all.
exptname_get = '';
switchnum_get = [];
Bregion = {'RA'};
plotneurzscore=0;
FFzscore =1;
onlyPlotTargNontarg=3; % 1 is targ/same; 3 is all [DEFAULT: 3]
saveFigs =0;
onlySingleDir =1; % if 1, then only does cases where all targs same dir
removebadsyl = 1;
lt_neural_v2_ANALY_Swtch_Binned(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg, saveFigs, onlySingleDir, Bregion, removebadsyl);


%% ======= BINNED LEARNING V2 (SONG BY SONG)
close all;
birdname_get = 'wh72pk12'; % JUST FOR PLOTTING
exptname_get = 'RALMANLearn3'; 
switchnum_get = [9];
Bregion = {'RA'};
plotneurzscore=0;
FFzscore =1;
onlyPlotTargNontarg=1;
saveFigs =0;
onlySingleDir =1; % if 1, then only does cases where all targs same dir
lt_neural_v2_ANALY_Swtch_Binned2(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
    onlyPlotTargNontarg, saveFigs, onlySingleDir, Bregion);


%% PLOT LEANRING (FF) TRAJECTORIES [and rasters...]

numbirds = length(SwitchStruct.bird);
for i=1:numbirds
    birdname = SwitchStruct.bird(i).birdname;
    numexpts = length(SwitchStruct.bird(i).exptnum);
    
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
        
        for iii=1:numswitch
            
            close all;
            birdname_get = birdname; % keep empty if want all.
            exptname_get = exptname;
            switchnum_get = [iii];
            plotneurzscore=0;
            FFzscore =1;
            onlyPlotTargNontarg=1;
            saveFigs =1;
            lt_neural_v2_ANALY_Swtch_Tcourse2(MOTIFSTATS_Compiled, SwitchStruct, ...
                birdname_get, exptname_get, switchnum_get, plotneurzscore, FFzscore, ...
                onlyPlotTargNontarg, saveFigs)
            
            
        end
        
        
    end
end



%%
% ========================== SPIKE COUNT CORRELATIONS, CHANGE DURING
% LEARNIG? Looks at a single switch.
% NOTE: all neurons must have same batch file for this to work.
lt_neural_v2_ANALY_Swtch_LearnNeurCorr(MOTIFSTATS_Compiled, SwitchStruct, ...
    birdname_get, exptname_get, switchnum_get, plotneurzscore)


% ============== SUMMARIZE FOR EACH EXPT (I.E. NEURAL AND FF CHANGE FOR
% TARG AND OTHER SYLS)
close all;
RemoveLowNumtrials = 1; % min number for both base and train(double)
MinTrials = 5; % for removing
skipMultiDir = 1;
usePeakLearn = 0; % assumes that care about learning for targ 1. (median time of max learning)
useTrialsRightAfterWNOn = 0; % note, if 1, then overrides usePeakLearn; if n is numtrials, takes the trials from n+1:2n from WN onset.


% ---- what metric to plot for main figs
% 1) split into half, corr those, repeat
fieldname_baseneur = 'AllNeurSplitCorrBase';
fieldname_trainneur = 'AllNeurSplitCorrTrain';
% fieldname_trainneur = 'AllNeurSplitCorrTrainvsTrain';
% 2) permute train and base, get corr, repeat - use that as base...
% fieldname_baseneur = 'AllNeurCorrShuffMean';
% fieldname_trainneur = 'AllNeurCorrDat';
% 3) old version, no permuting, each trial corr to base "template"
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
OnlyKeepSigLearn = 0; % either regression or end training has be significant.
learnsigalpha = 0.05; % for deciding that an experiment showed "no learning"

% ---- ONLY KEEP SWITCHES STARTING FROM WN OFF
OnlyKeepWNonset =0; % if 1, then yes, if 2, then only keeps those with WN transition (not onset); if 0, then takes all

% --- only use data on day of switch - i.e. don't go to next day for
% training end
OnlyUseDatOnSwitchDay=1; % NOTE: if use with "UsePeakLearn" then will constrain to be within switch day
% i.e. if peak learn is after switch day then will take end of first day...

[DATSylMot, ~] = lt_neural_v2_ANALY_Swtch_Summary(MOTIFSTATS_Compiled, SwitchStruct, RemoveLowNumtrials, ...
    MinTrials, UseZscoreNeural, neuralmetricname, fieldname_baseneur, fieldname_trainneur, ...
    skipMultiDir, usePeakLearn, plotLearnStatsOn, learnsigalpha, OnlyKeepSigLearn, ...
    OnlyKeepWNonset, OnlyUseDatOnSwitchDay, useTrialsRightAfterWNOn);




% ============================================ PAIRWISE CORRELATIONS DUR
% LEARNING 
% [IN PROGRESS!!! - USE DIFFERENT CODE, POPULATION STRUCTURE]
lt_neural_v2_LEARNING_NeurPairCorr(MOTIFSTATS_Compiled, SwitchStruct);



% ================ IN LINEAR MODEL IS THERE EFFECT OF SYLLABLE TYPE AFTER
% CONTROLLING FOR VARIOUS THINGS?
lt_neural_v2_ANALY_Swtch_LME(DATSylMot)



%% ########################################################################
%% ============= LEARNING CHANGE, SEPARATE BY BASELINE FR-PITCH CORRELATION
%% ========== first run preprocessgin stuff above
% ============== 1) LOAD MOTIFSTATS
load(['/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_06Jul2018_1227.mat']);
load(['/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/Params_06Jul2018_1227.mat']);

load(['/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_07Mar2019_0024.mat']);
load(['/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/Params_07Mar2019_0024.mat']);

load(['/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_17Mar2019_1032.mat']);
load(['/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/Params_17Mar2019_1032.mat']);


lt_neural_v2_ANALY_FRsmooth_Pre; % this just references functions from above.


%% ================== DIAGNOSTIC [TO CHECK RAW NEURAL FOR EACH SYLLABLWE]
% ========= WILL PLOT. WILL OVERRIDE PREVIOUS PLOTS.
% === for each experiment list all neruons, directories, song files...
% PLOTS IN DIRWCTORY, INSTANCES OF THE SYLALBLE...

% ======= IMPORTANT: MUST FIRST EXTRACT OUTDAT - WILL ONLY SAVE THINGS THAT
% ARE EXTRACTED THERE (I.E. IF ONLY GOOD SYLS THERE, THEN HERE ALSO...);
if (0)
onlyfirstswitch = 1;

for i=1:length(SwitchStruct.bird)
    
    for ii=1:length(SwitchStruct.bird(i).exptnum)
        for ss=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist)
            
            if onlyfirstswitch==1
                if ss>1
                    continue
                end
            end
            
            neurlist = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).goodneurons;
            if isempty(neurlist)
                continue
            end
            if isempty(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(neurlist(1)).DATA)
                continue
            end
            if ~any([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(neurlist(1)).DATA.motif.trainInds])
                continue
            end
            
           
            bname = SwitchStruct.bird(i).birdname;
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            
            motiflist = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).STATS_motifsyl.sylname};
            
            time1 = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            time2 = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            tdatenum_window = [time1 time2];
            
            % ======================== PLOT AND SAVE RAW
            close all;
            BirdToPlot = bname;
            % % ---- give it either
            % A) one neuron and a bunch of motifs or
            % B) bunch of neurons and one motif
            NeurToPlot = neurlist'; % 4 % vector (e.g. [5 7]) - if [] then plots all;
            % motiflist = {'a(b)', 'jbh(h)g'}; 
            % motiflist = {'(d)kcc', 'dk(c)c', '(n)hh', 'c(b)'};
%             motiflist = ngrams';
            
            % motifpredur = 0.15;
            % motifpostdur = 0.15;
            motifpredur = 0.125;
            motifpostdur = 0.075;
            preAndPostDurRelSameTimept = 1;
            
            % --- 1) directed song
            PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both
            
            saveON = 1;
            Nmax = 44;
            savedirmain = ['/bluejay0/bluejay2/lucas/analyses/neural/LEARN/RA/FIGS/DIAGN_PlotRawNeural/' bname '-' exptname];
            
            % ================= FIND TIME WINDOW OF DATA TO PLOT
            indsthis = OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                OUTDAT.All_swnum==ss;
            % ======= BASELINE
            twind_inds = cellfun(@(x)[x(1) x(end)], OUTDAT.AllBase_indsepoch(indsthis), ...
                'UniformOutput', 0);
            tvals = OUTDAT.All_FF_t(indsthis, 1);
            twind_all =[];
            for k=1:length(twind_inds)
                twind_all = [twind_all; tvals{k}(twind_inds{k})];
            end
            tdatenum_window = [min(twind_all(:,1)) max(twind_all(:,2))];

            saveON = 1;
            Nmax = 22;
            savedirmain = ['/bluejay0/bluejay2/lucas/analyses/neural/LEARN/RA/FIGS/DIAGN_PlotRawNeural/' bname '-' exptname '-BASE'];
            lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
                motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON, ...
                Nmax, savedirmain, tdatenum_window);
            
            % ================= WN
            twind_inds = cellfun(@(x)[x(epochtoplot) x(epochtoplot+1)-1], OUTDAT.AllWN_indsepoch(indsthis), ...
                'UniformOutput', 0);
            tvals = OUTDAT.All_FF_t(indsthis, 2);
            twind_all =[];
            for k=1:length(twind_inds)
                twind_all = [twind_all; tvals{k}(twind_inds{k})];
            end
            tdatenum_window = [min(twind_all(:,1)) max(twind_all(:,2))];
            
            saveON = 1;
            Nmax = 22;
            savedirmain = ['/bluejay0/bluejay2/lucas/analyses/neural/LEARN/RA/FIGS/DIAGN_PlotRawNeural/' bname '-' exptname '-WN'];
            lt_neural_DIAGN_PlotRawNeural(SummaryStruct, BirdToPlot, NeurToPlot, motiflist, ...
                motifpredur, motifpostdur, PlotDirSong, preAndPostDurRelSameTimept, saveON, ...
                Nmax, savedirmain, tdatenum_window);
            
        end
    end
end

end



%% ================== DISPLAY OVERVIEW OF DATA THAT IS EXTRACTED
onlyfirstswitch = 1;
fid = fopen('/bluejay0/bluejay2/lucas/analyses/neural/LEARN/RA/exptlog.txt', 'w');

for i=1:length(SwitchStruct.bird)
    
    for ii=1:length(SwitchStruct.bird(i).exptnum)
        for ss=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist)
            
            if onlyfirstswitch==1
                if ss>1
                    continue
                end
            end
            
            neurlist = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).goodneurons;
            
            if isempty(neurlist)
                continue
            end
            if isempty(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(neurlist(1)).DATA)
                continue
            end
            if ~any([SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(neurlist(1)).DATA.motif.trainInds])
                continue
            end
            
            nuerlist2 = find(MOTIFSTATS_Compiled.birds(i).exptnum(ii).neurIDOriginal_inorder);
            chansthis = [SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron.chan];
            assert(all([SummaryStruct.birds(i).neurons(nuerlist2).channel]==chansthis));

            bname = SwitchStruct.bird(i).birdname;
            exptname = SwitchStruct.bird(i).exptnum(ii).exptname;
            
            motiflist = {SwitchStruct.bird(i).exptnum(ii).switchlist(ss).STATS_motifsyl.sylname};
            
            time1 = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_previous;
            time2 = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).switchdnum_next;
            tdatenum_window = [time1 time2];
            
            % ================= FIND TIME WINDOW OF DATA TO PLOT
            indsthis = OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & ...
                OUTDAT.All_swnum==ss;

            % ======= BASELINE
            twind_inds = cellfun(@(x)[x(1) x(end)], OUTDAT.AllBase_indsepoch(indsthis), ...
                'UniformOutput', 0);
            tvals = OUTDAT.All_FF_t(indsthis, 1);
            twind_all =[];
            for k=1:length(twind_inds)
                twind_all = [twind_all; tvals{k}(twind_inds{k})];
            end
            tdatenum_window_base = [min(twind_all(:,1)) max(twind_all(:,2))];

            % ================= WN
            twind_inds = cellfun(@(x)[x(epochtoplot) x(epochtoplot+1)-1], OUTDAT.AllWN_indsepoch(indsthis), ...
                'UniformOutput', 0);
            tvals = OUTDAT.All_FF_t(indsthis, 2);
            twind_all =[];
            for k=1:length(twind_inds)
                twind_all = [twind_all; tvals{k}(twind_inds{k})];
            end
            tdatenum_window_wn = [min(twind_all(:,1)) max(twind_all(:,2))];
            
            
            % ============ DISPLAY DATA
            disp('==========================================');
            disp([bname '-' exptname '-sw' num2str(ss)]);
            disp(['NEURONS: ' num2str(nuerlist2)]);
            disp(['CHANS: ' num2str([SummaryStruct.birds(i).neurons(nuerlist2).channel])]);
            disp(['BATCH: ' SummaryStruct.birds(i).neurons(nuerlist2).batchfilename]);
            disp(['MOTIFS: ' motiflist]);
            disp(['TIME (bounds): ' datestr([tdatenum_window_base(1) ], 'ddmmmyyyy-HHMM')]);
            disp(['TIME (bounds): ' datestr([tdatenum_window_wn(end)], 'ddmmmyyyy-HHMM')]);
   
            % ============ write to text file
            fprintf(fid,'%s\n', '==========================================');
            fprintf(fid,'%s\n',[bname '-' exptname '-sw' num2str(ss)]);
            fprintf(fid,'%s\n',['NEURONS: ' num2str(nuerlist2)]);
            fprintf(fid,'%s\n',['CHANS: ' num2str([SummaryStruct.birds(i).neurons(nuerlist2).channel])]);
            fprintf(fid,'%s\n',['BATCH: ' SummaryStruct.birds(i).neurons(nuerlist2).batchfilename]);
            fprintf(fid,'%s\n',['TIME (bounds): ' datestr([tdatenum_window_base(1) ], 'ddmmmyyyy-HHMM')]);
            fprintf(fid,'%s\n',['TIME (bounds): ' datestr([tdatenum_window_wn(end)], 'ddmmmyyyy-HHMM')]);
            cellfun(@(x)fprintf(fid, '%s', [', ' x]), motiflist);
            fprintf(fid,'%s\n','');

        end
    end
end
fclose(fid);

%% ================== DIAGNOSTIC
if (0)
    sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, ss, motifthis);
end

%%
clear OUTDAT
% *************************************** EXTRACTION
onlyFirstSwitch = 1; % then only first switch...
onlyIfSameTarg = 0; % only if targ are all same sylalbles
onlyiftargsamedir = 1; % only if targs are all driving same directyion.
% BirdsToPlot = {'pu69wh78', 'wh44wh39'};
BrainLocation = {'RA'};
BirdsToPlot = {};
% BrainLocation = {};
throwoutlonggap = 0; % gap betwen end of base and start of train.
% ---- FILTERING DATA
removebadsyls =1;
removebadchan = 1;
removebadtrials = 1;
OUTDAT = lt_neural_v2_ANALY_FRsmooth(MOTIFSTATS_Compiled, SwitchStruct, onlyFirstSwitch, ...
    onlyIfSameTarg, BirdsToPlot, BrainLocation, throwoutlonggap, onlyiftargsamedir, ...
    removebadsyls, removebadchan, removebadtrials);


% =========== 1) PERFORM BASELINE SUBTRACT
usepercent = 3; % 0 = mean; 1 = pecent; 2 = zscore; 3 = dprime
% nbasetime = []; %
% nbasetime = 90; % minuts of baseline to take, from first train.
% nbasetime_ignoreswitch1 = 1; % defautl 1, useful to restrict baseline for switches only.
nbasetime = 60; % minuts of baseline to take, from first train.
nbasetime_ignoreswitch1 = 0; % defautl 1, useful to restrict baseline for switches only.
prctile_divs = [33 66 100]; % percentiles to divide up data by
% prctile_divs = [25 50 75 100]; % percentiles to divide up data by
% prctile_divs = [50 100]; % percentiles to divide up data by
OUTDAT = lt_neural_v2_ANALY_FRsmooth_MinusBase(OUTDAT, SwitchStruct,...
    usepercent, nbasetime, nbasetime_ignoreswitch1, prctile_divs);


% ============ [OPTIONAL] SCALE FIRING RATES TO MATCH WN TO BASELINE FOR
% EACH EXPT
% close all;
prewind = [-0.09 0.02]; % window to use to match FR 
ignorecontext = 1;
epochtoplot = 3;
meanver = 0; % 0: first take avearge, then take ratio; 1: first take ratios, then take average.
OUTDAT = lt_neural_v2_ANALY_FRsmooth_NormFR(OUTDAT, prewind, ignorecontext, ...
    epochtoplot, meanver, SwitchStruct, usepercent, nbasetime, ...
    nbasetime_ignoreswitch1, prctile_divs);


% =========== 2B) VARIOUS MEASURES COMPARED TO BASELINE
shuffSylType = 0;
epochtoplot = 3;
plotOn = 0;
OUTDAT = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, shuffSylType, ...
    epochtoplot, plotOn);


% ============= EXTRACT MAGNITUDE OF LEARNING
% --- in direction of training
OUTDAT.AllMinusBase_PitchZ_Ldir = ...
    OUTDAT.AllMinusBase_PitchZ(:, epochtoplot).*OUTDAT.AllTargLearnDir;


% ========== simplify future analyses, extract minus base, within epoch.
% tmp = cellfun(@(x)(x{3}), OUTDAT.AllMinusBase_FRmeanAll, 'UniformOutput', 0);
% tmp = cellfun(@(x)x', tmp, 'UniformOutput', 0);
% OUTDAT.AllMinusBase_FRmeanAll_epoch = tmp;


% *************************************** PLOTS
% ####################### PLOTTING RAW DATA
% ============== PLOT ALL NEURONS/MOTIFS [BASE, WN, OVERLAY]
close all;
% epochtoplot
dopause = 0;
lt_neural_v2_ANALY_FRsmooth_PlotSimple(OUTDAT, MOTIFSTATS_Compiled, SwitchStruct, ...
    epochtoplot, dopause)

% ============== PLOT ALL NEURONS/MOTIFS [FIRING RATE]
close all;
birdtoplot = []; % empty to p[lot all (will do pause)
exptoplot = [];
onlyplotTargs = 1;
lt_neural_v2_ANALY_FRsmooth_Plot(OUTDAT, MOTIFSTATS_Compiled, SwitchStruct, ...
    birdtoplot, exptoplot, onlyplotTargs)

% =========== 2) PLOT EACH DATAPOINT (i.e. raw dat)
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
close all
plotDevFromBase =0;
plotNormMeasures = 1;
lt_neural_v2_ANALY_FRsmooth_BasePlots(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, epochtoplot, plotDevFromBase, plotNormMeasures)

% ========== 4) PLOT RAW DAT, INCLUDING DERIVED MEASURES
close all;
% corrwindow = [-0.04 0.06];
% corrwindow = [-0.11 0.01];
corrwindow = [-0.08 0.03];
dontclosefig=1;
ignoreDiff = 0;
lt_neural_v2_ANALY_FRsmooth_Plot2(OUTDAT, MOTIFSTATS_Compiled, SwitchStruct, ...
    corrwindow, dontclosefig, ignoreDiff);


%%
% =========== SUMMARY OF DATA THAT IS COLELCTED (I.E. TIMING, NUM SWITCHES,
% ETC)



% ########################## SUMMARY ANALYSES
% ========== 3) ABSOLUTE VALUE OF LEARNING GREATER AT TARGET? 
% NOTE, CAN COMPARE TO SHUFFLE AS WELL
close all;
% -- what epoch to plot (dividing up learning into percentiles)
% analytype = 'AllOnlyMinusDiff_FRsmooth';
% analytype = 'AllDevAll_NotAbs';
% analytype = 'AllOnlyMinusBase_FRsmooth'; % just minus base, not abs.
% analytype = 'AllOnlyMinusBase_FRsmooth_abs'; % just minus base, not abs.
analytype = 'AllMinusAll_FRsmooth';
doShuff=1;
syltypesneeded = [1 0 1];
% premotorwind = [-0.08 0.03]; % ORIGINAL
% premotorwind = [-0.1 0.02];
premotorwind = [-0.07 0.015];
nshuffs = 200;
% minmotifs = 5; % min motif types to include.
minmotifs = 4; % min motif types to include.
lt_neural_v2_ANALY_FRsmooth_BaseMinu(OUTDAT, SwitchStruct, ...
    epochtoplot, analytype, doShuff, syltypesneeded, premotorwind, nshuffs, ...
    minmotifs);


% ########################## SUMMARY- DIRECTLY COMPARE TARG TO DIFF TYPE
close all;
% -- what epoch to plot (dividing up learning into percentiles)
% analytype = 'AllOnlyMinusDiff_FRsmooth';
% analytype = 'AllDevAll_NotAbs';
analytype = 'AllOnlyMinusBase_FRsmooth'; % just minus base, not abs.
% analytype = 'AllOnlyMinusBase_FRsmooth_abs'; % just minus base, not abs.
% analytype = 'AllMinusAll_FRsmooth';
syltypesneeded = [1 0 1];
% premotorwind = [-0.08 0.03]; % ORIGINAL
% premotorwind = [-0.1 0.02];
premotorwind = [-0.07 0.015];
minmotifs = 4; % min motif types to include.
datlevel = 'unit';
% datlevel = 'switch';
lt_neural_v2_ANALY_FRsmooth_Compare(OUTDAT, SwitchStruct, ...
    analytype, syltypesneeded, premotorwind, minmotifs, datlevel);

%% ============ PLTO AND SAVE RASTERS FOR ALL NEUROSN/MOTIFS
close all;
Ntrials = 40; % rasters, each base and wn

for i=1:length(SwitchStruct.bird)
    bname = SwitchStruct.bird(i).birdname;
    for ii=1:length(SwitchStruct.bird(i).exptnum)
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
        for ss=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist)
            
            for nn=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron)
                
                figcount=1;
                subplotrows=2;
                subplotcols=2;
                fignums_alreadyused=[];
                hfigs=[];
                hsplots = [];
                
                if ~isfield(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn), 'DATA')
                    continue
                end
                if isempty(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA)
                    continue
                end
                if isempty(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA)
                    continue
                end
                
                for mm=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif)
                    
                    bnum = i;
                    enum = ii;
                    
                    indsOUT = find(OUTDAT.All_birdnum==bnum & OUTDAT.All_exptnum==enum & OUTDAT.All_swnum==ss ...
                        & OUTDAT.All_motifnum==mm & OUTDAT.All_neurnum==nn);
                    
                    if isempty(indsOUT)
                        continue
                    end
                    
                    assert(length(indsOUT)==1);
                    disp([num2str(bnum) '-' num2str(enum) '-sw' num2str(ss) '-n' num2str(nn) '-m' num2str(mm)])
                    
                    motifname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.motif_regexpr_str{mm};
                    motif_predur = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.motif_predur;
                    
                    % ============ PLOT
                    indsbase = OUTDAT.All_indsbase_all{indsOUT};
                    indswn = OUTDAT.All_indstrain_all{indsOUT};
                    segextract = MOTIFSTATS_Compiled.birds(bnum).exptnum(enum).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
                    
                    if length(indsbase)>Ntrials
                        indsbase = indsbase(sort(randperm(length(indsbase), Ntrials)));
                    end
                    if length(indswn)>Ntrials
                        indswn = indswn(sort(randperm(length(indswn), Ntrials)));
                    end
                    
                    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                    title([bname '-' ename '-sw' num2str(ss) '-n' num2str(nn)]);
                    ylabel(motifname);
                    trialnum = [indsbase indswn];
                    spks = {segextract(trialnum).spk_Times};
                    for j=1:length(trialnum)
                        sp = spks{j};
%                         tr = trialnum(j);
                        %                         if any(ismember(indsbase, tr))
                        %                             pcol = 'k';
                        %                         elseif any(ismember(indswn, tr))
                        %                             pcol = 'r';
                        %                         end
                        
%                         lt_neural_PLOT_rasterline(sp, tr, 'k');
                        lt_neural_PLOT_rasterline(sp, j, 'k');
                    end
                    axis tight;
%                     line(xlim, [max(indsbase)+0.5 max(indsbase)+0.5]);
                    line(xlim, [length(indsbase)+0.5 length(indsbase)+0.5]);
                    line([motif_predur motif_predur], ylim);
                end
                
                % =============== SAVE ALL FIGS
                sdir = ['/bluejay0/bluejay2/lucas/analyses/neural/LEARN/RA/FIGS/PlotRast/' bname '-' ename '-sw' num2str(ss) '-n' num2str(nn)]; 
                mkdir(sdir);
                cd(sdir);
                try
                lt_save_all_figs;
                catch err
                end
                close all;
            end
        end
    end
end

%% ==== PLOT RASTERS FOR A SINGLE NEURON
Ntrials = 30; % rasters, each base and wn
i=4;
ii=3;
ss=1;
nn=1;


% ======================
figcount=1;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for mm=1:length(SwitchStruct.bird(i).exptnum(ii).switchlist(ss).neuron(nn).DATA.motif)
    
    bnum = i;
    enum = ii;
        bname = SwitchStruct.bird(i).birdname;
        ename = SwitchStruct.bird(i).exptnum(ii).exptname;
    indsOUT = find(OUTDAT.All_birdnum==bnum & OUTDAT.All_exptnum==enum & OUTDAT.All_swnum==ss ...
        & OUTDAT.All_motifnum==mm & OUTDAT.All_neurnum==nn);
    
    if isempty(indsOUT)
        continue
    end
    
    assert(length(indsOUT)==1);
    disp([num2str(bnum) '-' num2str(enum) '-sw' num2str(ss) '-n' num2str(nn) '-m' num2str(mm)])
    
    motifname = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.motif_regexpr_str{mm};
    motif_predur = MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params.motif_predur;
    
    % ============ PLOT
    indsbase = OUTDAT.All_indsbase_all{indsOUT};
    indswn = OUTDAT.All_indstrain_all{indsOUT};
    segextract = MOTIFSTATS_Compiled.birds(bnum).exptnum(enum).MOTIFSTATS.neurons(nn).motif(mm).SegmentsExtract;
    
    if length(indsbase)>Ntrials
        indsbase = indsbase(sort(randperm(length(indsbase), Ntrials)));
    end
    if length(indswn)>Ntrials
        indswn = indswn(sort(randperm(length(indswn), Ntrials)));
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(ss) '-n' num2str(nn)]);
    ylabel(motifname);
    trialnum = [indsbase indswn];
    spks = {segextract(trialnum).spk_Times};
    for j=1:length(trialnum)
        sp = spks{j};
        %                         tr = trialnum(j);
        %                         if any(ismember(indsbase, tr))
        %                             pcol = 'k';
        %                         elseif any(ismember(indswn, tr))
        %                             pcol = 'r';
        %                         end
        
        %                         lt_neural_PLOT_rasterline(sp, tr, 'k');
        lt_neural_PLOT_rasterline(sp, j, 'k');
    end
    axis tight;
    %                     line(xlim, [max(indsbase)+0.5 max(indsbase)+0.5]);
    line(xlim, [length(indsbase)+0.5 length(indsbase)+0.5]);
    line([motif_predur motif_predur], ylim);
end


        
%% ============ HEAT MAP PLOT. FOR EACH NEURON, PLOT ALL SYLALBLES FR MODULATION IN HEAT MAP

close all;
% -- what epoch to plot (dividing up learning into percentiles)
% analytype = 'AllOnlyMinusBase_FRsmooth_abs';
analytype = 'AllOnlyMinusBase_FRsmooth';
syltypesneeded = [1 0 1];
minmotifs = 4; % min motif types to include.
prewind = [-0.07 0.015];
lt_neural_v2_ANALY_FRsmooth_FRheatmap(OUTDAT, SwitchStruct, ...
    epochtoplot, analytype, syltypesneeded, minmotifs, prewind)

%%

% ========== 4) PREDICT FR CHANGE BASED ON LEARNING?
close all;
corrwindow = [-0.07 0.015];
% corrwindow = [-0.1 0.02];
syltypesneeded = [1 0 1];
% analytoplot = 'AllDevAll_NotAbs';
analytoplot = 'AllOnlyMinusBase_FRsmooth';
syltoplot = 'targ';
% syltypesneeded = [1 1 1];
docorrvsdiff = 1; % default is 1 
onlyPlotIfAllTargSameDir=1;
onlyIfLearnCorrectDir = 0; % onoy applies to summary plots
doregression = 0; % mixed effects model.
plotRaw = 1; % i.e. each experiment broken out.
plotSummary = 1;
OutStruct = lt_neural_v2_ANALY_FRsmooth_PredLearn(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, corrwindow, syltypesneeded, epochtoplot, analytoplot, ...
    syltoplot, docorrvsdiff, onlyPlotIfAllTargSameDir, onlyIfLearnCorrectDir, ...
    doregression, plotRaw, plotSummary);


% ------------  TO ITERATE OVER MULTIPLE TIME WINDOWS 
% [USES ABOVE CODE]
close all;
timeshift = 0.01;
windsize = 0.08;
syltoplot = 'targ';
syltypesneeded = [1 0 1];
analytoplot = 'AllDevDiff_NotAbs';
onlyIfLearnCorrectDir = 1;
lt_neural_v2_ANALY_FRsmooth_PredLMult(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, syltoplot, timeshift, windsize, syltypesneeded, ...
    epochtoplot, analytoplot, onlyIfLearnCorrectDir);

% --------- TO ITERATE OVER MULTIPLE TIME WINDOWS AND EPOCHS
% [USES ABOVE CODE]
lt_neural_v2_ANALY_FRsmooth_PredLMult2



% =============== [v1] CHANGE IN FR SIMILAR FOR SAME TYPE?
timewind = [-0.08 0.02];
onlyifonetarget = 1; % haven't coded up for two targets yet..
usediffFromBase = 1; % if 0, then also norms to global drift.
syltypesneeded = [1 1 1];
lt_neural_v2_ANALY_FRsmooth_CompSylTypes(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, timewind, onlyifonetarget, usediffFromBase, syltypesneeded);



% ###################### [v2] COMPARE SYLTYPES
close all;
nshuff = 500;
usediffFromBase = 1; % if 0, then also norms to global drift.
plotdifftype = 1;
[rhomean_dat, rhomean_shuff, rhomean_dat_diff, rhomean_shuff_diff] ...
    = lt_neural_v2_ANALY_FRsmooth_CompSylTypes3(OUTDAT,...
    SwitchStruct, MOTIFSTATS_Compiled, nshuff, usediffFromBase, plotdifftype, ...
    epochtoplot);


lt_figure; hold on;

% ==== 1) 
lt_subplot(2,2,1); hold on;
title('shuff and dat');
xlabel('targ-same corr');
lt_plot_histogram(rhomean_shuff)
line([rhomean_dat rhomean_dat], ylim, 'Color', 'r');
p = (sum(rhomean_shuff>=rhomean_dat)+1)/(nshuff+1);
lt_plot_pvalue(p, 'vs shuff', 1);

% ==== 2)
lt_subplot(2,2,2); hold on;
title('shuff and dat');
xlabel('[targ-smae corr] - [targ-diff corr]');
diffshuff = rhomean_shuff - rhomean_shuff_diff;
diffdat = rhomean_dat - rhomean_dat_diff;
lt_plot_histogram(diffshuff);
line([diffdat diffdat],  ylim, 'Color', 'r');
p = (sum(diffshuff>=diffdat)+1)./(nshuff+1);
lt_plot_pvalue(p, 'vs shuff', 1);



% ***************************************************************
% TROUBLESHOOTING - to plot time points for all trials.
if (0)
    figure; hold on;
    
    tbase = OUTDAT.All_FF_t{j,1};
    twn = OUTDAT.All_FF_t{j,2};
    
    plot(tbase, 0, 'ok');
    plot(twn, 0, 'or');
end


