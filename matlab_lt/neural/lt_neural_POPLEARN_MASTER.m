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


%% ==== MAKE SURE SUMMARY STRCUT IS THE CORRECT ONE
SummaryStruct = MOTIFSTATS_Compiled.SummaryStruct;

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
BirdExptPairsToPlot = {'pu69wh78', 'RALMANOvernightLearn1'};
motiftoplot = 'aa(b)';

lt_neural_POPLEARN_PlotLearnTraj(MOTIFSTATS_pop ,SwitchStruct, ...
    SummaryStruct, BirdExptPairsToPlot, motiftoplot);


%% #######################################################################
%% ############################# CROSS CORR CHANGE DURING LAERNING
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


%% ============= VARIOUS PLOTS OF COHERENCE
i=1;
ii=1;

numsets = length(MOTIFSTATS_pop.birds(i).exptnum(ii).Sets_neurons);

for ss = 1:numsets
    
    dat = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss);
    
    nummotifs = length(dat.motif);
    
    for mm=1:nummotifs
       
        % ======= GET LIST OF CHANNELS AND ASSOCIATED BRAIN REGIONS
        numneur = length(dat.motif(mm).SegExtr_neurfakeID);
        Chanlist = [];
        Bregionlist = {};
        for n=1:numneur
            nID = dat.motif(mm).SegExtr_neurfakeID(n).neurID_orig;
            chan = SummaryStruct.birds(i).neurons(nID).channel;
            bregion = SummaryStruct.birds(i).neurons(nID).NOTE_Location;
            
            Chanlist = [Chanlist; chan];
            Bregionlist = [Bregionlist; bregion];
        end
        
        % ====== FOR EACH PAIR OF CHAN COLLECT ALL TRIALS IN REGEXP
        
        segextract = dat.motif(mm).SegExtr_neurfakeID(1).SegmentsExtract;
        chan1 = Chanlist(1);
        chan2 = Chanlist(2);
        
        
        
    end
end
    
    
    
    
    
    
    



