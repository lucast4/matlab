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
