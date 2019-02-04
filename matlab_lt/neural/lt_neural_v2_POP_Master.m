%% ###############################################################
%% ##################################################### POPULATION
%% ======================== EXTRACT SEGMENTS FOR POPULATIONS
% ========== KEY DIFFERENCES,
% 1) MIGHT HAVE TO SKIP SONGS TO MAKE SURE ALL
% NEURONS ARE REPRESENTED IN ALL SONGS
% 2)


close all; clear MOTIFSTATS_Compiled;
lt_neural_ExtractMotifs_Regular;


%% ================ REMOVE DIR SONG 

MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);


%% ================ extraction continued
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;

%% ================ 
birdnum = 4;

lt_neural_POP_PlotRast


%% %%%%%%%%%%%%%%% VERSION 1
%% ================ PLOT [CORRELATION WITH FF]
close all;
MOTIFSTATS_pop = lt_neural_POP_FFcorr(MOTIFSTATS_pop, SummaryStruct);

close all;
lt_neural_POP_FFcorrPlot

%% %%%%%%%%%%%%%%% VERSION 2
%% ================ PLOT [CORRELATION WITH FF]
close all;
xcov_dattotake = [-0.08 0.030];
xcov_dattotake = [-0.1 0.04];
% xcovwindmax = 0.04;
% binsize_spk = 0.005;
xcovwindmax = 0.05;
binsize_spk = 0.0025;
MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk);

% 
% %% =============== PLOT, DISTRIBUTIONS ACROSS ALL MOTIFS/BIRDS
% numbirds = length(MOTIFSTATS_pop.birds);
% 
% for i=1:numbirds
%  numexptMOTIFSTATS_pop.birds(i)
% 
% end

%% =============== SUMMARY PLOT OF ALL CROSS-CORRELATIONS [GOOD - XCOV OF SPIKING]
% i.e. spike-spike xcov
plotRaw = 1;
birdstoplot = [];
% OUTSTRUCT = lt_neural_POP_PlotSummary(MOTIFSTATS_pop, SummaryStruct, plotRaw, ...
%     birdstoplot);
OUTSTRUCT = lt_neural_POP_PlotSummary_v2(MOTIFSTATS_pop, SummaryStruct, plotRaw, ...
    birdstoplot);

close all;
plotRaw =0;
% lt_neural_POP_PlotSummary2(MOTIFSTATS_pop, SummaryStruct, OUTSTRUCT, ...
%     plotRaw);
lt_neural_POP_PlotSummary2_v2(MOTIFSTATS_pop, SummaryStruct, OUTSTRUCT, ...
    plotRaw);

%% ================= PLOT DISTRIBUTIONS ACROSS ALL BIRDS



%% #############################################################
%% #############################################################
%% ######################## POPULATION - TAKE ENTIRE MOTIF

% ================ extract full motifs and time warp here
lt_neural_v2_ExtractFullMotifs;


%% ======================= REMOVE DIR SONG

MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);


%% summary struct
SummaryStruct = MOTIFSTATS_Compiled.SummaryStruct;

%% ======================== EXTRACT POPULATION
close all;
MOTIFSTATS_pop = lt_neural_v2_POP_ExtractMotifs(MOTIFSTATS_Compiled, SummaryStruct);
% clear MOTIFSTATS_Compiled;


%% ================ PLOT [CORRELATION WITH FF] [RUN THIS FIRST!!]
close all;
xcov_dattotake = [-0.02 0.02];
xcov_dattotake_anchorpoints = 2; % 
xcovwindmax = 0.15;
binsize_spk = 0.0025;
MOTIFSTATS_pop = lt_neural_POP_ExtractXCov(MOTIFSTATS_pop, SummaryStruct, ...
    xcov_dattotake, xcovwindmax, binsize_spk, xcov_dattotake_anchorpoints);


%% =============== SUMMARY PLOT OF ALL CROSS-CORRELATIONS
close all;
% i.e. spike-spike xcov
plotRaw = 1;
birdstoplot = [];
OUTSTRUCT = lt_neural_POP_PlotSummary(MOTIFSTATS_pop, SummaryStruct, plotRaw, ...
    birdstoplot);
% OUTSTRUCT = lt_neural_POP_PlotSummary_v2(MOTIFSTATS_pop, SummaryStruct, [], ...
%     birdstoplot);

close all;
plotRaw =0;
lt_neural_POP_PlotSummary2(MOTIFSTATS_pop, SummaryStruct, OUTSTRUCT, ...
    plotRaw);
% lt_neural_POP_PlotSummary2_v2(MOTIFSTATS_pop, SummaryStruct, OUTSTRUCT, ...
%     plotRaw);

