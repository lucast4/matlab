%% ===== EXAMPLES FOR PAPER

% ###################### pu69wh78
% pu69 - 17(LMAN) NO, not good, shows some drift (naa, gaa)
% pu69 - 1 (LMAN), NO, CMI is too high.
% pu69 -2(RA, ~), 

% L, R
% 13, 14 - not good
% 19, 18 - not good.

% ---- LMAN
BirdToPlot = 'pu69wh78';
NeurToPlot = [1]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
motiflist = {'a(b)h', 'j(b)h'};
plotbytime = 0; % links rasters for all motifs by time of song.
motifpredur = 0.25;
motifpostdur = 0.3;
plotIndivRaster = 0; % one raster for each neuron/motif
plotCombRast = 1; % one figure, all rasters
plotSmFR = 1; % all smoothed FR.
LearnKeepOnlyBase = 1;
plotSongSpec = 1;
Ntrialtoplot = 10; % will take random subset
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)


% ---- RA
NeurToPlot = [4]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)



%% ========== wh44

% 71 72 -- 73, 75, 74

% =================== LMAN
BirdToPlot = 'wh44wh39';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = [71]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
motiflist = {'j(n)h', 'a(n)h'};
% motiflist = {'a(r)dabh', 'r(r)dabh'};
plotbytime = 0; % links rasters for all motifs by time of song.

motifpredur = 0.25;
motifpostdur = 0.3;

plotIndivRaster = 0; % one raster for each neuron/motif
plotCombRast = 1; % one figure, all rasters
plotSmFR = 1; % all smoothed FR.

LearnKeepOnlyBase = 1;
plotSongSpec = 1;

Ntrialtoplot = 10; % will take random subset
% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)



% =============== RA
NeurToPlot = [74]; % 4 % vector (e.g. [5 7]) - if [] then plots all;

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)


%% ================= wh72

% 12, 13, 14 --- 16
% 13 or 14 good

% ========= LMAN
BirdToPlot = 'wh72pk12';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = [14]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
motiflist = {'o(b)h', 'l(b)h'};
% motiflist = {'a(r)dabh', 'r(r)dabh'};
plotbytime = 0; % links rasters for all motifs by time of song.

motifpredur = 0.2;
motifpostdur = 0.3;

plotIndivRaster = 0; % one raster for each neuron/motif
plotCombRast = 1; % one figure, all rasters
plotSmFR = 1; % all smoothed FR.

LearnKeepOnlyBase = 1;
plotSongSpec = 1;

Ntrialtoplot = 10; % will take random subset
% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)

% ============ RA
NeurToPlot = [16]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)


%% gr48

% 1 -- 5, arr, rrr

% 10, 11, rrr, ok

% 17 -- 16, 18, not good

BirdToPlot = 'gr48bu5';
% % ---- give it either
% A) one neuron and a bunch of motifs or
% B) bunch of neurons and one motif
NeurToPlot = [1]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
motiflist = {'a(r)r', 'r(r)r'};
% motiflist = {'a(r)dabh', 'r(r)dabh'};
plotbytime = 0; % links rasters for all motifs by time of song.

motifpredur = 0.2;
motifpostdur = 0.3;

plotIndivRaster = 0; % one raster for each neuron/motif
plotCombRast = 1; % one figure, all rasters
plotSmFR = 1; % all smoothed FR.

LearnKeepOnlyBase = 1;
plotSongSpec = 1;

Ntrialtoplot = 10; % will take random subset
% --- 1) directed song
PlotDirSong = 0; % 0 is only UNDIR, 1 is only DIR; 2 is both

lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)


% ========================== RA
NeurToPlot = [5]; % 4 % vector (e.g. [5 7]) - if [] then plots all;
lt_neural_v2_DIAGN_PlotRasterMotif(SummaryStruct, BirdToPlot, NeurToPlot, ...
    motiflist, plotbytime, motifpredur, motifpostdur, plotIndivRaster, ...
    plotCombRast, plotSmFR, PlotDirSong, LearnKeepOnlyBase, plotSongSpec, ...
    Ntrialtoplot)
