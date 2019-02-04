%% ########## take entier motif (stereotyped) and plot things. also compares RA and LMAN

% == TO DO:
% trial by trial xcov.
% also subtract shuffle corrector.


%% ========== 1) extract full motifs and time warp

lt_neural_v2_ExtractFullMotifs;


%% ==== load latest
clear all; close all;
load('/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/MOTIFSTATS_Compiled_31Mar2018_1153.mat');
load('/bluejay0/bluejay2/lucas/analyses/neural/MOTIFSTATS_Compiled/Params_31Mar2018_1153.mat');

%% ======= removed directed song
MOTIFSTATS_Compiled = lt_neural_QUICK_MotCom_RemoveDIR(MOTIFSTATS_Compiled);


%% TIME WARP ALL DATA (FOR A GIVEN MOTIF, SAME WARP ACROSS ALL NEURONS)
% close all;

MOTIFSTATS_Compiled = lt_neural_MOTIF_TimeWarpAll(MOTIFSTATS_Compiled);

%% === 2) plot all motifs, with smoothed firing, running std, for all cells
close all;
% birdtoplot = 'wh44wh39';
% birdtoplot = 'pu69wh78';
% birdtoplot = 'gr48bu5';
birdtoplot = 'wh72pk12';
% motiftoplot = '(n)hh';
% motiftoplot = '(j)jjbhhg';
% motiftoplot = '(a)jjbh';
motiftoplot = '(m)kd';
plotcv = 0; % if 1, then plots running cv. if 0 then plots fr and running mean.
bregionsToPlot = {'RA', 'LMAN'};

doshuffle = 1; % shuffles timing of spikes [the dumbest version]
whichhalf = []; % 1=early, 2 = late; [] = all data (i.e. motif split in half...)
onlywithinsyl = 1; % then splices out to only keep data during song (i.e. concatenates).


plotRaw=0;
lt_neural_FullMotif_SmFR(MOTIFSTATS_Compiled, birdtoplot, ...
    motiftoplot, plotcv, bregionsToPlot, plotRaw, doshuffle, whichhalf);


%% === FOR A GIVEN NEURON, PLOT RASTERS
close all;
i=2; % birdnum
mm = 2; % motif num
neurid = 11;
maxtrials = 50;

% ============= RUN
segextract = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.neurons(neurid).motif(mm).SegmentsExtract;

bregionthis = SummaryStruct.birds(i).neurons(neurid).NOTE_Location;
birdname = SummaryStruct.birds(i).birdname;
motifname = MOTIFSTATS_Compiled.birds(i).MOTIFSTATS.params.motif_regexpr_str{mm};

lt_figure; hold on;
title([birdname '-' motifname '-n' num2str(neurid) '[' bregionthis ']']);

if length(segextract)>maxtrials
    trialstoplot = randperm(length(segextract), maxtrials);
else
   trialstoplot = 1:length(segextract);
end

for t = 1:length(trialstoplot)
    tt = trialstoplot(t);
    spktimes = segextract(tt).spk_Times;
    assert(all(segextract(tt).spk_Clust == SummaryStruct.birds(i).neurons(neurid).clustnum));
    
    lt_neural_PLOT_rasterline(spktimes, t, 'k');
end

% ======== OVERLAY SYL ONSETS/OFFESTS
onsets = segextract(1).motifsylOnsets;
offsets = segextract(1).motifsylOffsets;

axis tight
YLIM = ylim;

for j=1:length(onsets)
    line([onsets(j) offsets(j)], [3 3], 'Color', 'm', 'LineWidth', 4);
    line([onsets(j) offsets(j)], [YLIM(2)-4 YLIM(2)-4], 'Color', 'm', 'LineWidth', 4);
    line([onsets(j) onsets(j)], ylim, 'Color', 'm');
    line([offsets(j) offsets(j)], ylim, 'Color', 'm');
end


%% #################################################
%% ############################### AUTOCORRELATIONS

close all;
lt_neural_AutoCorr_Calc(MOTIFSTATS_Compiled);
