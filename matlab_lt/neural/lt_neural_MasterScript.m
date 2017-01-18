%% TO DO:
% 1) save new time base for amplitude plotting (lin warped)
% 2) 


%% ===== plot song

songname='bk7_160809_110736.rhd';
songnotmat='bk7_160809_110736.rhd.not.mat';

lt_neural_extractaudio(songname);

% ---- overlay onsets and offsets [first label with evsonganaly]
notmatstruct=load(songnotmat);

for i=1:length(notmatstruct.onsets)
    
    line([notmatstruct.onsets(i) notmatstruct.offsets(i)]./1000, [1.65 1.65], 'Color','r','LineWidth', 7);
end


%% get a list of sampling rates of all files


%% ====== plot single file dat [align neural and song]
close all;
filename='r87gr4_4_20_2015_DV1200_150420_150330.rhd';
ChansToPlot.DigChans_zero=[0]; % make string "all" to plot all that exist. empty array to ignore
ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
% ChansToPlot.AmpChans_zero=[9 14 19];
% ChansToPlot.AmpChans_zero=[9];
ChansToPlot.AmpChans_zero=[8 9 11 16 17 20 21];
ChansToPlot.AmpChans_zero=[8 9 16 17 20 21];

% neuralFiltLow=500;
neuralFiltLow=300;

PlotWhat.raw=0;
PlotWhat.filt=1;
PlotWhat.rect_sm=0;
PlotWhat.raster=0;
PlotWhat.digital=0;

Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
Raster.ThrXNoise=6; % threshold for units, used for all channels, uses absolute data for peak detection
Raster.PosOrNeg=-1; % -1 or +1, for crossings.
lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster)

%% ======= clean song files




%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ANALYSIS PIPELINE
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ==================================================== VARIOUS THINGS
%% PCA to try to remove common noise

close all;
ChansToUse=[8 9 12 15 17 19 20 22];
ChansToUse=[8 12 19 20 22];
batchf='batch_epoch3';
plotOn = 0;
save_raw_dat = 0;
PlotAllChanSep = 1; % then makes one plot for each chan, comparing diff methods

lt_neural_commonNoise(batchf, ChansToUse, plotOn, save_raw_dat, PlotAllChanSep);

% NOTE (12/2/16) - currently regression method worked well for one set of
% songs and PCA method for another. Neither method worked great. 



%% ==================================================== PREPROCESSING
%% label songs using modified evsonganaly
clear all; close all;
% batchf='batch_test';
batchf='BatchChan14Late';

%% ==== exploratory - concat all audio and neural and plot for each neural channel
close all;
ChansToPlot=[8 9 12 15 17 19 20 22];
batchtoplot='batchall_sub';

% -----  v1, plots just raw neural, filtered
% batchtoplot=batchf;
lt_neural_concatExplore(batchtoplot, ChansToPlot);

% ----- v2, plots filtered neural, smoothed, and spike waveforms
PlotRectDat=0; % 1, plots, 0 skips.
PlotFiltDat=1; % usually 1, filt neural.

lt_neural_concatExplore_v2(batchtoplot, ChansToPlot, PlotRectDat, PlotFiltDat); % WAVEFORMS ONLY PLOTTED FOR ONE CHANNEL!!

%% ==== concatenate multiple files for given channel, [and saves]
% based on expectation of duration for single unit.
% -- saves into downstream folder
close all; clear all;
batchf='batchall_sub';
channel_board=20;
lt_neural_concatOneChan(batchf, channel_board)

%% ==== run wave_clus on this concatted data


%% ==== [SANITY CHECK] plot song files in entirety, aligned to extracted spikes
% -- makes multiple plots if too much dat.
close all;
plotcols={'m', 'r','c', 'b', 'g'};
lt_neural_AlgnWavclus(batchf, channel_board, plotcols);


%% ++++++++++ ANALYSIS (NEURONS ONE BY ONE) 
%% ==== EXTRACT SONG, LABEL,SongDat ONSETS, SPIKE DATA
close all;
batchf='BatchChan14good_v4';
channel_board=14;
[SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);


%% ===== Extract song segments that match the desired regexp 
% along with associated neural data
close all;

% - desired motifs
% regexpr_str='[^vb](b)'; % token determines where align. predur is relative to token.
% regexpr_str='(g)h'; % token determines where align. predur is relative to token.
% regexpr_str='n(h)h'; % token determines where align. predur is relative to token.
% regexpr_str='[^h](h)h'; % token determines where align. predur is relative to token.
regexpr_str='h(h)'; % token determines where align. predur is relative to token.
predur=0.1;
postdur=0.1;
alignByOnset=1;

% - entire motifs
regexpr_str='(S)[\w-]*?E'; % gets any S ---- E
predur=6; % sec
postdur=8; % sec
alignByOnset=1;

% - Extract song motifs (onset and offset) automatically)
regexpr_str='WHOLEBOUTS';
predur=8; % sec
postdur=3; % sec
alignByOnset=0;

[SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, regexpr_str, predur, postdur, alignByOnset);



%% ===== [OPTIONAL] extract segments automatically - motif on and off


%% ===== uniformly stretch durations of segments
% assumes that onset of syl 1 to offset of syl last should be aligned.
close all;
TargetDur=[]; 
[SegmentsExtract, Params] = lt_neural_LinTimeWarp(SegmentsExtract, Params, TargetDur);


%% ==== plot rasters

close all;
useRescaled=0; % 1, then need to run LinTimeWarp first (plots scaled spikes, not song dat)
plotAllSegs=1; % then plots each trial own plot.
[Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs);


%% === compare rasters and mean firing rate for differnet motifs
% will align them by their own respective alignment locations.
% - can either rescale them to all be the same dur, or not rescale.

close all; 
% MotifList_regexp={'(S)[\w-]*?E'};
% MotifList_regexp={'g(h)h', 'n(h)h'};
% MotifList_regexp={'v(b)b', '[^vb](b)b'};
% MotifList_regexp={'n(b)b', 'v(b)b'};
% MotifList_regexp={'[nv](b)', '[gm](b)'};
% MotifList_regexp={'(v)', '(b)'};
% % MotifList_regexp={'[^b](b)', 'b(b)'};
% MotifList_regexp={'(g)b', '(g)h'};
% MotifList_regexp={'gh(h)hh'};
MotifList_regexp={'g(h)', 'n(h)', 'h(h)'};

predur=0.1;
postdur=0.1;
LinearWarp=0; % 0=no; 1=yes, each motif individually; 2=yes, all motifs to global motif median dur
suppressplots=0; % diagnostic plots.
onlyPlotMean=0; % then does not plot rasters, just means.

lt_neural_motifRastMult(SongDat, NeurDat, Params, MotifList_regexp, ...
    predur, postdur, LinearWarp, suppressplots, onlyPlotMean)

% === TO DO, LIMIT TO A SINGLE CHANNEL


%% +++++++++++++++++++++++++++++++++++++++++++++++++
%% [debug] extract data for single channel of single file

filename='bk7_160809_110736.rhd';
chan=2;

[amplifier_data, ~, frequency_parameters, ...
    board_adc_data, board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(filename);

data=amplifier_data(2,:);

save('data.mat', 'data');

% ===== after using wave_clus
times_data=load('times_data.mat');
numclust=size(times_data.cluster_class,2);

% - plot this file
    ChansToPlot.DigChans_zero=[]; % make string "all" to plot all that exist. empty array to ignore
    ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
    ChansToPlot.AmpChans_zero=amplifier_channels(chan).chip_channel;
    % ChansToPlot.AmpChans_zero=[10 14 18 23];

    neuralFiltLow=500;

    PlotWhat.raw=0;
    PlotWhat.filt=1;
    PlotWhat.rect_sm=1;
    PlotWhat.raster=1;

    Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
    Raster.ThrXNoise=times_data.par.stdmin; % threshold for units, used for all channels, uses absolute data for peak detection
    Raster.PosOrNeg=-1; % -1 or +1, for crossings.
lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster)

% - overlay spikes detected from wave_clus
plotcols={'g', 'b', 'm', 'c'};
for i=1:numclust
   inds=find([times_data.cluster_class(:,1)]==i);
   
   for ii=inds
      
       spktime=times_data.cluster_class(ii, 2); % in ms
       
       line([spktime spktime], [-20 -40], 'Color', plotcols{i}, 'LineWidth',2);
   end
end




%% output
% for each song file, have a data structure containing various things







