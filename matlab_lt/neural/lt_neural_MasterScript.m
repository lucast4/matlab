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
filename='bk7_160809_102808.rhd';
ChansToPlot.DigChans_zero=[]; % make string "all" to plot all that exist. empty array to ignore
ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
% ChansToPlot.AmpChans_zero=[9 14 19];
ChansToPlot.AmpChans_zero=[14];
% ChansToPlot.AmpChans_zero=[10 14 18 23];

% neuralFiltLow=500;
neuralFiltLow=300;

PlotWhat.raw=0;
PlotWhat.filt=1;
PlotWhat.rect_sm=0;
PlotWhat.raster=0;

Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
Raster.ThrXNoise=6; % threshold for units, used for all channels, uses absolute data for peak detection
Raster.PosOrNeg=-1; % -1 or +1, for crossings.
lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster)



%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ANALYSIS PIPELINE
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ==================================================== PREPROCESSING
%% label songs using modified evsonganaly
clear all; close all;
% batchf='batch_test';
batchf='BatchTest';
batchf='BatchChan14LateTmp';
batchf='BatchChan14good';
batchf='BatchChan14good_v3';
batchf='BatchChan14good_v4';

%% ==== exploratory - concat all audio and neural and plot for each neural channel
close all;
ChansToPlot=[14];
batchtoplot='BatchTest';
% batchtoplot=batchf;
lt_neural_concatExplore(batchtoplot, ChansToPlot)


%% ==== concatenate multiple files for given channel, [and saves]
% based on expectation of duration for single unit.
% -- saves into downstream folder

channel_board=14;
lt_neural_concatOneChan(batchf, channel_board)

%% ==== run wave_clus on this concatted data


%% ==== [SANITY CHECK] plot song files in entirety, aligned to extracted spikes
% -- makes multiple plots if too much dat.
close all;
plotcols={'m', 'r','c', 'b', 'g'};
lt_neural_AlgnWavclus(batchf, channel_board, plotcols);


%% =================================================== ANALYSIS 
%% ==== EXTRACT SONG, LABEL,SongDat ONSETS, SPIKE DATA
close all;
[SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board);


%% ===== Extract song segments that match the desired regexp 
% along with associated neural data
close all;

% - desired motifs
regexpr_str='g(h)h'; % token determines where align. predur is relative to token.
regexpr_str='n(h)h'; % token determines where align. predur is relative to token.
regexpr_str='[^h](h)h'; % token determines where align. predur is relative to token.
predur=1;
postdur=2;

% - entire motifs
regexpr_str='(S)[\w-]*?E'; % gets any S ---- E
predur=6; % sec
postdur=8; % sec

[SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, regexpr_str, predur, postdur);

%% ===== [OPTIONAL] extract segments automatically - motif on and off


%% ===== uniformly stretch durations of segments
% assumes that onset of syl 1 to offset of syl last should be aligned.
close all;
TargetDur=[]; 
[SegmentsExtract, Params] = lt_neural_LinTimeWarp(SegmentsExtract, Params, TargetDur);


%% ==== plot rasters
close all;
useRescaled=0; % 1, then need to run LinTimeWarp first (plots scaled spikes, not song dat)
plotAllSegs=0; % then plots each trial own plot.
[Params]=lt_neural_motifRaster(SegmentsExtract, Params, useRescaled, plotAllSegs);


%% === compare rasters and mean firing rate for differnet motifs
% will align them by their own respective alignment locations.
% - can either rescale them to all be the same dur, or not rescale.

close all; 
% MotifList_regexp={'(S)[\w-]*?E'};
% MotifList_regexp={'gh(h)', 'nh(h)'};
MotifList_regexp={'v(b)b', '[^vb](b)b'};
MotifList_regexp={'n(b)b', 'v(b)b'};

% MotifList_regexp={'[^b](b)', 'b(b)'};
MotifList_regexp={'(g)b', '(g)h'};
predur=0.5;
postdur=1;
LinearWarp=1; % 0=no; 1=yes, each motif individually; 2=yes, all motifs to global motif median dur
suppressplots=1; % diagnostic plots.
onlyPlotMean=0; % then does not plot rasters, just means.
lt_neural_motifRastMult(SongDat, NeurDat, Params, MotifList_regexp, ...
    predur, postdur, LinearWarp, suppressplots, onlyPlotMean)



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







