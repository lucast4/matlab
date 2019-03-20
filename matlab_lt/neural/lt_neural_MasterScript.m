%% TO DO:
% 1) save new time base for amplitude plotting (lin warped)
% 2)


%% ===== plot song

songname='br92br54_170416_111151.rhd';
% songnotmat='bk7_160809_110736.rhd.not.mat';

lt_neural_extractaudio(songname);

% ---- overlay onsets and offsets [first label with evsonganaly]
notmatstruct=load(songnotmat);

for i=1:length(notmatstruct.onsets)
    
    line([notmatstruct.onsets(i) notmatstruct.offsets(i)]./1000, [1.65 1.65], 'Color','r','LineWidth', 7);
end


%% get a list of sampling rates of all files


%% ====== plot single file dat [align neural and song]
close all;
filename='wh44wh39_180324_171611.rhd';
ChansToPlot.DigChans_zero=[0]; % make string "all" to plot all that exist. empty array to ignore
ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
% ChansToPlot.AmpChans_zero=[9 14 19];
ChansToPlot.AmpChans_zero=[9 11 12 14 1 18];
% ChansToPlot.AmpChans_zero=[8 9 11 16 17 20 21];

ChansToPlot.AmpChans_zero=[9 14 17 18 21];
ChansToPlot.AmpChans_zero=[9 14 21];
ChansToPlot.AmpChans_zero=[14 15 17 18 20 21];
% ChansToPlot.AmpChans_zero=16:31;
ChansToPlot.AmpChans_zero=0:15;
% ChansToPlot.AmpChans_zero=[8 10 15];
ChansToPlot.AmpChans_zero=0:31;
ChansToPlot.AmpChans_zero=[14 15 18 20 21 22];

% neuralFiltLow=500;
neuralFiltLow=300;

PlotWhat.lfp = 0;
PlotWhat.raw=0;
PlotWhat.filt=1;
PlotWhat.rect_sm=0;
PlotWhat.raster=0;
PlotWhat.digital=0;

Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
Raster.ThrXNoise=3; % threshold for units, used for all channels, uses absolute data for peak detection
Raster.PosOrNeg=-1; % -1 or +1, for crossings.

numsubplots = 6;

lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster, ...
    numsubplots)

%% ======= clean song files




%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ANALYSIS PIPELINE
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%% ==================================================== VARIOUS THINGS
%% PCA to try to remove common noise

% - SAVES ONE CHANNEL, READY FOR WAVE_CLUS (HERE SAVES ALREADY BANDPASS
% FILTERED, WHILE WAVECLUS ITSELF ALSO RUNS A FILTER)
% - DENOISES THAT CHANNEL CASED ON LINEAR REGRESSION AGAINST THE OTHER
% CHANNELS + CHANNEL MEDIAN (TO IGNORE SPIKES - I.E. IF NOISE THEN ALL
% CHANNELS SHOULD BE DOING SAME THING)


close all;
% ChansToUse=[9 15 18 11 14];
ChansToUse=[9 15 18 11];
batchf='Batch1710to1925';
plotOn = 1;
PlotAllChanSep = 0; % then makes one plot for each chan, comparing diff methods

onlyPlotFinalRegressionFig = 1; % overwrites other plot params, only plots final summary fig.
skipPCAmethod = 1; % then only does regression method - useful if just want to run fast and save
save_raw_dat = 1;
ChanToSave = 9;
Subsample_NumSec = 160; % number of seconds to use for building model.

lt_neural_commonNoise(batchf, ChansToUse, plotOn, save_raw_dat, ...
    PlotAllChanSep, ChanToSave, onlyPlotFinalRegressionFig, skipPCAmethod, Subsample_NumSec)

% NOTE (12/2/16) - currently regression method worked well for one set of
% songs and PCA method for another. Neither method worked great.

% NOTE (2/22/17) - regression method seems to work well. An issue that I am
% worried about is introducing spikes where there should not be any. Is not
% easy to confirm that this method is not producing artifacts. So for now
% will avoid using this. But by eye looks like dramatically reduces one of
% the sources of noise, but there is another source of noise in many files
% that is not removed, and is actually enhanced in some cases - i.e.
% different weights depending on source.



%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==================================================== PREPROCESSING
clear all; close all; fclose all;
% channel_board = [8 11 14 18 20]; % wh6
% channel_board = [9 11 12 14 15 18]; % bu77
% channel_board = [14];
channel_board = [9 14 17 18 21];
channel_board = 0:31;
channel_board = [8 9 14 21];
channel_board = [20];
batchf = 'batchtmp';

%% ==== exploratory - concat all audio and neural and plot for each neural channel
close all;

% -----  v1, plots just raw neural, filtered
% batchtoplot=batchf;
if (0)
    lt_neural_concatExplore(batchf, channel_board);
end

% ----- v2, plots filtered neural, smoothed, and spike waveforms
PlotRectDat=0; % 1, plots, 0 skips.
PlotFiltDat=1; % usually 1, filt neural.
PosAndNeg =0; % then gets both. if 0, then just downwards
PlotLFP = 0; % takes precedence...

lt_neural_concatExplore_v2(batchf, channel_board, PlotRectDat, PlotFiltDat, ...
    PlotLFP, PosAndNeg); % WAVEFORMS ONLY PLOTTED FOR ONE CHANNEL!!

% ------- v3 -plots song by song
if (0)
    fid = fopen(batchf);
    filename = fgetl(fid);
    while ischar(filename)
        
        
        ChansToPlot.DigChans_zero=[0]; % make string "all" to plot all that exist. empty array to ignore
        ChansToPlot.AnalogChans_zero=[0]; % assumes that this is audio
        % ChansToPlot.AmpChans_zero=[9 14 19];
        ChansToPlot.AmpChans_zero=channel_board;
        % ChansToPlot.AmpChans_zero=[8 9 11 16 17 20 21];
        % ChansToPlot.AmpChans_zero=[8 9 16 17 20 21];
        
        % neuralFiltLow=500;
        neuralFiltLow=300;
        
        PlotWhat.raw=0;
        PlotWhat.filt=1;
        PlotWhat.rect_sm=0;
        PlotWhat.raster=0;
        PlotWhat.digital=0;
        
        Rect_sm.windowsize=0.03; % in sec, size of window, equals -2 to +2 sd.
        Raster.ThrXNoise=3; % threshold for units, used for all channels, uses absolute data for peak detection
        Raster.PosOrNeg=-1; % -1 or +1, for crossings.
        
        numsubplots = 2;
        
        lt_neural_alignRawDat(filename, ChansToPlot, neuralFiltLow, PlotWhat, Rect_sm, Raster, numsubplots)
        
        disp([' ======== ' filename]);
        
        pause;
        filename = fgetl(fid);
    end
    fclose all;
end


%% ==== concatenate multiple files for given channel, [and saves]
% based on expectation of duration for single unit.
% -- saves into downstream folder
close all;
lt_neural_concatOneChan(batchf, channel_board)
cd(['Chan' num2str(channel_board) 'amp-' batchf]);


%% ==== run wave_clus on this concatted data



%% ======= CREATE NOT.MAT. FILE FOR ALL SONG FILES IN BATCH [AUTO]

lt_neural_AutoMakeNotmat(batchf);


%% ==== [SANITY CHECK] plot song files in entirety, aligned to extracted spikes
% -- makes multiple plots if too much dat.
close all;
PlotSecondChan = 1;
SecondChan = 15;
plotcols={'m', 'r','c', 'b', 'g'};

% want to plot 2nd channel to compare noise?
if (0)
    lt_neural_AlgnWavclus(batchf, channel_board, plotcols, PlotSecondChan, SecondChan);
end

% === VERSION 2 - PLOTS EACH SONG FILE ONE AT A TIME - CAN CLOSE THEM ONE
% AT A TIME
maxfiguuresopen = 20;
figsstepsize = 5; % num figs to jump thru, useful if many song. 1 is default. DOES NOT WORK
lt_neural_AlgnWavclus_v2(batchf, channel_board, plotcols, PlotSecondChan, SecondChan, maxfiguuresopen, figsstepsize);


%% ===== HAND REMOVE NOISY SPIKES
% IMPORTANT!!: only run this once you are satisfied with wave-clus cluster
% assignments. If reassign clusters in wave_clus, then will delete all
% hand-removed clusters (cluster -1)
% NOTE: can check which spikes removed using lt_neural_AlgnWavclus_v2

close all;

SecondChan = 18;
lt_neural_CheckClust(channel_board, batchf, SecondChan);


% ========== TROUBLESHOOTING
% i.e. made mistake of hand checking before finalized (consolidated)
% clusters. THen run below (MODIFIED!) to manually change cluster numbers,
% whithout wiping out hand checked (cluster -1)

if (0)
    tmp = load('times_data_TMPBACKUP.mat');
    
    % to convert cluster 2 to 0
    inds = tmp.cluster_class(:,1)==2;
    tmp.cluster_class(inds,1) = 0;
    
    % to convert clusters 3 and 4 to 1
    inds = tmp.cluster_class(:,1) == 3 | tmp.cluster_class(:,1) == 4;
    tmp.cluster_class(inds,1) = 1;
    
    % to save
    ver = '-v6';
    save('times_data.mat', '-struct', 'tmp' , ver);
    disp('SAVED');
end
%% +++++++++++++++++++++++++++++++++++++== NOTES THINGS TO UPDATE WITH ANALYSIS PIPELINE

% save all things to a common folder. Levels:
% bird,
% batch - metadata (real times of files. names of files, etc).
%     do not save neural data. [if eventually need this can easily write
%     script to get]
%     do not save song data
% labels. onsets, offsets. [write script to update all labels]. or just use
% labels from notmat files
% spike times
% PROCESSED DATA: e.g. pitch ...




%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
%% ++++++++++ ANALYSIS (NEURONS ONE BY ONE)
%% ==== EXTRACT SONG, LABEL,SongDat ONSETS, SPIKE DATA
close all; clear all;
batchf='Batch0922to1251';
channel_board=14;
extractsound = 1;
[SongDat, NeurDat, Params] = lt_neural_ExtractDat(batchf, channel_board, extractsound);


%% ===== Extract song segments that match the desired regexp
% along with associated neural data
close all;

% - desired motifs
% regexpr_str='[^vb](b)'; % token determines where align. predur is relative to token.
% regexpr_str='(g)h'; % token determines where align. predur is relative to token.
% regexpr_str='n(h)h'; % token determines where align. predur is relative to token.
% regexpr_str='[^h](h)h'; % token determines where align. predur is relative to token.
regexpr_str='nk(h)'; % token determines where align. predur is relative to token.
predur=0.25;
postdur=0.1;
alignByOnset=1;

% - entire motifs
% regexpr_str='(S)[\w-]*?E'; % gets any S ---- E
% predur=6; % sec
% postdur=8; % sec
% alignByOnset=1;

% - Extract song motifs (onset and offset) automatically)
% regexpr_str='WHOLEBOUTS';
% predur=8; % sec
% postdur=3; % sec
% alignByOnset=0;

keepRawSongDat = 1;
[SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
    regexpr_str, predur, postdur, alignByOnset, [], [], keepRawSongDat);



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


%% #######################################################################
%% ######################## LOOK AT LFP
clear all; close all;
fname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/pu69wh78_171101_164255.rhd';
% fname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/pu69wh78_171101_164255.rhd';
chanstoget = [9 14 21];
chansbregions = {'LMAN', 'LMAN', 'RA'};
freqwind = [25 25];
dozscore = 1; % for lfp amplitude and coerehce.
disp('NOTE: by default opnly gets coehrence btween LAMN and RA');



% ################################ load file

[amplifier_data,board_dig_in_data,frequency_parameters, board_adc_data, ...
    board_adc_channels, amplifier_channels, board_dig_in_channels, t_amplifier] =...
    pj_readIntanNoGui(fname);

% ############################### RUNA ND PLOT

figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ==== 1) LOAD AUDIO
dataudio = board_adc_data(1,:);
fs = frequency_parameters.amplifier_sample_rate;
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title(fname);
lt_plot_spectrogram(dataudio, fs, 1, 0);


% ====== IONE COLOR FOR EACH CHANEL
pcols = lt_make_plot_colors(length(chanstoget), 0,0);


% ====== 1) RAW NEURAL (UNFILTERED)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('raw neural');
indschan = ismember([amplifier_channels.chip_channel], chanstoget);
assert(all([amplifier_channels(indschan).chip_channel] == chanstoget));
datneur = amplifier_data(indschan, :);
tbin = 1:length(datneur);
tbin = tbin./fs;
for i=1:size(datneur,1)
    plot(tbin, datneur(i,:), '-', 'Color', pcols{i});
    disp(i)
end


% ===== 2) PLOT LFP (UNFILTERED)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('lfp (0-400hz');
hsplots = [hsplots; hsplot];
datlfp = load([fname, '.lfp'], '-mat');

indschan = ismember(datlfp.lfpstruct.chanlist, chanstoget);
datlfp2 = datlfp.lfpstruct.dat(:, indschan);
assert(all(datlfp.lfpstruct.chanlist(indschan) == chanstoget'));
tbin = datlfp.lfpstruct.t;

for i=1:size(datlfp2,2)
    y = datlfp2(:,i);
    plot(tbin, y, '-', 'Color', pcols{i});
end


% ==== 2) PLOT LFP FOR EACH CHANNEL
% [AND DIFF FILTERED BANDS]
datfilt = load([fname '.filt'], '-mat');
indsthis = ismember(datfilt.filtstruct.chanlist, chanstoget);
assert(all(datfilt.filtstruct.chanlist(indsthis) == chanstoget'));

indsf = datfilt.filtstruct.freqvals>=freqwind(1) & datfilt.filtstruct.freqvals<=freqwind(2);
fbins = datfilt.filtstruct.freqvals(indsf);

datfilt2 = cellfun(@(x)x(indsf), datfilt.filtstruct.datfilt_chans(indsthis), 'UniformOutput', 0);
datfilt2 = cellfun(@cell2mat, cellfun(@transpose, datfilt2, 'UniformOutput', 0), 'UniformOutput', 0);
% datfilt2 = cell2mat(tmp);

tbins = linspace(datfilt.filtstruct.t_edges(1), datfilt.filtstruct.t_edges(2), size(datfilt2{1},1));


% - one plot for each frequency
for i=1:length(fbins)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([num2str(fbins(i)) ' hz']);
    % === once for each chan
    for cc=1:length(chanstoget)
        chanthis = chanstoget(cc);
        y = datfilt2{cc}(:, i);
        plot(tbins, real(y), '-', 'Color', pcols{cc});
    end
end
% --- note down chans/bregions relative to color
for cc=1:length(chanstoget)
    chanthis = chanstoget(cc);
    bregionthis = chansbregions{cc};
    y = datfilt2{cc}(:, i);
    lt_plot_text(tbins(end), real(y(end)), ['ch' num2str(chanthis) ',' bregionthis], pcols{cc});
    %         plot(tbins, real(y), '-', 'Color', pcols{cc});
end


% ==== 3) PLOT POWER (AMPLITUDE OF HILBERT TRANSFORM)
% - one plot for each frequency
for i=1:length(fbins)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([num2str(fbins(i)) ' hz']);
    ylabel('amplitude');
    % === once for each chan
    for cc=1:length(chanstoget)
        chanthis = chanstoget(cc);
        y = datfilt2{cc}(:, i);
        %         y = abs(y).^2;
        y = abs(y).^2;
        if dozscore==1
            y = (y-mean(y))./std(y);
        end
        plot(tbins, real(y), '-', 'Color', pcols{cc});
%         end
    end
end


% ==== 3) PLOT COHERENCE
% --- pairwise coherence for all pairs
lt_switch_chronux(1);
movingwin = [0.1 0.01];
params = struct;
params.fpass = [1/movingwin(1) 150];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title([num2str(fbins(i)) ' hz']);
ylabel('coherence');

for i=1:size(datlfp2,2)
    for ii=i+1:size(datlfp2,2)
        bregion1 = chansbregions{i};
        bregion2 = chansbregions{ii};
        if strcmp(bregion1, 'LMAN') & strcmp(bregion2, 'RA')
            dat1 = datlfp2(:, i);
            dat2 = datlfp2(:, ii);
            chan1 = chanstoget(i);
            chan2 = chanstoget(ii);
        elseif strcmp(bregion1, 'RA') & strcmp(bregion2, 'LMAN')
            dat1 = datlfp2(:, ii);
            dat2 = datlfp2(:, i);
            chan1 = chanstoget(ii);
            chan2 = chanstoget(i);
        else
            % then skip./
            continue
        end
        
        [C,phi,S12,S1,S2,t,f] = cohgramc(dat1, dat2, movingwin, params);
        colthis = [rand rand rand];
        
        % == take mean over those freqency b ands
        if freqwind(2) ==freqwind(1)
            % then get closest single band
            [~, indf] = min(abs(f-freqwind(1)));
            fthis = f(indf);
            y = mean(C(:, indf),2);
        else
            indf = f>freqwind(1) & f<freqwind(2);
            y = mean(C(:, indf),2);
        end
        
        if dozscore==1
            y = (y-mean(y))./std(y);
        end
        plot(t, y, '-', 'Color', colthis);
        lt_plot_text(t(end), y(end), ['ch:' num2str(chan1) '-' num2str(chan2)], colthis);
        
        % ====== smooth
        tsm = lt_running_stats(t, round(1/(t(2)-t(1))));
        ysm = lt_running_stats(y, round(1/(t(2)-t(1))));
        shadedErrorBar(tsm.Mean, ysm.Mean, ysm.STD, {'Color', colthis}, 1);
        
    end
end
lt_switch_chronux(0);


linkaxes(hsplots, 'x');


