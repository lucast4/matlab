%% LT 11/4/17 - rough analyses on populations


%% go thru all chans/clusts, extract data
clear all; close all;

% =============== PU69
% Batchname = 'BatchTest'; % must all have used the same batch
% ChansToGet = [9 14 17 18 21];
% ClustToGet = [1 1 1 1 1]; % should be aligned with ChansToGet
% BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'RA'};
% dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/NOTSINGING';

% Batchname = 'Batch1219to1304'; % must all have used the same batch
% ChansToGet = [9 14 17 21];
% ClustToGet = [1 1 1 1]; % should be aligned with ChansToGet
% BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'RA'};
% dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/110117_RALMANlearn1/NOTSINGING';

% Batchname = 'Batch1208to2205'; % must all have used the same batch
% ChansToGet = [17 18 21];
% ClustToGet = [1 1 1]; % should be aligned with ChansToGet
% dirname = '/bluejay5/lucas/birds/pu69wh78/NEURAL/103017_RAlearn1/';
% 

% ================================= WH44
% ----- pre WN
% Batchname = 'Batch1146to1158'; % must all have used the same batch
% ChansToGet = [14 15 17 18 20 20];
% ClustToGet = [1 1 1 1 1 2]; % should be aligned with ChansToGet
% BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'RA', 'RA'};
% dirname = '/bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2/LIGHTSOFF';

% ----- during WN
% Batchname = 'Batch1718to1734'; % must all have used the same batch
% ChansToGet = [14 15 17 20 20 22];
% ClustToGet = [1 1 1 1 2 1]; % should be aligned with ChansToGet
% BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'RA', 'RA'};
% dirname = '/bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2/LIGHTSOFF';

% --- Sleep
Batchname = 'Batch0009to0037'; % must all have used the same batch
ChansToGet = [15 15 17 17 14 18 20 21 22];
ClustToGet = [1 2 1 2 1 1 1 1 1]; % should be aligned with ChansToGet
BrainRegion = {'LMAN', 'LMAN', 'RA', 'RA', 'LMAN', 'RA', 'RA', 'RA', 'RA'};
dirname = '/bluejay5/lucas/birds/wh44wh39/NEURAL/031418_RALMANlearn2/NIGHT';



% -- defaults
extractsound = 0;

NeurDat = struct;
numsampsAll = []; % to confirm that all data are aligned (samme num samps);
fsAll = [];
filenamesAll = {};

for i = 1:length(ChansToGet)
    cd(dirname);
    chan = ChansToGet(i);
    clust = ClustToGet(i);
    
    [SongDatTMP, NeurDatTMP, ParamsTMP] = lt_neural_ExtractDat(Batchname, chan, ...
        extractsound, clust);
    
    % - convert from ms to sec
    NeurDatTMP.spikes_cat.cluster_class(:,2) = ...
    NeurDatTMP.spikes_cat.cluster_class(:,2)./1000;
        
    
    % -
    NeurDat.metaDat = NeurDatTMP.metaDat;
    NeurDat.unit(i).spikes_cat = NeurDatTMP.spikes_cat;
    NeurDat.unit(i).chan = chan;
    NeurDat.unit(i).clust = clust;
    NeurDat.dirname = dirname;
    NeurDat.unit(i).BrainRegion = BrainRegion{i};
    
    % ====== sanity checks - all data are from same files
    numsampsAll = [numsampsAll; [NeurDatTMP.metaDat.numSamps]];
    fsAll = [fsAll; [NeurDatTMP.metaDat.fs]];
    filenamesAll = [filenamesAll; [NeurDatTMP.metaDat.filename]];
end

% === sanity checks
assert(size(unique(numsampsAll, 'rows'),1) ==1, 'problem');
assert(size(unique(fsAll, 'rows'),1)==1, 'problem');
assert(numel(unique(filenamesAll)) ==1, 'problem');

%%

NumUnits = length(NeurDat.unit);

%% PLOT ALL CHANNELS

% ============== SPIKES
lt_figure; hold on;
title('all simult. units');

chantoplot = []; % empty for all;

if ~isempty(chantoplot)
NeurNums = find(ChansToGet==chantoplot);
else
    NeurNums = 1:length(ChansToGet);
end

for i=NeurNums
    
    spktimes = NeurDat.unit(i).spikes_cat.cluster_class(:,2);
    spktimes = spktimes;
    if strcmp(NeurDat.unit(i).BrainRegion, 'LMAN')
        plotcol = [0.2 0.4 0.4];
    elseif strcmp(NeurDat.unit(i).BrainRegion, 'RA')
        plotcol = [0.8 0.3 0.3];
    end
    lt_neural_PLOT_rasterline(spktimes, i, plotcol);    
end


% =============== VOLTAGE (filtered)
lt_neural_PLOT_RawMultChan(NeurDat.dirname, Batchname, ChansToGet);


%% =========================== PAIRWISE CROSS-COR BETWEEN ALL PAIRS OF CHANS

% ============= 1) CROSS-CORRELATION OF BINNED SPIKES
binsize = 0.005; % in sec. shift will be identical
windowmax = 0.2; % in sec, extent of cross-corr (+/- windowmax)

% ----- 1) bin all spike trains
fs = NeurDat.metaDat(1).fs;
nsamps = NeurDat.metaDat(1).numSamps;
trialdur = nsamps/fs; % sec
binedges = 0:binsize:trialdur;

FRmat = nan(length(binedges)-1, length(NeurDat.unit)); % timebin x channel
for i=1:length(NeurDat.unit)
    
    x = histc(NeurDat.unit(i).spikes_cat.cluster_class(:,2), binedges); % binned spikes   
    
    FRmat(:,i) = x(1:end-1);
end

% ----- 2) calcualte all cross-corrs
lt_figure; hold on;
% [CCall, lags] = xcov(FRmat, ceil(windowmax/binsize), 'coeff');
% [CCall, lags] = xcorr(FRmat, ceil(windowmax/binsize));

NeurID1 = [];
NeurID2 = [];
CCall = [];
for i=1:length(NeurDat.unit)
    for ii=i:length(NeurDat.unit)
        
        [cc, lags] = xcov(FRmat(:,i), FRmat(:,ii), ceil(windowmax/binsize), 'coeff');
        %         [cc, lags] = xcorr(FRmat(:,i), FRmat(:,ii), ceil(windowmax/binsize));
        
        % ---- flip CC so that alphabetically earlier brain region is to
        % left
        [~, indtmp] = sort({NeurDat.unit(i).BrainRegion, NeurDat.unit(ii).BrainRegion});
        if indtmp(1)>indtmp(2)
            % -- then need to flip
            cc = flipud(cc);
            CCall = [CCall cc];
            NeurID1 = [NeurID1 ii];
            NeurID2 = [NeurID2 i];
            
        else
            CCall = [CCall cc];
            NeurID1 = [NeurID1 i];
            NeurID2 = [NeurID2 ii];
            
        end
        
    end
end

    
% ------ 3) plot CC, separating by type of pair
figcount=1;
subplotrows=8;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ------------ a) autocorrelations
inds = find(NeurID1==NeurID2);
plotcol = 'k';

for indtmp=inds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    
    
    % -- plot
    plot(lags.*binsize, CCall(:,indtmp), 'Color', plotcol);
    
    % -- label traces   
    neur1 = NeurID1(indtmp);
    neur2 = NeurID2(indtmp);
    
    label = [num2str(ChansToGet(neur1)) '-' num2str(ChansToGet(neur2))];
    title(label);
%     lt_plot_text(lags(end-1).*binsize, CCall(end-1, indtmp), label, 'r');
    
    xlim([-windowmax windowmax]);
%     ylim([-0.3 0.8]);
    xlabel([BrainRegion{neur1} '<---->' BrainRegion{neur2}])
    line([0 0], ylim, 'LineStyle', '--');
end



% -------------- b) same region
inds = find(strcmp(BrainRegion(NeurID1), BrainRegion(NeurID2)) & NeurID1~=NeurID2);
plotcol = 'b';

for indtmp=inds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    
    % -- plot
    plot(lags.*binsize, CCall(:,indtmp), 'Color', plotcol);
    
    % -- label traces   
    neur1 = NeurID1(indtmp);
    neur2 = NeurID2(indtmp);
    
    label = [num2str(ChansToGet(neur1)) '-' num2str(ChansToGet(neur2)) '; ' BrainRegion{neur1} '-' BrainRegion{neur2}];
    title(label);
%     lt_plot_text(lags(end-1).*binsize, CCall(end-1, indtmp), label, 'r');
    
    xlim([-windowmax windowmax]);
%     ylim([-0.3 0.8]);
    xlabel([BrainRegion{neur1} '<---->' BrainRegion{neur2}])
    line([0 0], ylim, 'LineStyle', '--');
end



% -------------- b) diff region
inds = find(~strcmp(BrainRegion(NeurID1), BrainRegion(NeurID2)));
plotcol = 'r';

ccall = []; % for plotting mean;
for indtmp=inds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
  
    % -- plot
    plot(lags.*binsize, CCall(:,indtmp), 'Color', plotcol);
    
    % -- label traces   
    neur1 = NeurID1(indtmp);
    neur2 = NeurID2(indtmp);
    
    label = [num2str(ChansToGet(neur1)) '-' num2str(ChansToGet(neur2)) '; ' BrainRegion{neur1} '-' BrainRegion{neur2}];
    title(label);
%     lt_plot_text(lags(end-1).*binsize, CCall(end-1, indtmp), label, 'r');
    
    xlim([-windowmax windowmax]);
%     ylim([-0.3 0.8]);
    xlabel([BrainRegion{neur1} '<---->' BrainRegion{neur2}])
    
    % --- collect
    ccall = [ccall CCall(:,indtmp)];
    line([0 0], ylim, 'LineStyle', '--');
end
%-- plot mean
ccmean = mean(ccall,2);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean');
hsplots = [hsplots hsplot];

plot(lags.*binsize, ccmean, 'Color', plotcol);
    xlim([-windowmax windowmax]);
%     ylim([-0.3 0.8]);
        line([0 0], ylim, 'LineStyle', '--');


linkaxes(hsplots, 'x');


%% ========= calculate coherency between channels

chan1 = 9;

% load file
pj_readIntanNoGui('pu69wh78_171101_102913.rhd');

%%

indtmp = [amplifier_channels.chip_channel] == chan1;
x1 = amplifier_data(indtmp, :)';


params = struct;
params.Fs = 30000;
params.fpass = [0 1000];
% params.err = [2 0.05];
% [S, f, Serr] = mtspectrumc(x1, params);
[S, f] = mtspectrumc(x1);

lt_figure; hold on;
plot(f, S, '-ok');
% shadedErrorBar(f, S, Serr, {'Color', 'k'},1);

%% ###############################################
%% ############################################### LFP, SINGING
%%  extract raw dat for each channel
% here working thru single song files
clear all; close all;
lt_switch_chronux(1);
sfile = '/bluejay5/lucas/birds/wh44wh39/NEURAL/021518_Rec3/wh44wh39_180215_105438.rhd';
% ChansToPlot = [15 17 18 20 21 22];
% BrainRegions = {'LMAN', 'RA','RA','RA','RA', 'RA'};
ChansToPlot = [15 17 21];
BrainRegions = {'LMAN', 'RA', 'RA'};

% ============ 1) extract dat
[amplifier_data,board_dig_in_data,frequency_parameters, board_adc_data, ...
    board_adc_channels, amplifier_channels, board_dig_in_channels, t_amplifier] =...
    pj_readIntanNoGui(sfile);
fs = frequency_parameters.amplifier_sample_rate;

% --- for each channel extract raw dat
DatAll = struct;
for j = 1:length(ChansToPlot)
    
    chan = ChansToPlot(j);
    
    DatAll(j).datraw = amplifier_data([amplifier_channels.chip_channel] == chan, :)';
    DatAll(j).chan_amp = chan;
    DatAll(j).bregion = BrainRegions{j};
end


%% ====== SPECTROGRAM OF LFP
% ---- for each region get spectrogram and plot
% for j=1:length(DatAll)

datall = [DatAll.datraw];
params = struct;
params.Fs = fs;
params.fpass = [5 150];
movingwin = [0.1 0.005];


[S, t, f] = mtspecgramc(datall, movingwin, params);

% ====================================  plot [spectrograms]
figcount=1;
subplotrows=length(DatAll)+1;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];
% --- first plot sound spectrogram
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
lt_plot_spectrogram(board_adc_data, fs, 1, 0);

% --- then plot each neural spectrogram
for j=1:length(DatAll)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['ch' num2str(DatAll(j).chan_amp) '[' DatAll(j).bregion ']']);
    imagesc(t, f, 10*log10(S(:,:, j)'));
    axis tight;
    ylim([1/(movingwin(1)) params.fpass(2)]);
end
linkaxes(hsplots, 'x');
colormap('jet');
        
%% ======================= SPECTRUM
chantmp = 1;
spectrum = mean(S(:,:,chantmp),1);
lt_figure; hold on;
plot(f, 10*log10(spectrum), '-k');


%% ====================================  plot [raw unfiltered dat]
figcount=1;
subplotrows=length(DatAll)+1;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];
% --- first plot sound spectrogram
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
lt_plot_spectrogram(board_adc_data, fs, 1, 0);

% --- then plot each neural spectrogram
for j=1:length(DatAll)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['ch' num2str(DatAll(j).chan_amp) '[' DatAll(j).bregion ']']);
t = [1:length(DatAll(j).datraw)]./fs;
dat = DatAll(j).datraw;

% --- filter neural data
dat = lt_neural_filter(dat, fs, 0, 1, 300);
plot(t, dat, '-k');
dat = lt_neural_filter(dat, fs, 0, 25, 45);
plot(t, dat, '-r');
axis tight;
end
linkaxes(hsplots, 'x');


%% ============================= COHEROGRAM
PairToGet = {'LMAN', 'RA'}; 
        params = struct;
        params.Fs = fs;
        params.fpass = [5 1000];
        params.tapers = [3 5];
        movingwin = [0.15 0.01];

        DatCoh = struct;
for i=1:length(DatAll)
   
    region1 = DatAll(i).bregion;
    
    if ~strcmp(region1, PairToGet{1})
        continue
    end
    
    % --- go thru all other channels
    for ii=1:length(DatAll)
    
        region2 = DatAll(ii).bregion;
        
        if ~strcmp(region2, PairToGet{2})
            continue
        end
        
        disp([region1 '-' region2]);
        
       % ------- get coherence
       
       dat1 = DatAll(i).datraw;
       dat2 = DatAll(ii).datraw;
       
       [C,phi,S12,S1,S2,t,f] = cohgramc(dat1, dat2, movingwin, params);
       
       DatCoh(i,ii).C = C;
       DatCoh(i,ii).phi = phi;
       
       DatCoh(i,ii).t = t;
       DatCoh(i,ii).f = f;
       
       
        
    end
        
end

%% ============================ PLOT COHERENCE FOR ALL PAIRS
figcount=1;
subplotrows=6;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];
% --- first plot sound spectrogram
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
lt_plot_spectrogram(board_adc_data, fs, 1, 0);

% --- then plot each coherogram
for j=1:length(DatCoh(:))
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
% title(['ch' num2str(DatAll(j).chan_amp) '[' DatAll(j).bregion ']']);
dat = DatCoh(j).C;
t = DatCoh(j).t;
f = DatCoh(j).f;
imagesc(t, f, dat');
axis tight;
ylim([1/(movingwin(1)) params.fpass(2)]);
colormap('jet');
end
linkaxes(hsplots, 'x');


%% ============================= COHERENCE (SPECTRUM) DURING SONG
% ============== collect part of data that is during song




%%
lt_switch_chronux(0);

%% ###############################################
%% ========================================== overnight analysis
batchf = 'batchtmp';
chantype = 'ana'; % dig or ana
chan = 2; % 0, 1, ...
threshold = 1; % voltage, or rise and fall detection.
% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchf);

fline = fgetl(fid);
cumusamps = 0; 
rtAll = [];
ftAll = [];
datAll = [];
while ischar(fline)

    % ==================== extract signal
    [~,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels] = pj_readIntanNoGui(fline);

        
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==chan;
        dat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==chan;
        dat = board_adc_data(ind, :);
    end
    
    
    
    % --- detect threshold crossings
    rt = risetime(dat, 'StateLevels', [0 threshold]);
    ft = falltime(dat, 'StateLevels', [0 threshold]);
    
    % --- rise and fall times as cumu time
    rt = rt+cumusamps;
    ft = ft+cumusamps;
    
    % -- collect across all files
    rtAll = [rtAll rt];
    ftAll = [ftAll ft];
    
    if (1)
    datAll = [datAll dat];
    end
    
    % --- cumulateive time up to now
    cumusamps = cumusamps + length(dat);
    
    fline = fgetl(fid);
end



%% ############### plot batch songs, compare dig pulse with song
close all;

batchf = 'batchall';
chantype = 'ana'; % dig or ana
chan = 0; % 0, 1, ... % for pulse
audchan = 1;
threshold = 1.5; % voltage, or rise and fall detection.
% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchf);

fline = fgetl(fid);
cumusamps = 0; 
rtAll = [];
ftAll = [];
datAll = [];

figcount=1;
subplotrows=8;
subplotcols=1;
fignums_alreadyused=[];
hfigs=[];


hsplots = [];
while ischar(fline)

    % ==================== extract signal
    [~,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels] = pj_readIntanNoGui(fline);

        fs = frequency_parameters.amplifier_sample_rate;
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==chan;
        dat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==chan;
        dat = board_adc_data(ind, :);
    end

    songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
    
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
hsplots = [hsplots hsplot];
title(fline);

t= [1:length(songdat)]./fs;
plot(t, songdat,'k');
plot(t, dat, 'r');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hold on;
hsplots = [hsplots hsplot];

lt_plot_spectrogram(songdat, fs, 0, 0);

%     % --- detect threshold crossings
%     rt = risetime(dat, 'StateLevels', [0 threshold]);
%     ft = falltime(dat, 'StateLevels', [0 threshold]);
%     
%     % --- rise and fall times as cumu time
%     rt = rt+cumusamps;
%     ft = ft+cumusamps;
%     
%     % -- collect across all files
%     rtAll = [rtAll rt];
%     ftAll = [ftAll ft];
%     
%     if (1)
%     datAll = [datAll dat];
%     end
%     
%     % --- cumulateive time up to now
%     cumusamps = cumusamps + length(dat);
%     
linkaxes(hsplots, 'x');
    fline = fgetl(fid);
    hsplots = [];
end
