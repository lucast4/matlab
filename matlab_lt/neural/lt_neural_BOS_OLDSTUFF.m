%% ====================== [OLD] EXTRACT UNCLUSTERED NEURAL DATA
% -- neural
chanstoplot = [9 14 17];


% --- 1) extract timestamps of digital signals indicating bos playback
fid = fopen(batchfile);

fline = fgetl(fid);

% --------- for debugging, plotting all pulse data
figcount=1;
subplotrows=6;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

while ischar(fline)
    
    % ---------- load file
    [amplifier_data,board_dig_in_data,frequency_parameters, ...
        board_adc_data, board_adc_channels, amplifier_channels, ...
        board_dig_in_channels, t_amplifier] = pj_readIntanNoGui(fline);
    
    fs = frequency_parameters.amplifier_sample_rate;
    
    % ================== extract syl pulses
    if strcmp(chantype, 'dig')
        ind = [board_dig_in_channels.chip_channel]==pulsechan;
        pulsedat = board_dig_in_data(ind, :);
    elseif strcmp(chantype, 'ana')
        ind = [board_adc_channels.chip_channel]==pulsechan;
        pulsedat = board_adc_data(ind, :);
    end
    
    % ---------------------
    if (1)
        % --- plot pulse dat
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        
        plot(pulsedat);
        line(xlim, [threshold threshold], 'Color','r');
    end
    
    
    % ================================ IS THERE BOS IN THIS FILE? HOW MANY?
    % GET ONSETS/OFFSETS
    risetime
    
    
    
    % ====================== extract sound
    songdat = board_adc_data([board_adc_channels.chip_channel]==audchan, :);
    %     t= [1:length(songdat)]./fs;
    %     lt_plot_spectrogram(songdat, fs, 0, 0);
    %     plot(songdat, 'm');
    plot((songdat-median(songdat)).^2);
    
    % ================ extract neural
    if isempty(chanstoplot)
        chans = [amplifier_channels.chip_channel];
    else
        chans = chanstoplot;
    end
    
    
    for i=chans
        neurdat = amplifier_data([amplifier_channels.chip_channel]==i, :);
        
        
        
    end
    
    
    
    fline = fgetl(fid);
end

%% ============================= EXTRACT PRE-CLUSTERED NEURAL DATA
% ---

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

