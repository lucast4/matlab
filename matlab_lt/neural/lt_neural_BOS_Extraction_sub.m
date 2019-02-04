function [LFP_dat, LFP_tedges, LFP_chanlist] = ...
    lt_neural_BOS_Extraction_sub(SummaryBOS, fline, i, amplifier_data, ...
    amplifier_channels, fs, t_amplifier)
%%
[~, tmpb, ~] = fileparts(fline);

% -- do once for each channel
chanlist = unique(SummaryBOS.expt(i).channels);
% = for each channel, find raw data and filter to get lfp
LFP_dat = cell(1, length(chanlist));
LFP_tedges = [];
LFP_chanlist = chanlist;
for cc=1:length(chanlist)
    chanthis = chanlist(cc);
    dattmp = amplifier_data([amplifier_channels.chip_channel]==chanthis,:);
    
    % ================= LOWPASS
    Ntmp=4;
    fpass = 400;
    [filt_b,filt_a]=butter(Ntmp, [fpass*2/fs]);
    datfilt = filtfilt(filt_b, filt_a, dattmp'); % do multiple channels at once
    datfilt = datfilt';
    
    % ---- DOWNSAMPLE
    fs_new = 1500;
    t = t_amplifier;
    factor = floor(fs/fs_new);
    datfilt = downsample(datfilt', factor);
    t = downsample(t, factor);
    
    % ########################### FINAL STUFF.
    datfilt = single(datfilt);
    t = single(t);
    
    LFP_dat{cc} = datfilt;
    LFP_tedges = [t(1) t(end)];
    
    %             tmpdat = load([SummaryBOS.expt(i).dirname '/' SummaryBOS.expt(i).batchname '-Chan' num2str(chanthis) '/' tmpb '.ch' num2str(chanthis) '.mat']);
    %             all(amplifier_data([amplifier_channels.chip_channel]==chanthis,:) == tmpdat.data)
end
