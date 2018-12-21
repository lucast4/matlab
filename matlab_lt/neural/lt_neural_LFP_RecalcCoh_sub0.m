function [xcorr_max, xcorr_lag, freqsthis] = lt_neural_LFP_RecalcCoh_sub0(filtdatstruct, ...
    chan1, chan2, indstoget, lfpxcorr_twindtoget, lfpxcorr_ffwindtoget, ...
    onlyOnePeriod, fs, usecorrcoeff)

%%
% === 1) extract LFP for desired chan and frequencies
dat1 = cellfun(@(x)x(:,:, filtdatstruct.Chanlist==chan1), ...
    filtdatstruct.filtdat_t_f_chan(indstoget), 'UniformOutput', 0);
dat2 = cellfun(@(x)x(:,:, filtdatstruct.Chanlist==chan2), ...
    filtdatstruct.filtdat_t_f_chan(indstoget), 'UniformOutput', 0);

% --- convert to matrices.
dat1 = lt_neural_Coher_Cell2Mat(dat1);
dat2 = lt_neural_Coher_Cell2Mat(dat2);

% --- keep only the time iwnodw of interest
indtmp = find(filtdatstruct.t_relons>=lfpxcorr_twindtoget(1) & ...
    filtdatstruct.t_relons<=lfpxcorr_twindtoget(2));
dat1 = dat1(indtmp, :,:);
dat2 = dat2(indtmp, :,:);

% --- keep only the freq iwnodw of interest
indtmp = find(filtdatstruct.freqvals>=lfpxcorr_ffwindtoget(1) & ...
    filtdatstruct.freqvals<=lfpxcorr_ffwindtoget(2));
dat1 = dat1(:, indtmp,:);
dat2 = dat2(:, indtmp,:);
freqsthis = filtdatstruct.freqvals(indtmp);

% --- get real part of hilvert transomf
dat1 = real(dat1);
dat2 = real(dat2);

% --- for each trial/ff band, compute xcorr
plotindiv = 0;
[xcorr_max, xcorr_lag] = lt_neural_LFP_RecalcCoh_sub1(dat1, dat2, ...
    onlyOnePeriod, freqsthis, fs, plotindiv, usecorrcoeff);


