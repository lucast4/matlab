function [DatAll, Tall, Chanlist] = lt_neural_QUICK_Segextr_GetRawNeural(segextract, fnamebase, Chanlist, ...
    motifpredur, extrapad)
%% NOTE: This is actually slower than the old version (see note below), since this must load once for each trial...
%% lt 11/13/18 - extracts raw neural from origianl data
% NOTE: this is different from lt_neural_QUICK_GetRawNeural which extracts
% by going into the concatenated data file. that is probably much slower if have to
% do it one for each motif ...

% INPUTS:
%
% Chanlist =
%
%     14    21
%
% fnamebase =
%
%     '/bluejay5/lucas/birds/pu69wh78/NEURAL/111317_RALMANOvernightLearn1'

% OTUPUT:
% DatAll; channel x trial, cell array holding LFP
% Tall, timebins. (relative to syl onset, i.e. uses motifpredur and
% extrapad)


% NOTE: 
% makes sure all dat comes out with same duration by figuring out first
% what the length of vector ought to be, and then padding with last value
% if the length is one shorter than that.

%% ===== first sort chanlist
Chanlist = sort(Chanlist);

%% 

ntrials = length(segextract);
DatAll = cell(length(Chanlist), ntrials);
Tall = cell(1, ntrials);

for tt=1:ntrials
    
    fnamethis = segextract(tt).song_filename;
    tons = segextract(tt).WithinSong_TokenOns; % onset within song
    tons = tons-motifpredur; % take into account pre time desired.
    
    motifdur = segextract(tt).global_offtime_motifInclFlank - segextract(tt).global_ontime_motifInclFlank; % duration of data
    toffs = tons + motifdur;
    
    % ------------- extra padding (to account for window size)
    tons = tons-extrapad;
    toffs = toffs+extrapad;
    motifdur = motifdur+2*extrapad;
    
    % ========= extract data
    [amplifier_data, ~, frequency_parameters, ~, ...
        ~, amplifier_channels] =...
        pj_readIntanNoGui([fnamebase '/' fnamethis]);
    fs = frequency_parameters.amplifier_sample_rate;
    
    % ========= keep matrix of time x chan(desired)
    indstoget = ismember([amplifier_channels.chip_channel], Chanlist);
    assert(all([amplifier_channels(indstoget).chip_channel]==Chanlist), 'must be same order for output...');
    dat = amplifier_data(indstoget, :);
    t = [1:size(dat,2)]./fs;
    
    % --------------- max timebins to keep
    tbins_max = floor(motifdur/(t(2)-t(1)));
    
    % ========= CUT OFF TO TIME OF MOTIF
    ind_t = find(t>=tons & t<toffs);
    
    if length(ind_t)==tbins_max-1
        % --- then ok, pad with one value
        ind_t = [ind_t ind_t(end)];
%     elseif length(ind_t)==tbins_max+1
%         % -- then remove end value
%         ind_t(end) = [];
%     else
%         assert(length(ind_t)==tbins_max, 'safas');
    end
        
    dat = dat(:, ind_t(1:tbins_max));
    t = t(ind_t(1:tbins_max));
    
    % ------- convert t to time rel syl onset
    t = t-t(1)+(t(2)-t(1));
    t = t - (motifpredur + extrapad);
    
    % ========= SAVE OUTPUT
    for j=1:size(dat,1)
        DatAll{j, tt} = single(dat(j, :));
    end
    Tall{tt} = t;
end
