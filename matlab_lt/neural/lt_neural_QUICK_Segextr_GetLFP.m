function [LFPall, Tall, Chanlist] = lt_neural_QUICK_Segextr_GetLFP(segextract, fnamebase, Chanlist, ...
    motifpredur, extrapad)
%% lt 10/29/18 - given segextact, gets lfp aligned to each trial. need to have first extracted lfp...
%
% Chanlist =
%
%     14    21
%
% fnamebase =
%
%     '/bluejay5/lucas/birds/pu69wh78/NEURAL/111317_RALMANOvernightLearn1'

% OTUPUT:
% LFPall; channel x trial, cell array holding LFP
% Tall, timebins. (relative to syl onset, i.e. uses motifpredur and
% extrapad)



% NOTE: 
% makes sure all dat comes out with same duration by figuring out first
% what the length of vector ought to be, and then padding with last value
% if the length is one shorter than that.

%%


%%
ntrials = length(segextract);
LFPall = cell(length(Chanlist), ntrials);
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
    lfpstruct = load([fnamebase '/' fnamethis '.lfp'], '-mat');
    
    % ========= keep matrix of time x chan(desired)
    chansall = lfpstruct.lfpstruct.chanlist;
    indstoget = ismember(chansall', Chanlist);
    dat = lfpstruct.lfpstruct.dat(:, indstoget);
    assert(all(lfpstruct.lfpstruct.chanlist(indstoget)' == Chanlist), 'output chans are not in exact same order as input... (either dont have preextract for all chans or other error');
    t = lfpstruct.lfpstruct.t;
    
    % --------------- max timebins to keep
    tbins_max = motifdur/(t(2)-t(1));
    
    % ========= CUT OFF TO TIME OF MOTIF
    ind_t = find(t>=tons & t<=toffs);
    
    if length(ind_t)==tbins_max-1
        % --- then ok, pad with one value
        ind_t = [ind_t ind_t(end)];
    end
        
    dat = dat(ind_t(1:tbins_max), :);
    t = t(ind_t(1:tbins_max));
    
    % ------- convert t to time rel syl onset
    t = t-t(1)+(t(2)-t(1));
    t = t - (motifpredur + extrapad);
    
    % ========= SAVE OUTPUT
    for j=1:size(dat,2)
        LFPall{j, tt} = dat(:,j);
    end
    Tall{tt} = t;
end