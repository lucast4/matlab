function [filtdat, tout] = lt_neural_QUICK_Segextr_GetFiltBankv2(segextract, ...
    motifpredur, extrapad, FILTMATall, songfiles,tall)
%% v2 means first extract all files acorss motifs and trials, and then enter here to get all trials for a given motif.
% FILTMATall must already have desired f and chans extracted (but all t remian)

%% lt 12/14/18 - modified from LFP extract to now extract filtered LFP for freq range of inrterst.

% == OUT:
% filtdat = cell (trial x 1), each cell is , t x f x chans (mat);
% t seconds, relative to syl onset (alignement)
% freqsout and chansout are in correct order, matdched to data

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

% freqvals = [20:5:40]; [has to alread be extracted file by file]


% NOTE: 
% makes sure all dat comes out with same duration by figuring out first
% what the length of vector ought to be, and then padding with last value
% if the length is one shorter than that.

%%


%%
ntrials = length(segextract);
filtdat = cell(ntrials, 1);
tout = cell(ntrials, 1);

for tt=1:ntrials
    disp(tt);
    
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
    indthis = strcmp(songfiles, fnamethis);
    filtmatthis = FILTMATall{indthis};
    t = tall{indthis};
    assert(sum(indthis)==1);
         
    
    % ========= keep matrix of time x chan(desired) x freq(desired)
    % -- time...
    ind_t = find(t>=tons & t<=toffs);
    tbins_max = fix(motifdur/(t(2)-t(1))); % max timebins to keep
    if length(ind_t)==tbins_max-1
        % --- then ok, pad with one value
        ind_t = [ind_t ind_t(end)];
    elseif length(ind_t)>tbins_max
        ind_t = ind_t(1:tbins_max);
    end

    % ================================ KEEP TIME OF INTEREST    
    filtmatthis = filtmatthis(ind_t, :, :);
    t = t(ind_t);
    t = t-t(1)+(t(2)-t(1));
    t = t - (motifpredur + extrapad);    % ------- convert t to time rel syl onset

    
    % ============== OUTPUT
    filtdat{tt} = filtmatthis;
    tout{tt} = t;
end