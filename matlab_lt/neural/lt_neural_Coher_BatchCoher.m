function [C, t, f, phi,S12,S1,S2] = lt_neural_Coher_BatchCoher(lfpmat1, lfpmat2, version, t_LFP, ...
    ntapers, movingwin, tw)
%% lt 12/9/18 - calculates across trial coherence

% lfpmat, samples x trials

% version = 'mtaper_all'; % trials x tapers number of datapoints..
% version = 'mtaper_trials'; % one coh for each trial, then takes average
% over trials.

% t_LFP, time base for LFP;'

%%
% === chronux p[arams
lt_switch_chronux(1);
if isempty(movingwin)
movingwin = [0.1 0.01];
end

params = struct;
params.fpass = [1/movingwin(1) 150];
if isempty(tw)
tw = 3;
end
w = tw/movingwin(1); % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
% tw = movingwin(1)*w;
if isempty(ntapers)
    k = 2*tw-1;
else
    k = ntapers;
end
params.tapers = [tw k];
params.Fs = 1500; % hard coded fs for LFP;
params.trialave = 1;

if strcmp(version, 'mtaper_all')
    params.trialave = 1;
    [C,phi,S12,S1,S2,t,f] = cohgramc(lfpmat1, lfpmat2, movingwin, params);
    % [C,phi,S12,S1,S2,t,ffbins] = cohgramc(dat1, dat2, movingwin, params);
    %                     t = t - PARAMS.motif_predur - PARAMS.xtrapad;
    t = t+(t_LFP(1));
elseif strcmp(version, 'mtaper_trials')
    params.trialave = 0;
    [C,phi,S12,S1,S2,t,f] = cohgramc(lfpmat1, lfpmat2, movingwin, params);
    C = mean(C, 3);
    t = t+(t_LFP(1));
end

lt_switch_chronux(1);
