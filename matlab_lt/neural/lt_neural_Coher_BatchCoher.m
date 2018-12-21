function [C, t, f, phi,S12,S1,S2] = lt_neural_Coher_BatchCoher(lfpmat1, lfpmat2, version, t_LFP, ...
    ntapers, movingwin, tw)
%% lt 12/9/18 - calculates across trial coherence

% lfpmat, samples x trials

% version = 'mtaper_all'; % trials x tapers number of datapoints..
% version = 'mtaper_trials'; % one coh for each trial, then takes average
% version = 'welch_trials'; % each trial get coher over multiple subwindws
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
elseif strcmp(version, 'welch_trials')
    window = floor(movingwin(1)*params.Fs);
    % --- slide through manuall.
    ntot = size(lfpmat1,1); % total samples.
    shift = movingwin(2)*params.Fs;
    ons = 1:shift:ntot-window+1;
    offs = window:shift:ntot;
    windlist = [ons' offs'];
    % --- within each window, use these params (i.e. each is  welchs segment)
    windsmall = floor(window/2);
    olap = floor((3/4)*windsmall);
    nfft = max([256 2^(nextpow2(windsmall))]);
    % == iterate over this list
    Cgram = [];
    t = [];
    for j=1:size(windlist,1)
        on = windlist(j,1);
        off = windlist(j,2);
        [Cxy, f] = mscohere(lfpmat1(on:off,:), lfpmat2(on:off,:), windsmall, olap, ...
            nfft, 1500);
        % -- get average over trials
        Cxy = mean(Cxy, 2);
        Cgram = [Cgram Cxy];
        t = [t; mean([on off])];
    end
    t = t'./params.Fs;
    t = t+t_LFP(1);
    C = Cgram';
    
    indtmp = f>=params.fpass(1) & f<=params.fpass(2);
    C = C(:, indtmp);
    f = f(indtmp)';
    
    phi = [];
    S12 = [];
    S1 = [];
    S2 = [];
end

lt_switch_chronux(0);
