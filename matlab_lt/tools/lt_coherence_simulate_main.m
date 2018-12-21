function [Y, C,phi,S12,S1,S2,f] = lt_coherence_simulate_main(freqlist, amplist, ntrials, fs, twin)

t = linspace(1/fs, twin, fs*twin);

Y = {};

% === signal 1
indthis = 1;
% --- multiple trials
y_alltrials = [];
for tt = 1:ntrials
    y = lt_coherence_simulate_sub(freqlist, amplist, indthis, t);
    y_alltrials = [y_alltrials y'];
end
Y{indthis} = y_alltrials;

% === signal 2
indthis = 2;
% --- multiple trials
y_alltrials = [];
for tt = 1:ntrials
    y = lt_coherence_simulate_sub(freqlist, amplist, indthis, t);
    y_alltrials = [y_alltrials y'];
end
Y{indthis} = y_alltrials;


%% ==== coherence between these signals
movingwin = [0.1 0.01];
params = struct;
params.fpass = [1/movingwin(1) 200];
w = 30; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
% w = 20; % in hz, for desired frequency resolution of tapers. % note, t is set to movingwin(1)
tw = movingwin(1)*w;
params.tapers = [tw 2*tw-1];
params.Fs = 1500; % hard coded fs for LFP;

lt_switch_chronux(1);
[C,phi,S12,S1,S2,f]=coherencyc(Y{1},Y{2},params);

lt_switch_chronux(0);
