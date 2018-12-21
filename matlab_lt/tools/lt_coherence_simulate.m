%% generate random signal composed of sum of sinusoids.

freqlist{1} = [1 5 30 50];
freqlist{2} = [1 5 30 75];
amplist{1} = [1 1 0 1];
amplist{2} = [1 1 0 1];
ntrials = 50;

% === generate signals
fs = 1500;
twin = 0.1; % seconds.


[Y, C,phi,S12,S1,S2,f] = lt_coherence_simulate_main(freqlist, amplist, ntrials, fs, twin);
lt_figure; hold on; 
plot(f, C, '-k');

