%% the last params on 4/5/19. not exactly sure what they are for...

onlygoodexpt = 1;
% xtoplot = [-0.12 0.02];
% xtoplot = [-0.14 0.01];
xtoplot = [-0.1 0.0];
xtoplot = [-0.1 0.0];
plotraw = 0; % 1: plots each trial (examples...)
disp('NOTE: remember to include wh72 in this analysis');
zscore_lfp = 0; % if 1, then z-scores concatenating all trials (a given chan) within time segment (after filtering)
% fpass = [40 100]; % for bandpass filtering LFP.
% fpass = [20 35]; % for bandpass filtering LFP.
% fpass = [5 350]; % for bandpass filtering LFP.
fpass = [25 60]; % for bandpass filtering LFP.
% sta_wind = [-0.05 0.05]; % relative to spike, in sec % will only get spikes that are within data...
sta_wind = [-0.03 0.03]; % relative to spike, in sec % will only get spikes that are within data...
% kernelSD = 0.015; % empyt for default (for spike smoothgin)
kernelSD = 0.005; % empyt for default (for spike smoothgin)
removeBadLFP = 1; % then things that look like large fluctuations..
