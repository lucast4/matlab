function [FRmat, t] = lt_neural_QUICK_Spk2FRmat(Spks, mintime, maxtime)
%% lt 12/2/18 - from Spks, outputs matrix of smoothed FR (over trials)
% NOTE: lt_neural_motifRaster has old version, which is running boxcar mean over trials.

%%

% Spks is cell array, each cell spktimes (trials x 1)
% mintime is in sec, is timepoint data starts
% maxtime is in sec.

assert(size(Spks,2)==1);

if size(Spks{1},1)>1
    % then flip
    Spks = cellfun(@transpose, Spks, 'UniformOutput', 0);
end

%% first comvert spks cell array to segextract structure (for compativitlty with older code)

segextract = struct;

for j=1:length(Spks)
    
    % -- make sure spks start at 0. needed for the code below. will then
    % adjust to make sure times match input
    segextract(j).spk_Times = Spks{j}-mintime;
    segextract(j).datdur = maxtime-mintime;
    
end


%%

segextract = lt_neural_SmoothFR(segextract, [], [], [], 0, [], mintime);

FRmat = [segextract.FRsmooth_rate_CommonTrialDur];
t = segextract(1).FRsmooth_xbin_CommonTrialDur;

assert(size(FRmat,1) == length(t), 'if not then need to take transpoe within trials')

% === sanity check
if (0)
    % --- actual spikes
    spkthis = Spks{2};
    
    % --- smoothed
    x = segextract(2).FRsmooth_xbin_CommonTrialDur;
    y = segextract(2).FRsmooth_rate_CommonTrialDur;
    
    figure; hold on; 
    plot(x,y);
    plot(spkthis, 1, 'or');
    
end

