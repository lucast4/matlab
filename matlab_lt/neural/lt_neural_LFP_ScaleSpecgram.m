function [S_base, S_WN] = lt_neural_LFP_ScaleSpecgram(S_base, S_WN, scaletype)
%% lt 2/27/19 - takes matrix of spectrograms (many trials). rescales in different ways
% spectrogram should be power.



%                     S_base = S1mat(:,:, indsbase_epoch, j);
%                     S_WN = S1mat(:,:, indsWN_epoch, j);



%% RUN

if scaletype ==1
    % --- do zscore for each t,ff bin (relative to baseline
    % distribution)
    [S_base, S_WN] = fn_zscorespec(S_base, S_WN);
elseif scaletype ==2
    % --- for each slice and trial normalize as power
    % proportion (i.e. divide by sum across
    % frequencies)
    S_base = lt_neural_LFP_NormSpec(S_base);
    S_WN = lt_neural_LFP_NormSpec(S_WN);
elseif scaletype ==3
    % each ff window z-scored
    [S_base, S_WN] = fn_zscorespec_withinFF(S_base, S_WN);
elseif scaletype==4
    [S_base, S_WN] = fn_zscorespec_withinFF(S_base, S_WN);
    S_base = lt_neural_LFP_NormSpec(S_base);
    S_WN = lt_neural_LFP_NormSpec(S_WN);
end

end

%% functions
function [S_base_z, S_WN_z] = fn_zscorespec(S_base, S_WN)
% will zscore each t,ff bin based on baseline
% input sould be power (i.e not log). code will take log for you.

% -- take log [makes distrubtions more gaussina]
S_base = 10*log10(S_base);
S_WN = 10*log10(S_WN);

% get baseline mean and std
S_mean = mean(S_base,3);
S_std = std(S_base, [], 3);

% perform zscoring
S_base = S_base - repmat(S_mean, 1, 1, size(S_base,3));
S_base = S_base./repmat(S_std, 1, 1, size(S_base,3));
S_base_z = S_base;

S_WN = S_WN - repmat(S_mean, 1, 1, size(S_WN,3));
S_WN = S_WN./repmat(S_std, 1, 1, size(S_WN,3));
S_WN_z = S_WN;
end

function [S_base_z, S_WN_z] = fn_zscorespec_withinFF(S_base, S_WN)
% will zscore each t,ff bin based on baseline
% input sould be power (i.e not log). code will take log for you.

% -- take log [makes distrubtions more gaussina]
S_base = 10*log10(S_base);
S_WN = 10*log10(S_WN);

% get baseline mean and std
nffbin = size(S_base,2);
Smeans = nan(1, nffbin);
Sstds = nan(1, nffbin);
for n=1:nffbin
    tmp = S_base(:,n,:);
    s_mean = mean(tmp(:));
    s_std = std(tmp(:));
    
    Smeans(n) = s_mean;
    Sstds(n) = s_std;
end


S_base = S_base - repmat(Smeans, size(S_base,1), 1, size(S_base,3));
S_base_z = S_base./repmat(Sstds, size(S_base,1), 1, size(S_base,3));

S_WN = S_WN - repmat(Smeans, size(S_WN,1), 1, size(S_WN,3));
S_WN_z = S_WN./repmat(Sstds, size(S_WN,1), 1, size(S_WN,3));
end