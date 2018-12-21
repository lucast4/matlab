function [C, t, f] = lt_neural_Cohere_Recalcsub(thingstodo, lfp1_tmp, ...
    lfp2_tmp, cohversion, ntapers, movingwin, tw, t_LFP)
%% gets one mean for coh across trials
C = [];
t = [];
f = [];

% ====== COHERENCE
if any(strcmp(thingstodo, 'cohere'))
    [C, t, f, phi,S12,S1,S2] = lt_neural_Coher_BatchCoher(lfp1_tmp, lfp2_tmp, cohversion, ...
        t_LFP, ntapers, movingwin, tw);    
end

% ===== LFP XCORR
if any(strcmp(thingstodo, 'lfpxcorr'))
    
    
end

% ==== WAVELETS
if any(strcmp(thingstodo, 'waveletcoh'))
    % ======== calculate cohernec for each trial
    ntrial = size(lfp1_tmp, 2);
    for tt=1:ntrial
        [wC, ~, f] = wcoherence(lfp1_tmp(:,tt), lfp2_tmp(:,tt), 1500);
        wC = fliplr(wC');
        f = flipud(f);
        if tt==1
            C = nan(size(wC,1), size(wC,2), ntrial);
        end
        C(:,:,tt) = wC;
    end
    C = mean(single(C),3);
    t = [];
    f = f';
end

