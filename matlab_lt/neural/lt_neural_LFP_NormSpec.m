function specmat = lt_neural_LFP_NormSpec(specmat)
%% lt 11/6/18 - % --- for each slice and trial normalize as power proportion
% (i.e. divide by sum across
% frequencies)

% specmat, t x ff x trials.

ntrials = size(specmat,3);
for t=1:ntrials

    
    datthis = specmat(:,:,t);
    
    norm = repmat(sum(datthis, 2), 1, size(datthis,2));
    
    datthis = datthis./norm;
    
    if (0)
        figure; hold on;
        plot(1:size(datthis,2), datthis', '-');
    end
        
    specmat(:,:, t) = datthis;
end
