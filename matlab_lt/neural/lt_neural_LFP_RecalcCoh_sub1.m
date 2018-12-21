function [xcorr_max, xcorr_lag] = lt_neural_LFP_RecalcCoh_sub1(dat1, dat2, ...
    onlyOnePeriod, freqsthis, fs, plotindiv, usecorrcoeff)

%%


ntrials = size(dat1,3);
nfreq = size(dat1, 2);
xcorr_max = nan(nfreq, ntrials);
xcorr_lag = nan(nfreq, ntrials);
for nt = 1:ntrials
    for nf = 1:nfreq'
        if onlyOnePeriod==1
            maxlag = fix(0.75*(1500/freqsthis(nf)));
        end
        if usecorrcoeff==1
        [c, lags] = xcorr(dat1(:, nf, nt), dat2(:, nf, nt), maxlag, 'coeff');            
        else
        [c, lags] = xcorr(dat1(:, nf, nt), dat2(:, nf, nt), maxlag, 'unbiased');
        end
        
        % -- find peak
        [~, locs] = findpeaks([0; c; 0],  'SortStr', 'descend', 'NPeaks', 1); % find max peak
        locs = locs-1;
        
        xcorr_max(nf, nt) = c(locs);
        xcorr_lag(nf, nt) = lags(locs)/fs; % conver to sec
        
        % --- sanity check
        if rand<0.02 & plotindiv==1
            if usecorrcoeff==1
                [ctmp, lagstmp] = xcorr(dat1(:, nf, nt), dat2(:, nf, nt), maxlag, 'coeff');
            else
        [ctmp, lagstmp] = xcorr(dat1(:, nf, nt), dat2(:, nf, nt), maxlag, 'biased');
            end
            lt_figure; hold on;
            
            lt_subplot(2,2,1); hold on;
            xlabel('blue = 1');
            plot(dat1(:, nf, nt), 'b');
            plot(dat2(:, nf, nt), 'r');
            
            lt_subplot(2,2,2); hold on;
            title('unbiased (using this) ][dash = biased][');
            xlabel('blue leads -- red leads');
            plot(lags, c, '-k');
            plot(lags(locs), c(locs), 'sr');
            plot(lagstmp, ctmp, '--k');
            
            %
            %                                     lt_subplot(2,2,4); hold on;
            %                                     xlabel('blue leads -- red leads');
            %                                     title('biased')
            % %                                     plot(lags(locs), c(locs), 'sr');
            
        end
    end
end
