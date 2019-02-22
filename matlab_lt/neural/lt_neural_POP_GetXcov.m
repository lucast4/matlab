function [ccRealAll, ccShiftAll, lags_sec] = lt_neural_POP_GetXcov(dattmp1, dattmp2, xcov_dattotake, ...
    motifpredur, xcovwindmax, binsize_spk, plotSummary, plotThisCC, normmethod, ...
    normspiketoprob)
%% lt 1/9/19 - copied over from code for POPLEARN, XCov analyses
% NOTE: have only looked carefully at things that will run by default (i.e
% not things like plotSummary.


%% ================== CROSS CORRELATION [uses one of 2 methods (fr or spkbins)
% -------------------- PARAMS, INITIATE
ntrials = length(dattmp1.SegmentsExtract);

ccRealAll = [];
% ccPSTH = [];
ccShiftAll =[];

% ccAuto1 = [];
% ccAuto1Shift = [];

% ccAuto2 = [];
% ccAuto2Shift = [];

if ntrials<2
    return
end

if plotThisCC==1
    lt_figure; hold on;
    ypos = 1;
    
    lt_subplot(8, 2, 1:14); hold on;
    title([bregtmp{1} '-' bregtmp{2} '[' motifstr ']-set' ...
        num2str(iii) '-n' num2str(neurons_thisset(nn)) ',' ...
        num2str(neurons_thisset(nnn))]);
end



%% ============== NEW VERSION - bins spikes, cross-correlation,
% corrected against cross corr of PSTH
%     lt_figure; hold on;
%
% % ============= 1) FIGURE OUT WINDOW FOR DATA
% if ~isempty(xcov_dattotake)
%     % shorten to desired window
%     windtmp = nan(1,2);
%     windtmp(1) = motifpredur + xcov_dattotake(1);
%     windtmp(2) = motifendtime + xcov_dattotake(2);
% end
windtmp = motifpredur + xcov_dattotake;

%% ============= 2) PERFORM XCOV
for t = 1:ntrials
    % -------- extract spks binned
    spkbin1 = dattmp1.SegmentsExtract(t).spk_Binned;
    spkbin2 = dattmp2.SegmentsExtract(t).spk_Binned;
    x = dattmp1.SegmentsExtract(t).spk_Binned_x;
    
    % --- limit to data in "premotor" window
    indtmp = x>=windtmp(1) & x<windtmp(2);
    spkbin1 = spkbin1(indtmp);
    spkbin2 = spkbin2(indtmp);
    
    % ===
    %         plot(1:length(spkbin1), spkbin1+t-0.3, '-ok');
    %         plot(1:length(spkbin2), spkbin2+t+0.3, '-or');
    
    
    if normspiketoprob ==1
        % ========= THEN MAKE IT PROB OF FINDINGS SPIKE IN BIN...
        spkbin1 = spkbin1./sum(spkbin1);
        spkbin2 = spkbin2./sum(spkbin2);
    end
    
    % ------------------- calculate
    %                                 [cc, lags] = xcov(fr1_short, fr2_short, xcovwindmax/0.001, 'coeff');
    [cc, lags] = xcorr(spkbin1, spkbin2, xcovwindmax/binsize_spk, normmethod);
    %     if all(spkbin1==0) | all(spkbin2==0)
    %         cc = zeros(size(cc)); % i.e. will give nans is using coeff.
    %     end
    %     %     [cc, lags] = xcorr(spkbin1, spkbin2, xcovwindmax/binsize_spk, ');
    %     %                                 [cc, lags] = xcov(spkbin1, spkbin2, xcovwindmax/binsize_spk);
    %     assert(~any(isnan(cc)))
    
    ccRealAll = [ccRealAll; cc'];
    
    
    %% ############################  SHIFT
    % PREDICTOR
    if (0) % ===== OLD VERSION, ONLY GET N VS. N+1 IN ONE DIRECTION.
        if (1)
            if t<ntrials
                tshift = t+1;
            else
                tshift = t-2;
            end
        else
            if t>1
                tshift = t-1;
            else
                tshift = t+1;
            end
        end
        spkbin2shift = dattmp2.SegmentsExtract(tshift).spk_Binned;
        x = dattmp2.SegmentsExtract(tshift).spk_Binned_x;
        
        % ============================ COMPUTE CORR
        % --- limit to data in "premotor" window
        indtmp = x>=windtmp(1) & x<windtmp(2);
        spkbin2shift = spkbin2shift(indtmp);
        
        if normspiketoprob ==1
            % ========= THEN MAKE IT PROB OF FINDINGS SPIKE IN BIN...
            spkbin2shift = spkbin2shift./sum(spkbin2shift);
        end
        
        % ------------------------------ calculate
        [cc, ~] = xcorr(spkbin1, spkbin2shift, xcovwindmax/binsize_spk, normmethod);
        %         if all(spkbin1==0) | all(spkbin2shift==0)
        %             cc = zeros(size(cc));
        %         end
        %         %                                 [cc, lags] = xcov(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
        %         assert(~any(isnan(cc)))
        
        
        %%  ==================== OUTPUT
        ccShiftAll = [ccShiftAll; cc'];
        
    else % ========== NEW METHOD, DO BOTH DIRECTIONS AND COMBINE (SO OUTPUT SAMPLE SIZE OF SHUFFLE TRIALS WILL BE (N-1)*2)
        %  ###################### DIRECTION 1
        if t>1
            tshift = t-1;
            spkbin2shift = dattmp2.SegmentsExtract(tshift).spk_Binned;
            x = dattmp2.SegmentsExtract(tshift).spk_Binned_x;
            
            % ============================ COMPUTE CORR
            % --- limit to data in "premotor" window
            indtmp = x>=windtmp(1) & x<windtmp(2);
            spkbin2shift = spkbin2shift(indtmp);
            
            if normspiketoprob ==1
                % ========= THEN MAKE IT PROB OF FINDINGS SPIKE IN BIN...
                spkbin2shift = spkbin2shift./sum(spkbin2shift);
            end
            
            % ------------------------------ calculate
            [cc, ~] = xcorr(spkbin1, spkbin2shift, xcovwindmax/binsize_spk, normmethod);
            %             if all(spkbin1==0) | all(spkbin2shift==0)
            %                 cc = zeros(size(cc));
            %             end
            %             %                                 [cc, lags] = xcov(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
            %             assert(~any(isnan(cc)))
            
            % ==================== OUTPUT
            ccShiftAll = [ccShiftAll; cc'];
        end
        
        
        %  ###################### DIRECTION 2
        if t<ntrials
            tshift = t+1;
            spkbin2shift = dattmp2.SegmentsExtract(tshift).spk_Binned;
            x = dattmp2.SegmentsExtract(tshift).spk_Binned_x;
            
            % ============================ COMPUTE CORR
            % --- limit to data in "premotor" window
            indtmp = x>=windtmp(1) & x<windtmp(2);
            spkbin2shift = spkbin2shift(indtmp);
            
            if normspiketoprob ==1
                % ========= THEN MAKE IT PROB OF FINDINGS SPIKE IN BIN...
                spkbin2shift = spkbin2shift./sum(spkbin2shift);
            end
            
            % ------------------------------ calculate
            [cc, ~] = xcorr(spkbin1, spkbin2shift, xcovwindmax/binsize_spk, normmethod);
            %             if all(spkbin1==0) | all(spkbin2shift==0)
            %                 cc = zeros(size(cc));
            %             end
            %             %                                 [cc, lags] = xcov(spkbin1, spkbin2shift, xcovwindmax/binsize_spk);
            %             assert(~any(isnan(cc)))
            
            % ==================== OUTPUT
            ccShiftAll = [ccShiftAll; cc'];
        end
        
        
    end
end


% ################################## OUTPUT
lags_sec = lags*binsize_spk;

%% ================ sanity check - plot CCs
if plotSummary==1 & any(strcmp(bregtmp, 'LMAN')) & any(strcmp(bregtmp, 'RA'))
    lt_figure; hold on;
    
    % --------- 1) firing rates
    lt_subplot(2,1,1); hold on;
    title(motifstr);
    ylabel([bregtmp{1} '(' num2str(neurons_thisset(nn)) ...
        ')-' bregtmp{2} '(' num2str(neurons_thisset(nnn)) ')']);
    
    y = [dattmp1.SegmentsExtract.FRsmooth_rate_CommonTrialDur];
    x = 0.001*(1:size(y,1));
    plot(x, y, '-', 'Color', [0.7 0.7 0.7]);
    plot(x, mean(y'), 'k', 'LineWidth', 2);
    ymax = max(y(:));
    
    y = [dattmp2.SegmentsExtract.FRsmooth_rate_CommonTrialDur];
    x = 0.001*(1:size(y,1));
    plot(x, y+1.2*ymax, '-', 'Color', [0.8 0.3 0.3]);
    plot(x, mean(y')+1.2*ymax, 'r', 'LineWidth', 2);
    
    axis tight
    
    % --- lines for syl onsets, offsets
    for kkk = 1:length(segextract_for_trialdur(1).motifsylOnsets)
        %                                     line([segextract_for_trialdur(1).motifsylOnsets(kkk) ...
        %                                         segextract_for_trialdur(1).motifsylOnsets(kkk)], ylim, 'Color', 'b');
        %                                     line([segextract_for_trialdur(1).motifsylOffsets(kkk) ...
        %                                         segextract_for_trialdur(1).motifsylOffsets(kkk)], ylim, 'Color', 'r');
        line([segextract_for_trialdur(1).motifsylOnsets(kkk) ...
            segextract_for_trialdur(1).motifsylOffsets(kkk)], [0 0], 'Color', 'b', ...
            'LineWidth', 2);
    end
    
    
    lt_subplot(2,2,3); hold on;
    title('r,dash=PSTH; b,dash=SHIFT');
    plot(lags*binsize_spk, ccRealAll', 'k-');
    plot(lags*binsize_spk, mean(ccRealAll,1), 'r', 'LineWidth', 2);
    plot(lags*binsize_spk, ccPSTH, 'r--', 'LineWidth', 2);
    plot(lags*binsize_spk, mean(ccShiftAll,1), 'b--', 'LineWidth', 2);
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    lt_subplot(2,2,4); hold on;
    plot(lags*binsize_spk, mean(ccRealAll,1) - ccPSTH', 'r-');
    plot(lags*binsize_spk, mean(ccRealAll,1) - mean(ccShiftAll,1), 'b-');
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end


if plotThisCC ==1
    
    % ---- line showing premotor window
    lt_subplot(8,2,1:14); hold on;
    axis tight;
    line([windtmp(1) windtmp(1)], ylim, 'Color', 'k');
    line([windtmp(2) windtmp(2)], ylim, 'Color', 'k');
    
    % ----
    lt_subplot(8, 2, 15); hold on;
    title('cross corr');
    
    shadedErrorBar(lags*0.001, mean(CCall_shift,1), lt_sem(CCall_shift), {'Color', [0.7 0.7 0.7]}, 1);
    shadedErrorBar(lags*0.001, mean(CCall,1), lt_sem(CCall), {'Color', 'k'}, 1);
    
    %                             lt_plot(lags*0.001, mean(CCall_shift,1), {'Errors', lt_sem(CCall_shift), 'Color', [0.7 0.7 0.7]});
    %                             lt_plot(lags*0.001, mean(CCall,1), {'Errors', lt_sem(CCall), 'Color', 'k'});
    
    
    pause
    close all;
end

