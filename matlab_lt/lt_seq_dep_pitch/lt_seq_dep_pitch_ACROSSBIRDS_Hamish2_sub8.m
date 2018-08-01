%% NOTE: 
% section 1: useful (plots raw contours for all syls)
% section 2 (PBS vs MUSC) and section 3 (PBS vs PBS) less useful, since
% subseumed by RAW plots below...


%% #################### PITCH CONTOUR SANITY CHECK [GOOD!!]
% OVERLAY RAW PC (PBS AND MUSC) AND TIME WINDOWS THAT WERE USED

N = length(DATSTRUCT.All_PitchCont_BASE_PBS);
nrend = 15;

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:N
    
    
    ndays = length(DATSTRUCT.All_PitchCont_BASE_PBS(i).All_PCmat);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([]);
    
    for dd = 1
        
        % ============== PBS
        PC = DATSTRUCT.All_PitchCont_BASE_PBS(i).All_PCmat{dd};
        t = DATSTRUCT.All_PitchCont_BASE_PBS(i).All_tbins{dd};
        twind = DATSTRUCT.All_PitchCont_BASE_PBS(i).All_twind(dd, :);
        
        if (1)
            if size(PC,1)>nrend
                PC = PC(randperm(size(PC,1), nrend),:);
            end
            plot(t, PC, '-k');
            
        else
            ymean = mean(PC,1);
            ystd = std(PC, [], 1);
            
            shadedErrorBar(t, ymean, ystd, {'Color', 'k'}, 1);
        end
        
        axis tight
        if ~any(isnan(twind))
            twind = t(twind);
            line([twind(1) twind(1)], ylim, 'Color', 'b');
            line([twind(2) twind(2)], ylim, 'Color', 'b');
        end
        
        % ============== musc
        PC = DATSTRUCT.All_PitchCont_BASE_MUSC(i).All_PCmat{dd};
        t = DATSTRUCT.All_PitchCont_BASE_MUSC(i).All_tbins{dd};
        twind = DATSTRUCT.All_PitchCont_BASE_MUSC(i).All_twind(dd, :);
        
        if (1)
            if size(PC,1)>nrend
                PC = PC(randperm(size(PC,1), nrend),:);
            end
            plot(t, PC, '-r');
            
        else
            ymean = mean(PC,1);
            ystd = std(PC, [], 1);
            
            shadedErrorBar(t, ymean, ystd, {'Color', 'r'}, 1);
        end
        
        if ~any(isnan(twind))
            twind = t(twind);
            line([twind(1) twind(1)], ylim, 'Color', 'b');
            line([twind(2) twind(2)], ylim, 'Color', 'b');
        end
        
        
        % -------------------
        if ~any(isnan(twind))
            xlim([0.016 twind(2)+0.02]);
        end
    end
end


%% #################### CONTOUR RELATIVE TO MUSCIMOL [PLOTS ALL]

% ================== BASELINE BIAS UP
Indstoplot = find(All_Istarg==1 & (All_FF_BASE_PBS - All_FF_BASE_MUSC)<0);
Nrends = 20; % renditions to take from edges
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(Indstoplot)
    
    indthis = Indstoplot(i);
    
    % ============================ BASELINE, PBS
    pcmat = All_PitchCont_BASE_PBS(indthis).All_PCmat{end};
    ffmat = All_PitchCont_BASE_PBS(indthis).All_ffvals{end};
    tthis = All_PitchCont_BASE_PBS(indthis).All_tbins{end};
    twind = All_PitchCont_BASE_PBS(indthis).All_twind(end,:);
    
    [~, indsort] = sort(ffmat);
    
    % ---- hi pitch
    pcthis_PBS_hi = pcmat(indsort(end-Nrends+1:end), :);
    % ---- lo pitch
    pcthis_PBS_lo = pcmat(indsort(1:Nrends), :);
    
    
    % ============================= BASELINE, MUSC
    pcmat = All_PitchCont_BASE_MUSC(indthis).All_PCmat{end};
    ffmat = All_PitchCont_BASE_MUSC(indthis).All_ffvals{end};
    
    % ---- all
    pcthis_MUSC = pcmat;
    
    
    % ######### PLOT, OVERLAY MEANS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('baseline');
    
    shadedErrorBar(tthis, mean(pcthis_PBS_lo,1), lt_sem(pcthis_PBS_lo), ...
        {'Color', [0.7 0.7 0.7]}, 1);
    shadedErrorBar(tthis, mean(pcthis_PBS_hi,1), lt_sem(pcthis_PBS_hi), ...
        {'Color', [0.7 0.7 0.7]}, 1);
    shadedErrorBar(tthis, mean(pcthis_MUSC,1), lt_sem(pcthis_MUSC), ...
        {'Color', 'r'}, 1);
    axis tight;
    
    % ######## SUBTRACT MEAN PBS
    
    % ######## SUBTRACT MEAN MUSC, PLOT MULTIPLE TRIALS
    hsplots = [];
    pcmean_MUSC = mean(pcthis_MUSC,1);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline (hi, minus MUSC)');
    plot(tthis, pcthis_PBS_hi - pcmean_MUSC, 'Color', 'k');
    lt_plot_zeroline;
    
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline (lo, minus MUSC)');
    plot(tthis, pcthis_PBS_lo - pcmean_MUSC, 'Color', 'k');
    lt_plot_zeroline;
    
    % ######## SUBTRACT MEAN MUSC (and subtract within trial mean), PLOT MULTIPLE TRIALS
    pcmean_MUSC = mean(pcthis_MUSC,1);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline (hi, minus MUSC, centered)');
    pcthis_norm = pcthis_PBS_hi - pcmean_MUSC;
    pcthis_norm = pcthis_norm - mean(pcthis_norm(:, twind(1):twind(2)), 2); % subtract within trial means (whtin window)
    plot(tthis, pcthis_norm, 'Color', 'k');
    lt_plot_zeroline;
    line([tthis(twind(1)) tthis(twind(1))], ylim);
    line([tthis(twind(2)) tthis(twind(2))], ylim);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline (low, minus MUSC, centered)');
    pcthis_norm = pcthis_PBS_lo - pcmean_MUSC;
    pcthis_norm = pcthis_norm - mean(pcthis_norm(:, twind(1):twind(2)), 2); % subtract within trial means (whtin window)
    plot(tthis, pcthis_norm, 'Color', 'k');
    lt_plot_zeroline;
    line([tthis(twind(1)) tthis(twind(1))], ylim);
    line([tthis(twind(2)) tthis(twind(2))], ylim);
    
    
    % ====================
    linkaxes(hsplots, 'xy');
end




%% ################## ALL CASES, COLLECT WIGGLE AND BASELINE BIAS
% [PLOTS CONTOUR RELATIVE TO MEAN DURING PBS]

% NOTE: calcualte baseline bias using pitch countours and update time
% windows

% SKIP THIS, IN TERMS OF COLLECTION IT DOES A SUBSET OF STUFF BELOPW, BUT
% NOT COMPLETELY WRITTEN. IN TERMS OF PLOTTING DOES USEFUL STUFF BUT IS
% REPLACED BY THE RAW PLOT FUNCTION BELOW
if (0)
    Indstoplot = find(All_Istarg==1);
    NrendsORIG = 30; % renditions to take from edges
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % =========== COLLECT THINGS
    WiggleAll = cell(length(Indstoplot), 2); % syls x [lo, hi]
    BaseBiasAll = nan(length(Indstoplot),1);
    
    for i=1:length(Indstoplot)
        
        indthis = Indstoplot(i);
        
        % ============================ BASELINE, PBS
        pcmat = All_PitchCont_BASE_PBS(indthis).All_PCmat{1};
        ffmat = All_PitchCont_BASE_PBS(indthis).All_ffvals{1};
        tthis = All_PitchCont_BASE_PBS(indthis).All_tbins{1};
        twind = All_PitchCont_BASE_PBS(indthis).All_twind(1,:);
        
        [~, indsort] = sort(ffmat);
        
        % =========================== minimum sample size
        if length(ffmat)<1*NrendsORIG
            continue
        end
        if length(ffmat)<2*NrendsORIG
            % then make nrends smaller
            Nrends = floor(length(ffmat)/2);
        else
            Nrends = NrendsORIG;
        end
        
        
        % ============================================
        % ---- hi pitch
        pcthis_PBS_hi = pcmat(indsort(end-Nrends+1:end), :);
        % ---- lo pitch
        pcthis_PBS_lo = pcmat(indsort(1:Nrends), :);
        
        
        
        % ######## NORMALIZE 1 - subtract x-trial mean
        % ---- all, mean
        pcmean_PBS_all = mean(pcmat, 1);
        pcmean_MUSC_all = mean(All_PitchCont_BASE_MUSC(indthis).All_PCmat{1},1);
        
        
        
        % ------------- HI
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title('baseline (hi, minus mean PBS)');
        
        plot(tthis, pcthis_PBS_hi - pcmean_PBS_all, '-k');
        plot(tthis, pcmean_MUSC_all - pcmean_PBS_all, '-r');
        lt_plot_zeroline;
        line([tthis(twind(1)) tthis(twind(1))], ylim);
        line([tthis(twind(2)) tthis(twind(2))], ylim);
        
        % ---------------- LO
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title('baseline (lo, minus mean PBS)');
        
        plot(tthis, pcthis_PBS_lo - pcmean_PBS_all, '-k');
        plot(tthis, pcmean_MUSC_all - pcmean_PBS_all, '-r');
        lt_plot_zeroline;
        line([tthis(twind(1)) tthis(twind(1))], ylim);
        line([tthis(twind(2)) tthis(twind(2))], ylim);
        
        
        
        % ######## NORMALIZE 1 - subtract x-trial mean and then within trial
        % mean
        % ----------- HI
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title('baseline (hi, normed)');
        
        pcthis_norm = pcthis_PBS_hi - pcmean_PBS_all; % subtract x-trail mean
        tmp = mean(pcthis_norm(:, twind(1):twind(end)), 2); % get within trial mean
        pcthis_norm = pcthis_norm - tmp; % subtract within trial mean
        
        plot(tthis, pcthis_norm, '-k');
        lt_plot_zeroline;
        line([tthis(twind(1)) tthis(twind(1))], ylim);
        line([tthis(twind(2)) tthis(twind(2))], ylim);
        
        % %%%%%%%%%% COLLECT
        % --- calcualte std of wiggle in time window
        wiggle = std(pcthis_norm(:, twind(1):twind(end)), [], 2);
        wiggle = mad(pcthis_norm(:, twind(1):twind(end)), [], 2);
        WiggleAll{i,2} = wiggle;
        
        % ----------- LO
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title('baseline (lo, normed)');
        
        pcthis_norm = pcthis_PBS_lo - pcmean_PBS_all; % subtract x-trail mean
        tmp = mean(pcthis_norm(:, twind(1):twind(end)), 2); % get within trial mean
        pcthis_norm = pcthis_norm - tmp; % subtract within trial mean
        
        plot(tthis, pcthis_norm, '-k');
        lt_plot_zeroline;
        line([tthis(twind(1)) tthis(twind(1))], ylim);
        line([tthis(twind(2)) tthis(twind(2))], ylim);
        % %%%%%%%%%% COLLECT
        % --- calcualte std of wiggle in time window
        wiggle = std(pcthis_norm(:, twind(1):twind(end)), [], 2);
        wiggle = mad(pcthis_norm(:, twind(1):twind(end)), [], 2);
        WiggleAll{i,1} = wiggle;
        
        
        
        % ######################### BASELINE AFP BIAS?
        pcmat_PBS = All_PitchCont_BASE_PBS(indthis).All_PCmat{1};
        pcmat_MUSC = All_PitchCont_BASE_MUSC(indthis).All_PCmat{1};
        
        ffPBS = mean(pcmat_PBS(:, twind(1):twind(end)),2);
        ffMUSC = mean(pcmat_MUSC(:, twind(1):twind(end)),2);
        afpbias = median(ffPBS) - median(ffMUSC);
        BaseBiasAll(i) = afpbias;
        
        %     if median(WiggleAll{i,2}) - median(WiggleAll{i,1})>12
        %         keyboard
        %     end
        % ====================
        linkaxes(hsplots, 'xy');
        
    end
    
    
    % ================= PLOT SUMMARY
    lt_figure; hold on;
    xlabel('base bias');
    ylabel('wiggle (HI) - wiggle (LO)');
    title('wiggle = median of std');
    
    indsgood = ~isnan(BaseBiasAll); % i.e. smalle sample size some
    
    wigglemat = cellfun(@median, WiggleAll, 'UniformOutput', 0);
    wigglemat = cell2mat(wigglemat(indsgood,:));
    biasmat = BaseBiasAll(indsgood);
    
    plot(biasmat, wigglemat(:,2) - wigglemat(:,1), 'ok')
    lt_regress(wigglemat(:,2) - wigglemat(:,1), biasmat, 1);
    lt_plot_zeroline;
end

