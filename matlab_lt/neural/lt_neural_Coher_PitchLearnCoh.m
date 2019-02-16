function lt_neural_Coher_PitchLearnCoh(OUTSTRUCT, PARAMS, SwitchCohStruct, ...
    SwitchStruct, nsegs)
%% lt 1/21/19 - plots learning trajectories. Both pitch and cohscalar


% nsegs = 2; % quartiles, then 4... (for both coh and learning)
disp('NOTE: assumes that, if multiple targets, then both train same direction... (not asserted in code)');
Nmin = 6; % min trials per bin
combineMultTarg = 1; % dedualt 1, takes average. (diff motifs for given switch)

%%
OUTSTRUCT = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct);

%%
%
% [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.motifID_unique});
% % [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.motifID_unique});
%
%
% % ----- get timecourse of learning (into quartiles)
% all_learndir = [];
% all_ffbinned = [];
% all_ffstd = [];
% all_cohbinned = [];
% all_cohstd = [];
% all_learnslope_byhour = []; % slope, [Ci]
% for i=1:length(indsgrpU)
%
%     % ===== Target
%     indsthis = find(indsgrp==indsgrpU(i) & OUTSTRUCT.istarg==1);
%     if ~any(indsthis)
%         continue
%     end
%
%     bnum = unique(OUTSTRUCT.bnum(indsthis));
%     enum = unique(OUTSTRUCT.enum(indsthis));
%     sw = unique(OUTSTRUCT.switch(indsthis));
%
%     assert(length(unique(OUTSTRUCT.motifID_unique(indsthis)))==1, 'do code for mult ta4rgets ...');
%     indsthis_allchan = indsthis;
%     indsthis = indsthis(1);
%
%     % --- things general across channels
%     ff = OUTSTRUCT.ffvals{indsthis};
%     tt = OUTSTRUCT.tvals{indsthis};
%     indsbase = OUTSTRUCT.indsbase{indsthis};
%     indsbase_epoch = OUTSTRUCT.indsbase_epoch{indsthis};
%     indsWN = OUTSTRUCT.indsWN{indsthis};
%     learndir = OUTSTRUCT.learndirTarg(indsthis);
%
%     % --- things to get across channels
%     cohscal = OUTSTRUCT.cohscal(indsthis_allchan);
%
%
%     if (0)
%        lt_figure; hold on;
%        plot(tt, ff, 'ok');
%     end
%
%
%     % ====================== GET WN/BASELINE MEASURES
%
%
%
%     % ====================== GET QUARTILES OF ELARNING AND COHSCAL
%     ntrialbin = sum(indsWN)/nsegs;
%
%     % --- get list of timeedges
%     tmp = find(indsWN);
%     trialedges = tmp(round(linspace(1, sum(indsWN), nsegs+1)));
%     trialedges(1) = trialedges(1)-1;
%
%     % --- add one bin for baseline
%     trialedges = [indsbase_epoch(1) trialedges];
%     trialedges(1) = trialedges(1)-1;
%
%     ff_binned = nan(1,nsegs+1);
%     ffstd_binned = nan(1,nsegs+1);
%     coh_binned = nan(1,nsegs+1);
%     cohstd_binned = nan(1,nsegs+1);
%     for ii=1:length(trialedges)-1
%        t1 = trialedges(ii)+1;
%        t2 = trialedges(ii+1);
%        time1 = tt(t1);
%        time2 = tt(t2);
%
%        % === collect stuff in this bin
%        ffthis = ff(t1:t2);
% %        ffstd_this = std(ff(t1:t2);
%        cohthis = cellfun(@(x)x(t1:t2), cohscal, 'UniformOutput', 0);
%        cohthis = mean(cell2mat(cohthis),1); % take mean acros chan pairs
%
%
% %        % === mean coh for other syls
% %        indstmp = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & ...
% %            OUTSTRUCT.switch==sw & OUTSTRUCT.istarg==0);
% %        for iii=1:length(indstmp)
% %           indstmp(iii)
% %        end
% %
%
%        % ============================ OUTPUT
%        ff_binned(ii) = mean(ffthis);
%        ffstd_binned(ii) = std(ffthis);
%        coh_binned(ii) = mean(cohthis);
%        cohstd_binned(ii) = std(cohthis);
%     end
%
%
%     % ========================= GET SLOPE OVER WN
%     ftmp = ff(indsWN);
%     ttmp = tt(indsWN);
%     ttmp = (ttmp-ttmp(1))*24; % convert to hours.
%     [~, ~, ~, ~, ~, stats] = lt_regress(ftmp, ttmp, 0, 0, 0, 0, '');
%     all_learnslope_byhour = [all_learnslope_byhour; [stats.slope stats.slopeCI]];
%
%
%     % ################## OUTPUT, ONE PER EXPERIMENT
%     all_learndir = [all_learndir; learndir];
%     all_ffbinned = [all_ffbinned; ff_binned];
%     all_ffstd = [all_ffstd; ffstd_binned];
%     all_cohbinned = [all_cohbinned; coh_binned];
%     all_cohstd = [all_cohstd; cohstd_binned];
% end
%
%
%
%
% %% ================== nromalize learning
% % == learning, hz minus base, in learn dir
% all_fflearnrelbase = all_ffbinned;
% all_fflearnrelbase = (all_fflearnrelbase-all_fflearnrelbase(:,1)).*all_learndir;
%
%
% % == learning, zscore minsu base, learn dir
% all_learnz_learndir = (all_ffbinned - all_ffbinned(:,1))./(all_ffstd(:,1));
% all_learnz_learndir = all_learnz_learndir.*all_learndir;
%
% %% ========================================== PLOT
% figcount=1;
% subplotrows=3;
% subplotcols=2;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
%
%
% % =========== 1) LEARNING TRAJECTORY
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('binned (1=bnase, 2:end=WN[by ind])');
% ylabel('ff');
%
% x=1:size(all_ffbinned,2);
% y = all_fflearnrelbase;
% plot(x, y', '-ok');
%
% ymean = mean(y,1);
% ysem = lt_sem(y);
% lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
% lt_plot_zeroline;
%
%
% % --
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('binned (1=bnase, 2:end=WN[by ind])');
% ylabel('learn(z)');
%
% x=1:size(all_ffbinned,2);
% y = all_learnz_learndir;
% plot(x, y', '-ok');
%
% ymean = mean(y,1);
% ysem = lt_sem(y);
% lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
% lt_plot_zeroline;
%
%
% % ====== slope during elarning
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('slope dur elarn (dir of learn)');
%
% y = all_learnslope_byhour(:,1).*all_learndir;
% lt_plot_histogram(y, [], 1, 0, []);
% lt_plot_zeroline_vert;
%
%
% % ===== COHERENCE INCREASE DURING LEARNING
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('binned (1=bnase, 2:end=WN[by ind])');
% ylabel('coh');
%
% x=1:size(all_ffbinned,2);
% y = all_cohbinned;
% plot(x, y', '-ok');
%
% ymean = mean(y,1);
% ysem = lt_sem(y);
% lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
% lt_plot_zeroline;
%
%
% % ===== COHERENCE INCREASE DURING LEARNING
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('binned (1=bnase, 2:end=WN[by ind])');
% ylabel('coh');
%
% x=1:size(all_ffbinned,2);
% y = all_cohbinned-all_cohbinned(:,1);
% plot(x, y', '-ok');
%
% ymean = mean(y,1);
% ysem = lt_sem(y);
% lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
% lt_plot_zeroline;
%
%
%
% % ============= OCHERENCE VERSUS LEARNING
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('learning (zscore)');
% ylabel('coherence change (coh)');
%
% x = all_learnz_learndir(:,end);
% y = all_cohbinned(:,end) - all_cohbinned(:,1);
% plot(x,y, 'ok');
%
%
% % ============= OCHERENCE VERSUS LEARNING
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('learning (zscore)');
% ylabel('coherence change (z)');
%
% x = all_learnz_learndir(:,end);
% y = (all_cohbinned(:,end) - all_cohbinned(:,1))./(all_cohstd(:,1));
% plot(x,y, 'ok');
%




%% NEW VERSION, EXTRACTS FOR ALL CASES, THEN ASKS WHAT IS TARGS...

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.motifID_unique});
% [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.motifID_unique});


% ----- get timecourse of learning (into quartiles)
all_learndir = [];
all_ffbinned = [];
all_ffstd = [];
all_cohbinned = [];
all_cohstd = [];
all_learnslope_byhour = []; % slope, [Ci]
all_istarg = [];
all_issame = [];
all_bnum = [];
all_enum =[];
all_swnum =[];
for i=1:length(indsgrpU)
    
    % ===== Target
    indsthis = find(indsgrp==indsgrpU(i));
    
    bnum = unique(OUTSTRUCT.bnum(indsthis));
    enum = unique(OUTSTRUCT.enum(indsthis));
    sw = unique(OUTSTRUCT.switch(indsthis));
    istarg = unique(OUTSTRUCT.istarg(indsthis));
    issame = unique(OUTSTRUCT.issame(indsthis));
    
    %     assert(length(unique(OUTSTRUCT.motifID_unique(indsthis)))==1, 'do code for mult ta4rgets ...');
    indsthis_allchan = indsthis;
    indsthis = indsthis(1);
    
    % --- things general across channels
    ff = OUTSTRUCT.ffvals{indsthis};
    tt = OUTSTRUCT.tvals{indsthis};
    indsbase = OUTSTRUCT.indsbase{indsthis};
    indsbase_epoch = OUTSTRUCT.indsbase_epoch{indsthis};
    indsWN = OUTSTRUCT.indsWN{indsthis};
    learndir = OUTSTRUCT.learndirTarg(indsthis);
    
    % --- things to get across channels
    cohscal = OUTSTRUCT.cohscal(indsthis_allchan);
    if any(isnan(cohscal{1}))
        disp('REMOIVNG EXPT - COH SCALARS ARE NAN!!');
        pause;
        continue
    end
    if sum(indsWN)<Nmin*nsegs | length(indsbase_epoch)<Nmin*nsegs
        continue
    end
    
    
    if (0)
        lt_figure; hold on;
        plot(tt, ff, 'ok');
    end
    
    
    % ====================== GET WN/BASELINE MEASURES
    
    
    
    % ====================== GET QUARTILES OF ELARNING AND COHSCAL
    % --- get list of timeedges
    tmp = find(indsWN);
    trialedges = tmp(round(linspace(1, sum(indsWN), nsegs+1)));
    trialedges(1) = trialedges(1)-1;
    
    % --- add one bin for baseline
    trialedges = [indsbase_epoch(1) trialedges];
    trialedges(1) = trialedges(1)-1;
    
    ff_binned = nan(1,nsegs+1);
    ffstd_binned = nan(1,nsegs+1);
    coh_binned = nan(1,nsegs+1);
    cohstd_binned = nan(1,nsegs+1);
    for ii=1:length(trialedges)-1
        t1 = trialedges(ii)+1;
        t2 = trialedges(ii+1);
        time1 = tt(t1);
        time2 = tt(t2);
        
        % === collect stuff in this bin
        ffthis = ff(t1:t2);
        %        ffstd_this = std(ff(t1:t2);
        cohthis = cellfun(@(x)x(t1:t2), cohscal, 'UniformOutput', 0);
        cohthis = mean(cell2mat(cohthis),1); % take mean acros chan pairs
        
        
        %        % === mean coh for other syls
        %        indstmp = find(OUTSTRUCT.bnum==bnum & OUTSTRUCT.enum==enum & ...
        %            OUTSTRUCT.switch==sw & OUTSTRUCT.istarg==0);
        %        for iii=1:length(indstmp)
        %           indstmp(iii)
        %        end
        %
        if isnan(mean(cohthis))
            keyboard
        end
        
        % ============================ OUTPUT
        ff_binned(ii) = mean(ffthis);
        ffstd_binned(ii) = std(ffthis);
        coh_binned(ii) = mean(cohthis);
        cohstd_binned(ii) = std(cohthis);
    end
    
    
    % ========================= GET SLOPE OVER WN
    if all(isnan(ff))
        all_learnslope_byhour = [all_learnslope_byhour; [nan nan nan]];
    else
        ftmp = ff(indsWN);
        ttmp = tt(indsWN);
        ttmp = (ttmp-ttmp(1))*24; % convert to hours.
        [~, ~, ~, ~, ~, stats] = lt_regress(ftmp, ttmp, 0, 0, 0, 0, '');
        all_learnslope_byhour = [all_learnslope_byhour; [stats.slope stats.slopeCI]];
    end
    
    
    % ################## OUTPUT
    all_learndir = [all_learndir; learndir];
    all_ffbinned = [all_ffbinned; ff_binned];
    all_ffstd = [all_ffstd; ffstd_binned];
    all_cohbinned = [all_cohbinned; coh_binned];
    all_cohstd = [all_cohstd; cohstd_binned];
    all_istarg = [all_istarg; istarg];
    all_issame = [all_issame; issame];
    
    all_bnum = [all_bnum; bnum];
    all_enum =[all_enum; enum];
    all_swnum =[all_swnum; sw];
    
end

%% % === also, if multiple targets, then combines them (average)


if combineMultTarg==1
    [indsgrp, indsgrpU] = lt_tools_grp2idx({all_bnum, all_enum, all_swnum});
    indstoremove = [];
    all_istarg_old = all_istarg;
    for i=1:length(indsgrpU)
        
        % ============ TARG
        indstarg = find(indsgrp==indsgrpU(i) & all_istarg_old==1);
        
        % ========================== COMBINE IF MULT TARGETS
        if length(indstarg)>1
            indstoremove = [indstoremove; indstarg]; % rmeomve later
            
            % append to end
            %             all_learnslope_byhour
            all_learnslope_byhour = [all_learnslope_byhour; mean(all_learnslope_byhour(indstarg,:),1)];
            
            all_learndir = [all_learndir; mean(all_learndir(indstarg,:),1)];
            all_ffbinned = [all_ffbinned; mean(all_ffbinned(indstarg,:),1)];
            all_ffstd = [all_ffstd; mean(all_ffstd(indstarg,:),1)];
            all_cohbinned = [all_cohbinned; mean(all_cohbinned(indstarg,:),1)];
            all_cohstd = [all_cohstd; mean(all_cohstd(indstarg,:),1)];
            all_istarg = [all_istarg; mean(all_istarg(indstarg,:),1)];
            all_issame = [all_issame; mean(all_issame(indstarg,:),1)];
            all_bnum = [all_bnum; mean(all_bnum(indstarg,:),1)];
            all_enum = [all_enum; mean(all_enum(indstarg,:),1)];
            all_swnum = [all_swnum; mean(all_swnum(indstarg,:),1)];
        end
    end
    all_learnslope_byhour(indstoremove,:) = [];
    all_learndir(indstoremove,:) = [];
    all_ffbinned(indstoremove,:) = [];
    all_ffstd(indstoremove,:) = [];
    all_cohbinned(indstoremove,:) = [];
    all_cohstd(indstoremove,:) = [];
    all_istarg(indstoremove,:) = [];
    all_bnum(indstoremove,:) = [];
    all_enum(indstoremove,:) = [];
    all_swnum(indstoremove,:) = [];
    all_issame(indstoremove,:) = [];
end
%% ================= FOR EACH TARGET, GET STATS RELATIVE TO ALL OTHER SYLS

[indsgrp, indsgrpU] = lt_tools_grp2idx({all_bnum, all_enum, all_swnum});
all_cohbinned_minusnontarg = nan(size(all_cohbinned));
all_cohbinned_zrelnontarg = nan(size(all_cohbinned));
for i=1:length(indsgrpU)
    
    % ============ TARG
    indstarg = indsgrp==indsgrpU(i) & all_istarg==1;
    indsnontarg = indsgrp==indsgrpU(i) & all_istarg==0;
    
    % ============ get for targ as difference from mean of nontargs
    cohmean_nontarg = mean(all_cohbinned(indsnontarg, :),1);
    cohstd_nontarg = std(all_cohbinned(indsnontarg, :),[], 1);
    
    all_cohbinned_minusnontarg(indstarg, :) = ...
        all_cohbinned(indstarg,:) - cohmean_nontarg;
    
    all_cohbinned_zrelnontarg(indstarg,:) = ...
        (all_cohbinned(indstarg,:) - cohmean_nontarg)./cohstd_nontarg;
end


%% === warning
disp('NOTE: z relative to nontarg has outliers becuase soem cases not many nontargs exist..');
pause;
%% ================== nromalize learning
% == learning, hz minus base, in learn dir
all_fflearnrelbase = all_ffbinned;
all_fflearnrelbase = (all_fflearnrelbase-all_fflearnrelbase(:,1)).*all_learndir;


% == learning, zscore minsu base, learn dir
all_learnz_learndir = (all_ffbinned - all_ffbinned(:,1))./(all_ffstd(:,1));
all_learnz_learndir = all_learnz_learndir.*all_learndir;



%% ============== one for each expt
[indsgrp, indsgrpU] = lt_tools_grp2idx({all_bnum, all_enum, all_swnum});

all_learning_zscore_targdir = nan(length(indsgrpU), size(all_ffbinned,2), 3); % (n,nbins, [t,s,d])

for i=1:length(indsgrpU)
    
    % targ
    indsthis = indsgrp==indsgrpU(i) & all_istarg==1;
    all_learning_zscore_targdir(i, :, 1) = nanmean(all_learnz_learndir(indsthis,:),1);
    
    % same
    indsthis = indsgrp==indsgrpU(i) & all_istarg==0 & all_issame==1;
    all_learning_zscore_targdir(i, :, 2) = nanmean(all_learnz_learndir(indsthis,:),1);
    
    % diff
    indsthis = indsgrp==indsgrpU(i) & all_istarg==0 & all_issame==0;
    all_learning_zscore_targdir(i, :, 3) = nanmean(all_learnz_learndir(indsthis,:),1);
end



%% ========================================== PLOT
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% =========== 1) LEARNING TRAJECTORY
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('ff');

x=1:size(all_ffbinned,2);
y = all_fflearnrelbase(all_istarg==1,:);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% --ZSCORE (TARGET)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('learn(z)');

x=1:size(all_ffbinned,2);
y = all_learnz_learndir(all_istarg==1,:);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;



% ========= ZSCORE (SAME TPE0
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME');
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('learn(z)');

x=1:size(all_ffbinned,2);
y = squeeze(all_learning_zscore_targdir(:,:,2));
plot(x, y', '-ok');

ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ========= ZSCORE (SAME TPE0
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF');
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('learn(z)');

x=1:size(all_ffbinned,2);
y = squeeze(all_learning_zscore_targdir(:,:,3));
plot(x, y', '-ok');

ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ===========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME (each syl)');
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('learn(z)');

x=1:size(all_ffbinned,2);
y = all_learnz_learndir(all_istarg==0 & all_issame==1,:);
plot(x, y', '-ok');

ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;

% ======
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF (each syl)');
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('learn(z)');

x=1:size(all_ffbinned,2);
y = all_learnz_learndir(all_istarg==0 & all_issame==0,:);
plot(x, y', '-ok');

ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ======
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF (each syl)');
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('learn(ff minus base)');

x=1:size(all_ffbinned,2);
y = all_fflearnrelbase(all_istarg==0 & all_issame==0,:);
plot(x, y', '-ok');

ymean = nanmean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ====== slope during elarning
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('slope dur elarn (dir of learn)');

y = all_learnslope_byhour(all_istarg==1,1).*all_learndir(all_istarg==1);
lt_plot_histogram(y, [], 1, 0, []);
lt_plot_zeroline_vert;


% ===== COHERENCE INCREASE DURING LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh');

x=1:size(all_ffbinned,2);
y = all_cohbinned(all_istarg==1,:);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ===== COHERENCE INCREASE DURING LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh');

x=1:size(all_ffbinned,2);
y = all_cohbinned-all_cohbinned(:,1);
y = y(all_istarg==1, :);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ===== COHERENCE INCREASE DURING LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh (minus nontarg)');

x=1:size(all_ffbinned,2);
y = all_cohbinned_minusnontarg(all_istarg==1,:);
y = y-y(:,1);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ===== COHERENCE INCREASE DURING LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh (minus z rel nontarg (in each bin)');

x=1:size(all_ffbinned,2);
y = all_cohbinned_zrelnontarg(all_istarg==1,:);
% y = y-y(:,1);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ===== COHERENCE INCREASE DURING LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh (z rel own base)');

x=1:size(all_ffbinned,2);
y = (all_cohbinned - all_cohbinned(:,1))./all_cohstd(:,1);
y = y(all_istarg==1,:);
% y = y-y(:,1);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ===== COHERENCE INCREASE DURING LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh (minus z rel nontarg (in each bin)');

x=1:size(all_ffbinned,2);
y = all_cohbinned_zrelnontarg(all_istarg==1,:);
y = y-y(:,1);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ========== TIMECOURSE OF COH AND LEARN CHANG
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
% ylabel('coh(zrelbase) minus learn(zrelbase)');
ylabel('learn(zrelbase) minus coh(zrelbase)');
y2 = (all_cohbinned - all_cohbinned(:,1))./all_cohstd(:,1);
y2 = y2(all_istarg==1,:);

y1 = all_learnz_learndir(all_istarg==1,:);

y = y1-y2;
x=1:size(all_ffbinned,2);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;


% ========== TIMECOURSE OF COH AND LEARN CHANG
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binned (1=bnase, 2:end=WN[by ind])');
ylabel('coh(minusnontarg) minus learn(zrelbase)');
y1 = all_cohbinned_minusnontarg-all_cohbinned_minusnontarg(:,1);
y1 = y1(all_istarg==1,:);

y2 = all_learnz_learndir(all_istarg==1,:);

y = y1-y2;
x=1:size(all_ffbinned,2);
plot(x, y', '-ok');

ymean = mean(y,1);
ysem = lt_sem(y);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
lt_plot_zeroline;



% % ============= OCHERENCE VERSUS LEARNING
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('learning (zscore)');
% ylabel('coherence change (coh)');
%
% x = all_learnz_learndir(all_istarg==1,end);
% y = all_cohbinned(:,end) - all_cohbinned(:,1);
% y = y(all_istarg==1,:);
% plot(x,y, 'ok');


% ============= OCHERENCE VERSUS LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (zscore)');
ylabel('coherence change (z)');

x = all_learnz_learndir(all_istarg==1,end);
y = (all_cohbinned(:,end) - all_cohbinned(:,1))./(all_cohstd(:,1));
y = y(all_istarg==1,:);
plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1, 'k');



% ============= OCHERENCE VERSUS LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (zscore)');
ylabel('coherence change (coh)');

x = all_learnz_learndir(all_istarg==1,end);
y = (all_cohbinned(:,end) - all_cohbinned(:,1));
y = y(all_istarg==1,:);
plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1, 'k');


% ============= OCHERENCE VERSUS LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (zscore)');
ylabel('coherence change (minus nontarg, change from base)');

x = all_learnz_learndir(all_istarg==1,end);
y = all_cohbinned_minusnontarg(all_istarg==1,:);
y = y(:,end)-y(:,1);
% y = y(:,end);
plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1, 'k');


% ============= OCHERENCE VERSUS LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('learning (zscore)');
ylabel('coherence (z rel nontarg)');

x = all_learnz_learndir(all_istarg==1,end);
y = all_cohbinned_zrelnontarg(all_istarg==1,:);
y = y(:,end)-y(:,1);
% y = y(:,end);
plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1, 'k');

