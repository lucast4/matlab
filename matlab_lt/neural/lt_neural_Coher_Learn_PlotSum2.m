function lt_neural_Coher_Learn_PlotSum2(OUTSTRUCT, PARAMS, SwitchStruct, ...
    sumplottype, plotAllSwitchRaw, clim, fieldtoplot, birdstoplot, ...
    timewindowtoplot, zscoreLFP, expttoplot, swtoplot, ffbinsedges, ...
    indtoget_b_e_s, useAbsVal)
%% lt 11/5/18 - simpler version that is more general
% CAn do any of the matrix data (e.g. S, coh...)

% useAbsVal = 1; then for coherence difference uses abs value of diff
% [takes the abs after doing teh averages across channels];

% sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs


if ~exist('useAbsVal', 'var')
    useAbsVal = 0;
end

%%
if (0)
    assert(all(strcmp(OUTSTRUCT.bregionpair, 'LMAN-RA')), 'assumes chan1 is LAMN, chang 2 is RA...');
end

% clim = [-0.15 0.15];
plotEachTrial = 0; % for the timecourses, if 0 then only plots means, sem.

if ~exist('indtoget_b_e_s', 'var')
    indtoget_b_e_s = [];
end

%% onluy plot certain birds?

OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT); % so followniog code works...

% === filter each dimension independenyl
if ~isempty(birdstoplot)
    indstokeep = ismember(OUTSTRUCT.bnum, birdstoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

if ~isempty(expttoplot)
    indstokeep = ismember(OUTSTRUCT.enum, expttoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

if ~isempty(swtoplot)
    indstokeep = ismember(OUTSTRUCT.switch, swtoplot);
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end

% == filter by specific types of switches.
if ~isempty(indtoget_b_e_s)
    indstokeep = ismember([OUTSTRUCT.bnum OUTSTRUCT.enum OUTSTRUCT.switch], indtoget_b_e_s, 'rows');
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);
end
%% ======= field type to plot
DATSTRUCT = struct;

if strcmp(fieldtoplot, 'coher')
    %% for coherence [diff]
    
    
    if strcmp(sumplottype, 'switches')
        % ============
        %     fieldtoget = 'CohMean_WNminusBase';
        %     [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
        %         lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        
        fieldtoget = 'CohMean_WN';
        [~, ~, ~, ~, allbnum, allenum, allswnum, allDat2] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        
        fieldtoget = 'CohMean_Base';
        [~, ~, ~, ~, allbnum, allenum, allswnum, allDat1] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        
        % ---- subtract WN from base
        cohdiff = nan(size(allDat1));
        for i=1:size(cohdiff,3)
            for ii=1:size(cohdiff,4)
                
                cohdiff(:,:,i,ii) = allDat2(:,:,i, ii) - allDat1(:,:,i, ii);
            end
        end
        
        % save dat
        DATSTRUCT.cohdiff = cohdiff;
        
    elseif strcmp(sumplottype, 'chanpairs')
        % ============
        %     fieldtoget = 'CohMean_WNminusBase';
        %     [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
        %         lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        
        fieldtoget = 'CohMean_WN';
        [allbnum, allenum, allswnum, allDat2] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        
        fieldtoget = 'CohMean_Base';
        [allbnum, allenum, allswnum, allDat1] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        
        % ---- subtract WN from base
        cohdiff = nan(size(allDat1));
        for i=1:size(cohdiff,3)
            for ii=1:size(cohdiff,4)
                
                cohdiff(:,:,i,ii) = allDat2(:,:,i, ii) - allDat1(:,:,i, ii);
            end
        end
        
        % save dat
        DATSTRUCT.cohdiff = cohdiff;
        
    end
    
    if useAbsVal==1
        DATSTRUCT.cohdiff = abs(DATSTRUCT.cohdiff);
    end
    
    
    
elseif strcmp(fieldtoplot, 'spec')
    
    %% EXTRACT
    
    if strcmp(sumplottype, 'switches')
        % ============ 1) chan1, base
        fieldtoget = 'Spec1Mean_Base';
        [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan1.base = allDat;
        %         prctile(allDat(:), [2.75 97.5])
        
        % ============ 1) chan1, wn
        fieldtoget = 'Spec1Mean_WN';
        [~, ~, ~, ~, allbnum2, allenum, allswnum2, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan1.wn= allDat;
        %         prctile(allDat(:), [2.75 97.5])
        
        % ============ 1) chan2, base
        fieldtoget = 'Spec2Mean_Base';
        [~, ~, ~, ~, allbnum3, allenum, allswnum3, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan2.base = allDat;
        %         prctile(allDat(:), [2.75 97.5])
        
        % ============ 1) chan2, WN
        fieldtoget = 'Spec2Mean_WN';
        [~, ~, ~, ~, allbnum4, allenum, allswnum4, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan2.wn = allDat;
        %         prctile(allDat(:), [2.75 97.5])
        
        % --- sanity check
        assert(all(all(diff([allswnum allswnum2 allswnum3 allswnum4], 1, 2)==0)), 'extractions nt identical..');
        
    elseif strcmp(sumplottype, 'chanpairs')
        % ============ 1) chan1, base
        fieldtoget = 'Spec1Mean_Base';
        [allbnum1, allenum, allswnum, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan1.base = allDat;
        
        % ============ 1) chan1, wn
        fieldtoget = 'Spec1Mean_WN';
        [allbnum1, allenum, allswnum2, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan1.wn= allDat;
        
        % ============ 1) chan2, base
        fieldtoget = 'Spec2Mean_Base';
        [allbnum1, allenum, allswnum3, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan2.base = allDat;
        
        % ============ 1) chan2, WN
        fieldtoget = 'Spec2Mean_WN';
        [allbnum1, allenum, allswnum4, allDat] = ...
            lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
        % save dat
        DATSTRUCT.chan2.wn = allDat;
        
        % --- sanity check
        assert(all(all(diff([allswnum allswnum2 allswnum3 allswnum4], 1, 2)==0)), 'extractions nt identical..');
    end
end


%%
if (0) % THIS IS TO INSERT COHEREHNCE DATA (WN AND BASE) INTO LFP DATA, TO SEE IF CALCULATING
    % COHERENCE DIFFERENCE USING WN AND BASE LEADS TO SAME RESULT AS
    % PREVIOUS ANALYSES CALCULATING WN MINUS BASE BEFORE. INDEED IT DOES.
    % ============ 1) chan1, base
    fieldtoget = 'CohMean_Base';
    [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
    % save dat
    DATSTRUCT.chan1.base = allDat;
    
    % ============ 1) chan1, wn
    fieldtoget = 'CohMean_WN';
    [~, ~, ~, ~, allbnum1, allenum, allswnum2, allDat] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
    % save dat
    DATSTRUCT.chan1.wn= allDat;
    
    % ============ 1) chan2, base
    fieldtoget = 'CohMean_Base';
    [~, ~, ~, ~, allbnum1, allenum, allswnum3, allDat] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
    % save dat
    DATSTRUCT.chan2.base = allDat;
    
    % ============ 1) chan2, WN
    fieldtoget = 'CohMean_WN';
    [~, ~, ~, ~, allbnum1, allenum, allswnum4, allDat] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);
    % save dat
    DATSTRUCT.chan2.wn = allDat;
    
    % --- sanity check
    assert(all(all(diff([allswnum allswnum2 allswnum3 allswnum4], 1, 2)==0)), 'extractions nt identical..');
end


%% ==========

tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

XLIM = [-0.1 -0.02];
YLIM = [20 120];




%% PLOT - ONE PLOT PER SYL TYPE, AVERAGING EITHER OVER SWITCH OR CHANNEL PAIR

if strcmp(fieldtoplot, 'coher')
    
    CohMatDat = DATSTRUCT.cohdiff;
    %     CohMatDat_base = allswitch_CohBaseMat;
    %     CohMatDat_WN = allswitch_CohWNMat;
    CohScalDat = nan;
    birddat = allbnum;
    
    indstmp = ffbins>5 & ffbins<150;
    CohMatDat = CohMatDat(:,indstmp, :, :);
    ffbins = ffbins(indstmp);
    
    %% === check (diff and [wn - base] are the same?)
    if (0)
        tmp1 = CohMatDat_WN(:) - CohMatDat_base(:);
        tmp2 = CohMatDat(:);
        
        figure; hold on;
        plot(tmp1(~isnan(tmp1)), tmp2(~isnan(tmp2)), 'ob')
        figure; hold on;
        lt_plot_histogram(tmp2(~isnan(tmp2))-tmp1(~isnan(tmp1)));
    end
    
    %% ============= DIFFERENCE MATRICES
    figcount=1;
    subplotrows=3;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ------------- TARG
    cohmat = squeeze(CohMatDat(:,:,1,:));
    plottit = 'TARG';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title(plottit);
    
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, 0, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    
    % ------------- SAME
    cohmat = squeeze(CohMatDat(:,:,2,:));
    plottit = 'SAME';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title(plottit);
    
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, 0, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    
    % ------------- diff
    cohmat = squeeze(CohMatDat(:,:,3,:));
    plottit = 'DIFF';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title(plottit);
    colorbar('East');
    % 2. ffband differences
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, 0, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    
    
    % ================== TARG MINUS SAME
    cohmat1 = squeeze(CohMatDat(:,:,1,:));
    cohmat2 = squeeze(CohMatDat(:,:,2,:));
    cohmat = cohmat1 - cohmat2;
    plottit = 'TARG - SAME';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title(plottit);
    % 2. ffband differences
    plotindivtraces =1;
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    % 2. ffband differences
    plotindivtraces =0;
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    
    
    % ================== TARG MINUS DIFF
    cohmat1 = squeeze(CohMatDat(:,:,1,:));
    cohmat2 = squeeze(CohMatDat(:,:,3,:));
    cohmat = cohmat1 - cohmat2;
    plottit = 'TARG - DIFF';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
    title(plottit);
    % 2. ffband differences
    plotindivtraces =1;
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    % 2. ffband differences
    plotindivtraces =0;
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces, ffbinsedges);
    lt_plot_zeroline;
    lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
    
    
    % ########################### PVALUE PLOTS
    figcount=1;
    subplotrows=2;
    subplotcols=5;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ------------- TARG
    cohmat = squeeze(CohMatDat(:,:,1,:));
    plottit = 'TARG';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 3, '');
    title(plottit);
    
    % ------------- SAME
    cohmat = squeeze(CohMatDat(:,:,2,:));
    plottit = 'SAME';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 3, '');
    title(plottit);
    
    
    % ------------- diff
    cohmat = squeeze(CohMatDat(:,:,3,:));
    plottit = 'DIFF';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 3, '');
    title(plottit);
    colorbar('East');
    
    
    % ================== TARG MINUS SAME
    cohmat1 = squeeze(CohMatDat(:,:,1,:));
    cohmat2 = squeeze(CohMatDat(:,:,2,:));
    cohmat = cohmat1 - cohmat2;
    plottit = 'TARG - SAME';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 3, '');
    title(plottit);
    
    
    % ================== TARG MINUS DIFF
    cohmat1 = squeeze(CohMatDat(:,:,1,:));
    cohmat2 = squeeze(CohMatDat(:,:,3,:));
    cohmat = cohmat1 - cohmat2;
    plottit = 'TARG - DIFF';
    % 1. cohgram of diff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(cohmat, tbins, ffbins, 3, '');
    title(plottit);
    
    
    % ########################### SCALAR PLOTS
    if any(~isnan(CohScalDat(:)))
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        xlabel('TARG - SAME - DIFF');
        ylabel('change in coh');
        title('all dat (including some expt no same');
        x = [1 2 3];
        y = CohScalDat;
        plot(x,y', '-o', 'Color', [0.7 0.7 0.7]);
        lt_plot(x+0.15, nanmean(y), {'Errors', lt_sem(y)});
        xlim([0 4]);
        lt_plot_zeroline;
        
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        xlabel('TARG - SAME - DIFF');
        ylabel('change in coh');
        title('(only if have all syl types)');
        y = CohScalDat(~any(isnan(CohScalDat)'), :);
        x = [1 2 3];
        plot(x,y', '-o', 'Color', [0.7 0.7 0.7]);
        lt_plot(x+0.15, nanmean(y), {'Errors', lt_sem(y)});
        xlim([0 4]);
        lt_plot_zeroline;
        
        % --- test vs. nontarg
        [~, p] = ttest(y(:,1), y(:,2));
        lt_plot_pvalue(p, 'targVSsame', 1);
        
        [~, p] = ttest(y(:,1), y(:,3));
        lt_plot_pvalue(p, 'targVSdiff', 2);
        
        % [~, p] = ttest(y(:,2), y(:,3));
        % lt_plot_pvalue(p, 'sameVSdiff', 2);
        
        
        % ########################### SCALAR PLOTS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        xlabel('TARG - DIFF');
        ylabel('change in coh');
        title('all dat (ignore whether has same type)');
        x = [1 2];
        y = CohScalDat(:, [1 3]);
        plot(x,y', '-o', 'Color', [0.7 0.7 0.7]);
        lt_plot(x+0.15, nanmean(y), {'Errors', lt_sem(y)});
        xlim([0 4]);
        lt_plot_zeroline;
        
        % --- test vs. nontarg
        [~, p] = ttest(y(:,1), y(:,2));
        lt_plot_pvalue(p, 'targVSdiff', 1);
        
        
        % ########################### SCALAR PLOTS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        xlabel('TARG - DIFF');
        ylabel('change in coh');
        title('[same, but separate by bird]');
        x = [1 2];
        pcols = lt_make_plot_colors(length(unique(birddat)), 0, 0);
        for j=unique(birddat)'
            indstmp = birddat==j;
            y = CohScalDat(indstmp, [1 3]);
            plot(x,y', '-', 'Color', pcols{j});
            lt_plot(x+0.2, nanmean(y), {'Errors', lt_sem(y), 'Color', pcols{j}});
            % --- test vs. nontarg
            [~, p] = ttest(y(:,1), y(:,2));
            lt_plot_text(2.5, 0+0.01*j, ['b#' num2str(j) 'p=' num2str(p)], pcols{j});
        end
        xlim([0 4]);
        lt_plot_zeroline;
    end
    
end

if strcmp(fieldtoplot, 'spec')
    %% ========== PLOTTING LFP
    figcount=1;
    subplotrows=6;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    if min(DATSTRUCT.chan1.base(:))>0 & max(DATSTRUCT.chan1.base(:))>0.5
        % then is coherograms (breakpoint evaluated)
        clim = [0.3 0.7];
        addval = 0;
    elseif min(DATSTRUCT.chan1.base(:))>0 & max(DATSTRUCT.chan1.base(:))<0.5
        clim = [0 max(DATSTRUCT.chan1.base(:))];
        addval = 0.2;
    elseif zscoreLFP==3
        clim = [-0.5 0.5];
        addval = 0.2;
    else
        clim = [-8 8]; % zscore
        addval = 1;
    end
    
    % ============ TARG
    indsyl = 1;
    ylabthis = 'TARG';
    % --
    dat1_base = squeeze(DATSTRUCT.chan1.base(:,:,indsyl,:));
    dat2_base = squeeze(DATSTRUCT.chan2.base(:,:,indsyl,:));
    dat1_wn = squeeze(DATSTRUCT.chan1.wn(:,:,indsyl,:));
    dat2_wn = squeeze(DATSTRUCT.chan2.wn(:,:,indsyl,:));
    
    lt_neural_Coher_Learn_PlotSum2_sub
    
    % ============ SAME
    indsyl = 2;
    ylabthis = 'SAME';
    % --
    dat1_base = squeeze(DATSTRUCT.chan1.base(:,:,indsyl,:));
    dat2_base = squeeze(DATSTRUCT.chan2.base(:,:,indsyl,:));
    dat1_wn = squeeze(DATSTRUCT.chan1.wn(:,:,indsyl,:));
    dat2_wn = squeeze(DATSTRUCT.chan2.wn(:,:,indsyl,:));
    
    lt_neural_Coher_Learn_PlotSum2_sub
    
    
    % ============ DIFF
    indsyl = 3;
    ylabthis = 'DIFF';
    % --
    dat1_base = squeeze(DATSTRUCT.chan1.base(:,:,indsyl,:));
    dat2_base = squeeze(DATSTRUCT.chan2.base(:,:,indsyl,:));
    dat1_wn = squeeze(DATSTRUCT.chan1.wn(:,:,indsyl,:));
    dat2_wn = squeeze(DATSTRUCT.chan2.wn(:,:,indsyl,:));
    
    lt_neural_Coher_Learn_PlotSum2_sub
    
    
    
    %% ============== DIFFERENCES (TARGET MINUS NONTARG)
    % ================ TARGET MINUS SAME
    indsyl1 = 1;
    indsyl2 = 2;
    ylabthis = 'TARG - SAME[WNmeans]';
    
    % -
    dat1=  squeeze(DATSTRUCT.chan1.wn(:,:,indsyl1,:) - DATSTRUCT.chan1.wn(:,:,indsyl2,:));
    dat2=  squeeze(DATSTRUCT.chan2.wn(:,:,indsyl1,:) - DATSTRUCT.chan2.wn(:,:,indsyl2,:));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat1, tbins, ffbins, 1, '', clim);
    title('chan1');
    ylabel(ylabthis);
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat1, tbins, ffbins, 3, '', clim);
    title('chan1');
    ylabel(ylabthis);
    xlim(XLIM);
    ylim(YLIM);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat2, tbins, ffbins, 1, '', clim);
    title('chan2');
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat2, tbins, ffbins, 3, '', clim);
    title('chan2');
    xlim(XLIM);
    ylim(YLIM);
    colorbar('East');
    
    % ------------ PLOT TIMECOURSE IN A FEW FREQUENCY BINS
    if plotEachTrial==1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(dat1, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
        title('chan1');
        axis tight;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(dat2, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
        title('chan2');
        axis tight;
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat1, tbins, ffbins, 2, '-', clim, 1, 0);
    title('chan1');
    ylabel(ylabthis);
    axis tight;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat2, tbins, ffbins, 2, '-', clim, 1, 0);
    title('chan2');
    ylabel(ylabthis);
    axis tight;
    
    
    
    % ===============
    dat3 = squeeze((DATSTRUCT.chan1.wn(:,:,indsyl1,:) - DATSTRUCT.chan1.base(:,:,indsyl1,:)) - ...
        (DATSTRUCT.chan1.wn(:,:,indsyl2,:) - DATSTRUCT.chan1.base(:,:,indsyl2,:)));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 1, '', clim);
    title('chan1');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 3, '', clim);
    title('chan1');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
    title('(wn-base) - (wn-base)]');
    ylabel('chan1');
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim, 1, 0);
    title('(wn-base) - (wn-base)]');
    
    % ===============
    dat3 = squeeze((DATSTRUCT.chan2.wn(:,:,indsyl1,:) - DATSTRUCT.chan2.base(:,:,indsyl1,:)) - ...
        (DATSTRUCT.chan2.wn(:,:,indsyl2,:) - DATSTRUCT.chan2.base(:,:,indsyl2,:)));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 1, '', clim);
    title('chan2');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 3, '', clim);
    title('chan2');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
    %     title('(wn-base) - (wn-base)]');
    %     ylabel('chan2');
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim, 1, 0);
    title('(wn-base) - (wn-base)]');
    
    
    
    % ================ TARGET MINUS DIFF
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    indsyl1 = 1;
    indsyl2 = 3;
    ylabthis = 'TARG - DIFF[WNmeans]';
    
    % -
    dat1=  squeeze(DATSTRUCT.chan1.wn(:,:,indsyl1,:) - DATSTRUCT.chan1.wn(:,:,indsyl2,:));
    dat2=  squeeze(DATSTRUCT.chan2.wn(:,:,indsyl1,:) - DATSTRUCT.chan2.wn(:,:,indsyl2,:));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat1, tbins, ffbins, 1, '', clim);
    title('chan1');
    ylabel(ylabthis);
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat1, tbins, ffbins, 3, '', clim);
    title('chan1');
    ylabel(ylabthis);
    xlim(XLIM);
    ylim(YLIM);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat2, tbins, ffbins, 1, '', clim);
    title('chan2');
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat2, tbins, ffbins, 3, '', clim);
    title('chan2');
    xlim(XLIM);
    ylim(YLIM);
%     colorbar('East');
    
    % ------------ PLOT TIMECOURSE IN A FEW FREQUENCY BINS
    if plotEachTrial==1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(dat1, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
        title('chan1');
        axis tight;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(dat2, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
        title('chan2');
        axis tight;
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat1, tbins, ffbins, 2, '-', clim, 1, 0);
    title('chan1');
    ylabel(ylabthis);
    axis tight;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat2, tbins, ffbins, 2, '-', clim, 1, 0);
    title('chan2');
    ylabel(ylabthis);
    axis tight;
    
    
    
    % ===============
    dat3 = squeeze((DATSTRUCT.chan1.wn(:,:,indsyl1,:) - DATSTRUCT.chan1.base(:,:,indsyl1,:)) - ...
        (DATSTRUCT.chan1.wn(:,:,indsyl2,:) - DATSTRUCT.chan1.base(:,:,indsyl2,:)));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 1, '', clim);
    title('chan1');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 3, '', clim);
    title('chan1');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
    title('(wn-base) - (wn-base)]');
    ylabel('chan1');
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim, 1, 0);
    title('(wn-base) - (wn-base)]');
    
    % ===============
    dat3 = squeeze((DATSTRUCT.chan2.wn(:,:,indsyl1,:) - DATSTRUCT.chan2.base(:,:,indsyl1,:)) - ...
        (DATSTRUCT.chan2.wn(:,:,indsyl2,:) - DATSTRUCT.chan2.base(:,:,indsyl2,:)));
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 1, '', clim);
    title('chan2');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 3, '', clim);
    title('chan2');
    ylabel('(wn-base) - (wn-base)]');
    xlim(XLIM);
    ylim(YLIM);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    lt_neural_Coher_Plot(dat3, tbins, ffbins, 2, '-', clim, 1, 0);
    title('(wn-base) - (wn-base)]');
    
    
end
