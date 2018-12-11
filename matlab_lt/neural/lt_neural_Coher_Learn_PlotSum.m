function lt_neural_Coher_Learn_PlotSum(OUTSTRUCT, PARAMS, SwitchStruct, ...
    sumplottype, plotAllSwitchRaw, clim, fieldtoplot)

% sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs

tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

% clim = [-0.15 0.15];

%% ======= field type to plot

if strcmp(fieldtoplot, 'Spec1Mean_WNminusBase')
    fieldtoplot_base = 'Spec1Mean_Base';
    fieldtoplot_WN = 'Spec1Mean_WN';
end


%% EXTRACT (DATAPOINT = CHANNEL PAIR) - IE.. AVERAGES OVER SIMILAR TYPE SYLLABLES
% ======= for each channel pair, get mean for targ, nontarg
% -- gets motifs for each channel pair
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
indsgrp_unique = unique(indsgrp);

allbnum = [];
allenum = [];
allswnum = [];
allCohDiffMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_unique)); % [t, ff, targ-same-diff, cases]

allCohBaseMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_unique)); % [t, ff, targ-same-diff, cases]
allCohWNMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_unique)); % [t, ff, targ-same-diff, cases]

allCohDiffScalar = nan(length(indsgrp_unique), 3);
for i=1:length(indsgrp_unique)
    j=indsgrp_unique(i);
    
    % for this channel pair ...
    % ----- info...
    allbnum = [allbnum; unique(OUTSTRUCT.bnum(indsgrp==j))];
    allenum = [allenum; unique(OUTSTRUCT.enum(indsgrp==j))];
    allswnum = [allswnum; unique(OUTSTRUCT.switch(indsgrp==j))];
    
    % ---------- target
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    % difference matrix
    cohmat = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot)(indsthis)), 3);
    allCohDiffMat(:,:, 1, i) = cohmat;
    % scalar of difference
    if isfield(OUTSTRUCT, 'CohMean_WNminusBase_scalar')
        cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    else
        cohscal = nan;
    end
    allCohDiffScalar(i, 1) = nanmean(cohscal);
    % base and WN means
    cohmat_base = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot_base)(indsthis)), 3);
    cohmat_WN = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot_WN)(indsthis)), 3);
    allCohBaseMat(:,:, 1, i) = cohmat_base;
    allCohWNMat(:,:,1, i) = cohmat_WN;
    
    
    % ---------- same-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    if any(indsthis)
        cohmat = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot)(indsthis)), 3);
        allCohDiffMat(:,:, 2, i) = cohmat;
        if isfield(OUTSTRUCT, 'CohMean_WNminusBase_scalar')
            cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        else
            cohscal = nan;
        end
        allCohDiffScalar(i, 2) = nanmean(cohscal);
    % base and WN means
    cohmat_base = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot_base)(indsthis)), 3);
    cohmat_WN = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot_WN)(indsthis)), 3);
    allCohBaseMat(:,:, 2, i) = cohmat_base;
    allCohWNMat(:,:,2, i) = cohmat_WN;
    end

    % -------- diff type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    if any(indsthis)
        cohmat = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot)(indsthis)), 3);
        allCohDiffMat(:,:, 3, i) = cohmat;
        
        if isfield(OUTSTRUCT, 'CohMean_WNminusBase_scalar')
            cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        else
            cohscal = nan;
        end
        allCohDiffScalar(i, 3) = nanmean(cohscal);
        % base and WN means
    cohmat_base = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot_base)(indsthis)), 3);
    cohmat_WN = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoplot_WN)(indsthis)), 3);
    allCohBaseMat(:,:, 3, i) = cohmat_base;
    allCohWNMat(:,:,3, i) = cohmat_WN;
end
end




%% EXTRACT (DATAPOINT = SWITCH) - IE.. AVERAGES OVER ALL CHANNEL PAIRS FOR A SWITCH

% ==================== PLOT, ONE FOR EACH BIRD (overlay all switches)
indsgrp_switch = lt_tools_grp2idx({allbnum, allenum, allswnum});
indsgrp_switch_unique = unique(indsgrp_switch);
allswitch_bnum = [];
allswitch_enum = [];
allswitch_swnum =[];
allswitch_CohDiffMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_switch_unique)); % [t, ff, targ-same-diff, cases]
allswitch_CohBaseMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_switch_unique)); % [t, ff, targ-same-diff, cases]
allswitch_CohWNMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_switch_unique)); % [t, ff, targ-same-diff, cases]
allswitch_CohDiffScalar = nan(length(indsgrp_switch_unique), 3);

for i=1:length(indsgrp_switch_unique)
    j=indsgrp_switch_unique(i);
    
    indsthis = indsgrp_switch==j; % all channel pairs for this switch
    % ---------------------- INFO
    allswitch_bnum = [allswitch_bnum; unique(allbnum(indsthis))];
    allswitch_enum = [allswitch_enum; unique(allenum(indsthis))];
    allswitch_swnum =[allswitch_swnum; unique(allswnum(indsthis))];
    
    % ================= target
    sylind = 1;
    cohmat = nanmean(squeeze(allCohDiffMat(:,:,sylind,indsthis)), 3);
    allswitch_CohDiffMat(:,:,sylind, i) = cohmat;
    
    cohscal = nanmean(allCohDiffScalar(indsthis,sylind));
    allswitch_CohDiffScalar(i,sylind) = cohscal;
    
    cohmatbase = nanmean(squeeze(allCohBaseMat(:,:,sylind,indsthis)), 3);
    allswitch_CohBaseMat(:,:,sylind, i) = cohmatbase;
    cohmatWN = nanmean(squeeze(allCohWNMat(:,:,sylind,indsthis)), 3);
    allswitch_CohWNMat(:,:,sylind, i) = cohmatWN;
    
    % ================= SAME
    sylind = 2;
    cohmat = nanmean(squeeze(allCohDiffMat(:,:,sylind,indsthis)), 3);
    allswitch_CohDiffMat(:,:,sylind, i) = cohmat;
    
    cohscal = nanmean(allCohDiffScalar(indsthis,sylind));
    allswitch_CohDiffScalar(i,sylind) = cohscal;
    
    cohmatbase = nanmean(squeeze(allCohBaseMat(:,:,sylind,indsthis)), 3);
    allswitch_CohBaseMat(:,:,sylind, i) = cohmatbase;
    cohmatWN = nanmean(squeeze(allCohWNMat(:,:,sylind,indsthis)), 3);
    allswitch_CohWNMat(:,:,sylind, i) = cohmatWN;
    
    % ================= DIFF
    sylind = 3;
    cohmat = nanmean(squeeze(allCohDiffMat(:,:,sylind,indsthis)), 3);
    allswitch_CohDiffMat(:,:,sylind, i) = cohmat;
    
    cohscal = nanmean(allCohDiffScalar(indsthis,sylind));
    allswitch_CohDiffScalar(i,sylind) = cohscal;
    
    cohmatbase = nanmean(squeeze(allCohBaseMat(:,:,sylind,indsthis)), 3);
    allswitch_CohBaseMat(:,:,sylind, i) = cohmatbase;
    cohmatWN = nanmean(squeeze(allCohWNMat(:,:,sylind,indsthis)), 3);
    allswitch_CohWNMat(:,:,sylind, i) = cohmatWN;
end


%% ========================= PLOT ALL DAT ONE FOR EACH SWITCH
if plotAllSwitchRaw ==1
    figcount=1;
    subplotrows=3;
    subplotcols=10;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    indsgrp_switch = lt_tools_grp2idx({allbnum, allenum, allswnum});
    indsgrp_switch_unique = unique(indsgrp_switch);
    for j=indsgrp_switch_unique'
        
        indsthis = indsgrp_switch==j;
        % ---------------------- INFO
        bname = SwitchStruct.bird(unique(allbnum(indsthis))).birdname;
        ename = SwitchStruct.bird(unique(allbnum(indsthis))).exptnum(unique(allenum(indsthis))).exptname;
        swnum = unique(allswnum(indsthis));
        
        % ================= target
        sylind = 1;
        plotitle = 'TARG';
        
        % ----- RUN
        cohmat = squeeze(allCohDiffMat(:,:,sylind,indsthis));
        % 1. cohgram of diff
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
        title(plotitle);
        ylabel({[bname '-' ename '-sw' num2str(swnum)]});
        % 2. ffband differences
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
        lt_plot_zeroline;
        lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
        
        % ================= SAME
        sylind = 2;
        plotitle = 'SAME';
        
        % ----- RUN
        cohmat = squeeze(allCohDiffMat(:,:,sylind,indsthis));
        % 1. cohgram of diff
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
        title(plotitle);
        % 2. ffband differences
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
        lt_plot_zeroline;
        lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
        
        % ================= DIFF
        sylind = 3;
        plotitle = 'DIFF';
        
        % ----- RUN
        cohmat = squeeze(allCohDiffMat(:,:,sylind,indsthis));
        % 1. cohgram of diff
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
        title(plotitle);
        % 2. ffband differences
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
        lt_plot_zeroline;
        lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
        
        
        % ================= target - same
        % ----- RUN
        cohmat1 = nanmean(squeeze(allCohDiffMat(:,:,1,indsthis)),3);
        cohmat2 = nanmean(squeeze(allCohDiffMat(:,:,2,indsthis)), 3);
        cohmat = cohmat1 - cohmat2;
        % 1. cohgram of diff
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
        title('TARG - SAME');
        % 2. ffband differences
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
        lt_plot_zeroline;
        lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
        
        
        % ================= target - DIFF
        % ----- RUN
        cohmat1 = nanmean(squeeze(allCohDiffMat(:,:,1,indsthis)),3);
        cohmat2 = nanmean(squeeze(allCohDiffMat(:,:,3,indsthis)), 3);
        cohmat = cohmat1 - cohmat2;
        % 1. cohgram of diff
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
        title('TARG - DIFF');
        % 2. ffband differences
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
        lt_plot_zeroline;
        lt_plot_text(0, clim(2)-0.05, ['n=' num2str(size(cohmat,3))], 'b');
        
    end
    
end

%% PLOT - ONE PLOT PER SYL TYPE, AVERAGING EITHER OVER SWITCH OR CHANNEL PAIR

if strcmp(sumplottype, 'switches')
    CohMatDat = allswitch_CohDiffMat;
    CohMatDat_base = allswitch_CohBaseMat;
    CohMatDat_WN = allswitch_CohWNMat;
    CohScalDat = allswitch_CohDiffScalar;
    birddat = allswitch_bnum;
elseif strcmp(sumplottype, 'chanpairs')
    CohMatDat = allCohDiffMat;
    CohMatDat_base = allCohBaseMat;
    CohMatDat_WN = allCohWNMat;
    CohScalDat = allCohDiffScalar;
    birddat = allbnum;
end

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
subplotcols=6;
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
colorbar
% 2. ffband differences
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
lt_plot_zeroline;
lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');

% ------------- SAME
cohmat = squeeze(CohMatDat(:,:,2,:));
plottit = 'SAME';
% 1. cohgram of diff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
title(plottit);
colorbar
% 2. ffband differences
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
lt_plot_zeroline;
lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');

% ------------- diff
cohmat = squeeze(CohMatDat(:,:,3,:));
plottit = 'DIFF';
% 1. cohgram of diff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 1, '', clim);
title(plottit);
colorbar
% 2. ffband differences
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
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
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces);
lt_plot_zeroline;
lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
% 2. ffband differences
plotindivtraces =0;
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces);
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
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces);
lt_plot_zeroline;
lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
% 2. ffband differences
plotindivtraces =0;
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1, plotindivtraces);
lt_plot_zeroline;
lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');

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


%% ========== PLOTTING LFP
figcount=1;
subplotrows=3;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

if strcmp(sumplottype, 'switches')
%     LFPmat_base = allswitch_CohDiffMat;
    CohMatDat_base = allswitch_CohBaseMat;
    CohMatDat_WN = allswitch_CohWNMat;
    CohScalDat = allswitch_CohDiffScalar;
    birddat = allswitch_bnum;
elseif strcmp(sumplottype, 'chanpairs')
    CohMatDat = allCohDiffMat;
    CohMatDat_base = allCohBaseMat;
    CohMatDat_WN = allCohWNMat;
    CohScalDat = allCohDiffScalar;
    birddat = allbnum;
end

% ============ TARG
ind3 = 1;
% -- BASE
lfpmat = nanmean(squeeze(CohMatDat_base(:,:, ind3,:)), 3);
plottit = 'TARG';


% 1. cohgram 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title([plottit '[BASE]']);
imagesc(tbins, ffbins, cohmat')
colorbar
% --- WN
cohmat = nanmean(squeeze(CohMatDat_WN(:,:, ind3,:)), 3);
plottit = 'TARG';
% 1. cohgram 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title([plottit '[WN]']);
imagesc(tbins, ffbins, 10*log10(cohmat'))
colorbar


% % 2. ffband differences
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', clim, 1);
% lt_plot_zeroline;
% lt_plot_text(0, clim(2)-0.05, ['n=' num2str(sum(squeeze(~isnan(cohmat(1,1,:)))))], 'b');
% 

