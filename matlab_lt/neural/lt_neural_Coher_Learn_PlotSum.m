function lt_neural_Coher_Learn_PlotSum(OUTSTRUCT, PARAMS, SwitchStruct, ...
    sumplottype, plotAllSwitchRaw, clim)

% sumplottype = 'switches'; % i.e. what is datapoint?
% switches
% chanpairs

tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

% clim = [-0.15 0.15];

%% EXTRACT (DATAPOINT = CHANNEL PAIR) - IE.. AVERAGES OVER SIMILAR TYPE SYLLABLES
% ======= for each channel pair, get mean for targ, nontarg
% -- gets motifs for each channel pair
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
indsgrp_unique = unique(indsgrp);

allbnum = [];
allenum = [];
allswnum = [];
allCohDiffMat = nan(length(tbins), length(ffbins), 3, length(indsgrp_unique)); % [t, ff, targ-same-diff, cases]
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
    cohmat = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis)), 3);
    allCohDiffMat(:,:, 1, i) = cohmat;
        
    cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
    allCohDiffScalar(i, 1) = nanmean(cohscal);
    
    % ---------- same-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    if any(indsthis)
        cohmat = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis)), 3);
        allCohDiffMat(:,:, 2, i) = cohmat;

        cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        allCohDiffScalar(i, 2) = nanmean(cohscal);
end
    
    % -------- diff type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    if any(indsthis)
        cohmat = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis)), 3);
        allCohDiffMat(:,:, 3, i) = cohmat;
        
        cohscal = OUTSTRUCT.CohMean_WNminusBase_scalar(indsthis);
        allCohDiffScalar(i, 3) = nanmean(cohscal);
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
    
    % ================= SAME
    sylind = 2;
    cohmat = nanmean(squeeze(allCohDiffMat(:,:,sylind,indsthis)), 3);
    allswitch_CohDiffMat(:,:,sylind, i) = cohmat;
    
    cohscal = nanmean(allCohDiffScalar(indsthis,sylind));
    allswitch_CohDiffScalar(i,sylind) = cohscal;
    
    % ================= DIFF
    sylind = 3;
    cohmat = nanmean(squeeze(allCohDiffMat(:,:,sylind,indsthis)), 3);
    allswitch_CohDiffMat(:,:,sylind, i) = cohmat;
    
    cohscal = nanmean(allCohDiffScalar(indsthis,sylind));
    allswitch_CohDiffScalar(i,sylind) = cohscal;
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
    CohScalDat = allswitch_CohDiffScalar;
    birddat = allswitch_bnum;
elseif strcmp(sumplottype, 'chanpairs')
    CohMatDat = allCohDiffMat;
    CohScalDat = allCohDiffScalar;
    birddat = allbnum;
end
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




