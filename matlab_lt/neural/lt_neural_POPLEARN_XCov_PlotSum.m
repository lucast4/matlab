function lt_neural_POPLEARN_XCov_PlotSum(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    indtoget_b_e_s, datlevel, fracchange, plotNotminShuff, clim, alignto)
%% lt 1/9/19 - plot summary of xcov during learning.

if ~exist('clim', 'var')
    clim = [-0.05 0.05];
end


%%
if ~isempty(indtoget_b_e_s)
    OUTSTRUCT_XCOV = lt_structure_RmvEmptyField(OUTSTRUCT_XCOV);
    indstokeep = ismember([OUTSTRUCT_XCOV.bnum OUTSTRUCT_XCOV.enum OUTSTRUCT_XCOV.swnum], indtoget_b_e_s, 'rows');
    OUTSTRUCT_XCOV = lt_structure_subsample_all_fields(OUTSTRUCT_XCOV, indstokeep, 1);
end

%% ==== ALIGN TO SYL OR WN?
if strcmp(alignto, 'syl')
    % -- do nothing, deault
elseif strcmp(alignto, 'wn')
   OUTSTRUCT_XCOV.XcovgramBase =  OUTSTRUCT_XCOV.XcovgramBase_alignWN;
   OUTSTRUCT_XCOV.XcovgramWN =  OUTSTRUCT_XCOV.XcovgramWN_alignWN;
   PARAMS.xcenters_gram = PARAMS.xcenters_gram_alignWN;
end

%% COLLECT DATA ACROSS EXPTS
if strcmp(datlevel, 'switch')
    [indsgrp, indsgrpUni] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
        OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.swnum});
elseif strcmp(datlevel, 'neurpair')
    [indsgrp, indsgrpUni] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
        OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.swnum, OUTSTRUCT_XCOV.neurpairnum});
end
%
% figcount=1;
% subplotrows=6;
% subplotcols=2;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
%

XCovDiffAll = nan(size(OUTSTRUCT_XCOV.XcovWN,2), 4, length(indsgrpUni));
XCovBaseAll = nan(size(OUTSTRUCT_XCOV.XcovWN,2), 4, length(indsgrpUni));
XCovWNAll = nan(size(OUTSTRUCT_XCOV.XcovWN,2), 4, length(indsgrpUni));

XCovWNGramAll = cell(length(indsgrpUni), 4);
XCovBaseGramAll = cell(length(indsgrpUni), 4);
XCovDiffGramAll = cell(length(indsgrpUni), 4);
for i=1:length(indsgrpUni)
    
    % ========== expt details
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsgrp==indsgrpUni(i)));
    enum = unique(OUTSTRUCT_XCOV.enum(indsgrp==indsgrpUni(i)));
    swnum = unique(OUTSTRUCT_XCOV.swnum(indsgrp==indsgrpUni(i)));
    
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    %     Yall = []; % targ, same, diff, rest [WN minus base, xcov]
    
    % ################################ TARGET
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==1;
    colthis = 1;
    
    %     ptit = 'TARG';
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
        try
            covgram_base = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
            covgram_wn = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
        catch err
            covgram_base = [];
            covgram_wn = [];
        end
    end
    
    % === for each pair get difference and take average
    if fracchange==1
        y = mean((covWN - covbase)./covbase,1);
    else
        y = mean(covWN - covbase,1);
    end
    XCovDiffAll(:,colthis,i) = y';
    
    % ---- baseline and WN
    ybase = mean(covbase,1);
    ywn = mean(covWN,1);
    
    XCovBaseAll(:,colthis,i) = ybase';
    XCovWNAll(:,colthis,i) = ywn';
    
    % ========== COV GRAM
    ygram = mean(covgram_wn-covgram_base, 3);
    XCovDiffGramAll{i, colthis} = ygram;
    
    XCovBaseGramAll{i,colthis} = mean(covgram_base,3);
    XCovWNGramAll{i, colthis} = mean(covgram_wn, 3);
    
    
    
    % ################################ SAME
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & ...
        OUTSTRUCT_XCOV.issame==1;
    colthis = 2;
    
    %     ptit = 'TARG';
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
        
        try
            covgram_base = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
            covgram_wn = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
        catch err
            covgram_base = [];
            covgram_wn = [];
        end
    end
    
    % === for each pair get difference and take average
    if fracchange==1
        y = mean((covWN - covbase)./covbase,1);
    else
        y = mean(covWN - covbase,1);
    end
    XCovDiffAll(:,colthis,i) = y';
    
    % ---- baseline and WN
    ybase = mean(covbase,1);
    ywn = mean(covWN,1);
    
    XCovBaseAll(:,colthis,i) = ybase';
    XCovWNAll(:,colthis,i) = ywn';
    
    % ========== COV GRAM
    ygram = mean(covgram_wn-covgram_base, 3);
    XCovDiffGramAll{i, colthis} = ygram;
    
    XCovBaseGramAll{i,colthis} = mean(covgram_base,3);
    XCovWNGramAll{i, colthis} = mean(covgram_wn, 3);
    
    
    
    
    % ################################ DIFF
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & ...
        OUTSTRUCT_XCOV.issame==0;
    colthis = 3;
    
    %     ptit = 'TARG';
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
        
        try
            covgram_base = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
            covgram_wn = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
        catch err
            covgram_base = [];
            covgram_wn = [];
        end
    end
    
    % === for each pair get difference and take average
    if fracchange==1
        y = mean((covWN - covbase)./covbase,1);
    else
        y = mean(covWN - covbase,1);
    end
    XCovDiffAll(:,colthis,i) = y';
    
    % ---- baseline and WN
    ybase = mean(covbase,1);
    ywn = mean(covWN,1);
    
    XCovBaseAll(:,colthis,i) = ybase';
    XCovWNAll(:,colthis,i) = ywn';
    
    % ========== COV GRAM
    ygram = mean(covgram_wn-covgram_base, 3);
    XCovDiffGramAll{i, colthis} = ygram;
    
    XCovBaseGramAll{i,colthis} = mean(covgram_base,3);
    XCovWNGramAll{i, colthis} = mean(covgram_wn, 3);
    
    
    
    
    % ################################ NONTARG
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0;
    
    colthis = 4;
    
    %     ptit = 'TARG';
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
        
        try
            covgram_base = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
            covgram_wn = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
        catch err
            covgram_base = [];
            covgram_wn = [];
        end
    end
    
    % === for each pair get difference and take average
    if fracchange==1
        y = mean((covWN - covbase)./covbase,1);
    else
        y = mean(covWN - covbase,1);
    end
    XCovDiffAll(:,colthis,i) = y';
    
    % ---- baseline and WN
    ybase = mean(covbase,1);
    ywn = mean(covWN,1);
    
    XCovBaseAll(:,colthis,i) = ybase';
    XCovWNAll(:,colthis,i) = ywn';
    
    % ========== COV GRAM
    ygram = mean(covgram_wn-covgram_base, 3);
    XCovDiffGramAll{i, colthis} = ygram;
    
    XCovBaseGramAll{i,colthis} = mean(covgram_base,3);
    XCovWNGramAll{i, colthis} = mean(covgram_wn, 3);
    
    
    
    % ================ SAVE OUTPUT
    
end


%% sampel sizes
[~, tmp] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.swnum});
Nsw = length(tmp);


[~, tmp] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.swnum, OUTSTRUCT_XCOV.neurpairnum});

Nneur = length(tmp);


%% ======= [PLOT] Baseline and WN, overlaid
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ############################################################### TARG
colthis = 1;
x = PARAMS.Xcov_ccLags;

% ======== BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [BASE]');
ymat = squeeze(XCovBaseAll(:,colthis,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== WN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [WN]');
ymat = squeeze(XCovWNAll(:,colthis,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== OVERLAY WN AND BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [WN(r), base(k)]');

% -- base
ymat = squeeze(XCovBaseAll(:,colthis,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
end
% -- WN
ymat = squeeze(XCovWNAll(:,colthis,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
% -- formating
lt_plot_zeroline;



% ############################################################### TARG
colthis = 1;
x = PARAMS.Xcov_ccLags;

% ======== BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [BASE] [ONLY IF HAVE SAME]');
indstmp = ~isnan(squeeze(XCovBaseAll(1,2,:)));
ymat = squeeze(XCovBaseAll(:,colthis,indstmp));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== WN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [WN]');
indstmp = ~isnan(squeeze(XCovWNAll(1,2,:)));
ymat = squeeze(XCovWNAll(:,colthis,indstmp));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== OVERLAY WN AND BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [WN(r), base(k)]');

% -- base
ymat = squeeze(XCovBaseAll(:,colthis,indstmp));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
end
% -- WN
ymat = squeeze(XCovWNAll(:,colthis,indstmp));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
% -- formating
lt_plot_zeroline;

% ############################################################### SAME
colthis = 2;
x = PARAMS.Xcov_ccLags;

% ======== BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('SAME [BASE]');
ymat = squeeze(XCovBaseAll(:,colthis,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== WN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('[WN]');
ymat = squeeze(XCovWNAll(:,colthis,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== OVERLAY WN AND BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('[WN(r), base(k)]');

% -- base
ymat = squeeze(XCovBaseAll(:,colthis,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
end
% -- WN
ymat = squeeze(XCovWNAll(:,colthis,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
% -- formating
lt_plot_zeroline;



% ############################################################### DIFF
colthis = 3;
x = PARAMS.Xcov_ccLags;

% ======== BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('DIFF [BASE]');
ymat = squeeze(XCovBaseAll(:,colthis,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== WN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('[WN]');
ymat = squeeze(XCovWNAll(:,colthis,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== OVERLAY WN AND BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('[WN(r), base(k)]');

% -- base
ymat = squeeze(XCovBaseAll(:,colthis,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
end
% -- WN
ymat = squeeze(XCovWNAll(:,colthis,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
% -- formating
lt_plot_zeroline;


linkaxes(hsplots, 'xy');


%% ===== TARG GREATER THAN DIFF, WN END

% ======== OVERLAY WN AND BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('[WNtarg(r), WNdiff(k)]');

% --
ymat = squeeze(XCovWNAll(:,1,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
% -- WN
ymat = squeeze(XCovWNAll(:,3,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
end
% -- formating
lt_plot_zeroline;

% ======== OVERLAY WN AND BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('[WNtarg(r)-WNdiff(k)]');

% --
ymat = squeeze(XCovWNAll(:,1,:)) - squeeze(XCovWNAll(:,3,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'b'}, 1);
end
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), mean(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- formating
lt_plot_zeroline;
% ---





linkaxes(hsplots, 'xy');



%% ===== [PLOT] BASELINE AND WN, COVGRAM
try
    figcount=1;
    subplotrows=3;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    indsgood = ~cellfun(@isempty, XCovDiffGramAll(:,2));
    
    % #################################### TARG
    colthis = 1;
    
    % == base
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('TARG [BASE]');
    
    indsgood = ~cellfun(@isempty, XCovDiffGramAll(:,colthis));
    y = lt_neural_Coher_Cell2Mat(XCovBaseGramAll(indsgood,colthis)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    % == WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('[WN]');
    
    indsgood = ~cellfun(@isempty, XCovWNGramAll(:,colthis));
    y = lt_neural_Coher_Cell2Mat(XCovWNGramAll(indsgood,colthis)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    
    % #################################### SAME
    colthis = 2;
    
    % == base
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('SAME [BASE]');
    
    indsgood = ~cellfun(@isempty, XCovDiffGramAll(:,colthis));
    y = lt_neural_Coher_Cell2Mat(XCovBaseGramAll(indsgood,colthis)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    % == WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('[WN]');
    
    indsgood = ~cellfun(@isempty, XCovWNGramAll(:,colthis));
    y = lt_neural_Coher_Cell2Mat(XCovWNGramAll(indsgood,colthis)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    
    % #################################### DIFF
    colthis = 3;
    
    % == base
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('DIFF [BASE]');
    
    indsgood = ~cellfun(@isempty, XCovDiffGramAll(:,colthis));
    y = lt_neural_Coher_Cell2Mat(XCovBaseGramAll(indsgood,colthis)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    % == WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('[WN]');
    
    indsgood = ~cellfun(@isempty, XCovWNGramAll(:,colthis));
    y = lt_neural_Coher_Cell2Mat(XCovWNGramAll(indsgood,colthis)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
catch err
end

%% ==== PLOT

lt_figure; hold on;
hsplots = [];

% ======== TARG
hsplot = lt_subplot(3,2,1); hold on;
hsplots = [hsplots; hsplot];
title('TARG');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,1,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end
% -- sampel size
lt_plot_annotation(1, ['N=' num2str(Nsw) 'expt, ' num2str(Nneur) 'neur'], 'b');


% ======== SAME
hsplot = lt_subplot(3,2,2); hold on;
hsplots = [hsplots; hsplot];
title('SAME');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,2,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end



% ======== DIFF
hsplot = lt_subplot(3,2,3); hold on;
hsplots = [hsplots; hsplot];
title('DIFF');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,3,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end



% ======== NONTARG
hsplot = lt_subplot(3,2,4); hold on;
hsplots = [hsplots; hsplot];
title('NONTARG');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,4,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end




% ======== TARG - NONTARG
hsplot = lt_subplot(3,2,5); hold on;
hsplots = [hsplots; hsplot];
title('TARG - NONTARG');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,1,:) - XCovDiffAll(:,4,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end


% ======== TARG - DIFF
hsplot = lt_subplot(3,2,6); hold on;
hsplots = [hsplots; hsplot];
title('TARG - DIFF');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,1,:) - XCovDiffAll(:,3,:));

plot(x, ymat, '-k');
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.15
        lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
    end
end


linkaxes(hsplots, 'xy');
%% ==== PLOT [JKUST MEANS]

lt_figure; hold on;
hsplots = [];

% ======== TARG
hsplot = lt_subplot(3,2,1); hold on;
hsplots = [hsplots; hsplot];
title('TARG');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,1,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
axis tight;
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
    end
end

% ======== TARG
hsplot = lt_subplot(3,2,2); hold on;
hsplots = [hsplots; hsplot];
title('TARG');
x = PARAMS.Xcov_ccLags;
indstmp = ~isnan(squeeze(XCovDiffAll(1,2,:)));
ymat = squeeze(XCovDiffAll(:,1,indstmp));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
axis tight;
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
    end
end



% ======== SAME
hsplot = lt_subplot(3,2,3); hold on;
hsplots = [hsplots; hsplot];
title('SAME');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,2,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
axis tight;
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
    end
end



% ======== DIFF
hsplot = lt_subplot(3,2,4); hold on;
hsplots = [hsplots; hsplot];
title('DIFF');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,3,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
axis tight;
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
    end
end



% ======== NONTARG
% hsplot = lt_subplot(3,2,4); hold on;
% hsplots = [hsplots; hsplot];
% title('NONTARG');
% x = PARAMS.Xcov_ccLags;
% ymat = squeeze(XCovDiffAll(:,4,:));
% ymean = nanmean(ymat,2);
% ysem = lt_sem(ymat');
% if length(ysem)>1
%     shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
% end
% axis tight;
% lt_plot_zeroline;
% % --- calcualte p value in each time bin
% for j=1:size(ymat, 1)
%     p = signrank(ymat(j,:));
%     if p<0.1
%         lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
%     end
% end



%
% % ======== TARG - NONTARG
% lt_subplot(3,2,5); hold on;
% title('TARG - NONTARG');
% x = PARAMS.Xcov_ccLags;
% ymat = squeeze(XCovDiffAll(:,1,:) - XCovDiffAll(:,4,:));
% ymean = nanmean(ymat,2);
% ysem = lt_sem(ymat');
% shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
% lt_plot_zeroline;
% % --- calcualte p value in each time bin
% for j=1:size(ymat, 1)
%     p = signrank(ymat(j,:));
%     if p<0.15
%         lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
%     end
% end

% ======== TARG - SAME
hsplot = lt_subplot(3,2,5); hold on;
hsplots = [hsplots; hsplot];
title('TARG - SAME');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,1,:) - XCovDiffAll(:,2,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
axis tight;
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
    end
end


% ======== TARG - DIFF
hsplot = lt_subplot(3,2,6); hold on;
hsplots = [hsplots; hsplot];
title('TARG - DIFF');
x = PARAMS.Xcov_ccLags;
ymat = squeeze(XCovDiffAll(:,1,:) - XCovDiffAll(:,3,:));
ymean = nanmean(ymat,2);
ysem = lt_sem(ymat');
if length(ysem)>1
    shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
end
axis tight;
lt_plot_zeroline;
% --- calcualte p value in each time bin
for j=1:size(ymat, 1)
    p = signrank(ymat(j,:));
    if p<0.1
        lt_plot_text(x(j), ymean(j), ['p=' num2str(p)], 'm');
    end
end


% % ======== TARG - NONTARG
% lt_subplot(3,2,6); hold on;
% title('TARG - DIFF');
% x = PARAMS.Xcov_ccLags;
% ymat = squeeze(XCovDiffAll(:,1,:) - XCovDiffAll(:,3,:));
% ymean = nanmean(ymat,2);
% ysem = lt_sem(ymat');
% shadedErrorBar(x, ymean, ysem, {'Color', 'r'}, 1);
% lt_plot_zeroline;
% % --- calcualte p value in each time bin
% for j=1:size(ymat, 1)
%     p = signrank(ymat(j,:));
%     if p<0.15
%         lt_plot_text(x(j), max(ymat(j,:)), ['p=' num2str(p)], 'r');
%     end
% end
%


linkaxes(hsplots, 'xy');
%% ============ [PLOTS] xcovgram
try
    lt_figure; hold on;
    
    
    indsgood = ~cellfun(@isempty, XCovDiffGramAll(:,2));
    
    % ======== TARG
    lt_subplot(3,2,1); hold on;
    title('TARG');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(:,1)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    lt_subplot(3,2,2); hold on;
    title('TARG');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(:,1)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    % ======== SAME
    lt_subplot(3,2,3); hold on;
    title('SAME');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(indsgood,2)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    lt_subplot(3,2,4); hold on;
    title('SAME');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(indsgood,2)); % time, lags, n
    if ~isempty(y)
        lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    end
    axis tight;
    lt_plot_zeroline;
    
    % ======== DIFF
    lt_subplot(3,2,5); hold on;
    title('DIFF');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(:,3)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    lt_subplot(3,2,6); hold on;
    title('DIFF');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(:,3)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    
    
    lt_figure; hold on;
    % ======== TARG MINUS SAME
    lt_subplot(3,2,1); hold on;
    title('TARG - SAME');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(indsgood,1)) - ...
        lt_neural_Coher_Cell2Mat(XCovDiffGramAll(indsgood,2)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    lt_subplot(3,2,2); hold on;
    title('TARG - SAME');
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
    
    
    % ======== TARG MINUS DIFF
    lt_subplot(3,2,3); hold on;
    title('TARG - DIFF');
    y = lt_neural_Coher_Cell2Mat(XCovDiffGramAll(:,1)) - ...
        lt_neural_Coher_Cell2Mat(XCovDiffGramAll(:,3)); % time, lags, n
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    lt_subplot(3,2,4); hold on;
    title('TARG - DIFF');
    lt_neural_Coher_Plot(y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, []);
    axis tight;
    lt_plot_zeroline;
catch err
end



