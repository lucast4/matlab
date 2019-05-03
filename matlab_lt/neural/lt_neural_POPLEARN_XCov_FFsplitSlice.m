function lt_neural_POPLEARN_XCov_FFsplitSlice(OUTSTRUCT_XCOV, PARAMS, epochstoplot, ploteachsyl)
%%
if ~exist('ploteachsyl', 'var')
    ploteachsyl = 0; % default is 0, only plot targets
end

%% lt 3/2019 -

% epochstoplot = [1:3] then will avearge over them,..

%% extract data

covbase = OUTSTRUCT_XCOV.Xcovslice_ffsplits_base;
covWN = OUTSTRUCT_XCOV.Xcovslice_ffsplits_epochs;

assert(length(covbase{1})==2, 'assumes that split into 2 ff (i.e hi and lo)');

% ---- convert into n x 2 cell
covbase = [cellfun(@(x)x{1}, covbase, 'UniformOutput', 0) ...
    cellfun(@(x)x{2}, covbase, 'UniformOutput', 0)];

covWN = [cellfun(@(x)mean(x{1}(:,:, epochstoplot), 3), covWN, 'UniformOutput', 0) ...
    cellfun(@(x)mean(x{2}(:,:, epochstoplot), 3), covWN, 'UniformOutput', 0)];


%% =================== FOR EACH CASE FLIP SO IS IN ADAPTIVE DIRECTION
ldir_targ = OUTSTRUCT_XCOV.learndirTarg;

covbase(ldir_targ==-1, :) = fliplr(covbase(ldir_targ==-1, :));
covWN(ldir_targ==-1, :) = fliplr(covWN(ldir_targ==-1, :));


%% ===================== PLOT EACH EXPERIMENT
if ploteachsyl==1
    [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum ...
        OUTSTRUCT_XCOV.switch, OUTSTRUCT_XCOV.motifID_unique});
else
    [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum ...
        OUTSTRUCT_XCOV.switch});
end

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



for i=1:length(indsgrpU)
    
    if ploteachsyl ==1
        indsthis = indsgrp==indsgrpU(i);
        motifname = OUTSTRUCT_XCOV.motifname(indsthis); motifname = motifname{1};
    elseif ploteachsyl==0
        indsthis = indsgrp==indsgrpU(i) & OUTSTRUCT_XCOV.istarg==1;
    end
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis));
    enum = unique(OUTSTRUCT_XCOV.enum(indsthis));
    sw = unique(OUTSTRUCT_XCOV.switch(indsthis));
    
    
    % ============ 1)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline');
    ylabel('(k=nonad, r=adaptive)');
    if ploteachsyl==1
        xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw) '-' motifname]);
    elseif ploteachsyl ==0
        xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw)]);
    end
    x = PARAMS.Xcov_ccLags;
    
    y1 = cell2mat(covbase(indsthis, 1));
    y2 = cell2mat(covbase(indsthis, 2));
    plot(x, y1', '-k');
    plot(x, y2', '-r');
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ============ 2) WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('WN');
    
    x = PARAMS.Xcov_ccLags;
    
    y1 = cell2mat(covWN(indsthis, 1));
    y2 = cell2mat(covWN(indsthis, 2));
    plot(x, y1', '-k');
    plot(x, y2', '-r');
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ============ 2) WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('WN-base');
    
    x = PARAMS.Xcov_ccLags;
    
    y1 = cell2mat(covWN(indsthis, 1)) - cell2mat(covbase(indsthis, 1));
    y2 = cell2mat(covWN(indsthis, 2)) - cell2mat(covbase(indsthis, 2));
    plot(x, y1', '-k');
    plot(x, y2', '-r');
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    % ============ 2) WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('WN-base');
    
    x = PARAMS.Xcov_ccLags;
    
    y1 = cell2mat(covWN(indsthis, 1)) - cell2mat(covbase(indsthis, 1));
    y2 = cell2mat(covWN(indsthis, 2)) - cell2mat(covbase(indsthis, 2));
    if size(y1,1)>1
        shadedErrorBar(x, mean(y1), lt_sem(y1), {'Color', 'k'},1);
        shadedErrorBar(x, mean(y2), lt_sem(y2), {'Color', 'r'},1)
    end
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
end

%% =============== PLOT OVERALL SUMMARY
if ploteachsyl==1
indsthis = 1:length(OUTSTRUCT_XCOV.bnum);
elseif ploteachsyl==0
indsthis = OUTSTRUCT_XCOV.istarg==1;
end

% ============ 1)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('baseline');
ylabel('(k=nonad, r=adaptive)');
%     xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw)]);
x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covbase(indsthis, 1));
y2 = cell2mat(covbase(indsthis, 2));
shadedErrorBar(x, mean(y1), lt_sem(y1), {'Color', 'k'},1);
shadedErrorBar(x, mean(y2), lt_sem(y2), {'Color', 'r'},1)
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ============ 2) WN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('WN');

x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covWN(indsthis, 1));
y2 = cell2mat(covWN(indsthis, 2));
shadedErrorBar(x, nanmean(y1), lt_sem(y1), {'Color', 'k'},1);
shadedErrorBar(x, nanmean(y2), lt_sem(y2), {'Color', 'r'},1)
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ============ 2) WN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('WN-base');

x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covWN(indsthis, 1)) - cell2mat(covbase(indsthis, 1));
y2 = cell2mat(covWN(indsthis, 2)) - cell2mat(covbase(indsthis, 2));
if size(y1,1)>1
    shadedErrorBar(x, mean(y1), lt_sem(y1), {'Color', 'k'},1);
    shadedErrorBar(x, mean(y2), lt_sem(y2), {'Color', 'r'},1)
end
plot(x, mean((y1+y2)./2), '-b')
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;

% =================== NONADAPTIVE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('NONADAPTIVE TRIALS');
ylabel('r = WN; k= base')
x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covbase(indsthis, 1));
y2 = cell2mat(covWN(indsthis, 1));
shadedErrorBar(x, nanmean(y1), lt_sem(y1), {'Color', 'k'},1);
shadedErrorBar(x, nanmean(y2), lt_sem(y2), {'Color', 'r'},1)
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;

% =================== NONADAPTIVE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('NONADAPTIVE TRIALS');
ylabel('(wn minus base)')
x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covbase(indsthis, 1));
y2 = cell2mat(covWN(indsthis, 1));
y = y2-y1;
shadedErrorBar(x, nanmean(y), lt_sem(y), {'Color', 'b'},1);
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;


% =================== ADAPTIVE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('ADAPTIVE TRIALS');
ylabel('r = WN; k= base')
x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covbase(indsthis, 2));
y2 = cell2mat(covWN(indsthis, 2));
shadedErrorBar(x, nanmean(y1), lt_sem(y1), {'Color', 'k'},1);
shadedErrorBar(x, nanmean(y2), lt_sem(y2), {'Color', 'r'},1)
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;


% =================== ADAPTIVE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('ADAPTIVE TRIALS');
ylabel(' (wn minus base)')
x = PARAMS.Xcov_ccLags;

y1 = cell2mat(covbase(indsthis, 2));
y2 = cell2mat(covWN(indsthis, 2));
y = y2-y1;
shadedErrorBar(x, nanmean(y), lt_sem(y), {'Color', 'b'},1);
axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;

linkaxes(hsplots, 'xy');

xlim([-0.04 0.04]); ylim([-0.3 0.5]);
linkaxes(hsplots, 'xy');


%% ===================== ADAPTIVE MINUS NONADAPTIVE
if ploteachsyl==1
    [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum ...
        OUTSTRUCT_XCOV.switch, OUTSTRUCT_XCOV.motifID_unique});
else
    [indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum ...
        OUTSTRUCT_XCOV.switch});
end

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



for i=1:length(indsgrpU)
    
    if ploteachsyl ==1
        indsthis = indsgrp==indsgrpU(i);
        motifname = OUTSTRUCT_XCOV.motifname(indsthis); motifname = motifname{1};
    elseif ploteachsyl==0
        indsthis = indsgrp==indsgrpU(i) & OUTSTRUCT_XCOV.istarg==1;
    end
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis));
    enum = unique(OUTSTRUCT_XCOV.enum(indsthis));
    sw = unique(OUTSTRUCT_XCOV.switch(indsthis));
    
    
    % ============ 1)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline');
    ylabel('adaptive - nonad');
    if ploteachsyl==1
        xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw) '-' motifname]);
    elseif ploteachsyl ==0
        xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw)]);
    end
    x = PARAMS.Xcov_ccLags;
    
    y1 = cell2mat(covbase(indsthis, 1));
    y2 = cell2mat(covbase(indsthis, 2));
    y = y2-y1;
    plot(x, y', '-b');
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    % ============ 1)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('baseline');
    ylabel('adaptive - nonad');
    xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw)]);
    x = PARAMS.Xcov_ccLags;
    
    if size(y,1)>1
        shadedErrorBar(x, mean(y,1), lt_sem(y), {'Color', 'b'}, 1);
    end
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    
    % ============ 1)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('WN');
    ylabel('adaptive - nonad');
    xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw)]);
    x = PARAMS.Xcov_ccLags;
    
    y1 = cell2mat(covWN(indsthis, 1));
    y2 = cell2mat(covWN(indsthis, 2));
    y = y2-y1;
    plot(x, y', '-b');
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    % ============ 1)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('WN');
    ylabel('adaptive - nonad');
    xlabel([num2str(bnum) '-' num2str(enum) '-' num2str(sw)]);
    x = PARAMS.Xcov_ccLags;
    
    if size(y,1)>1
        shadedErrorBar(x, mean(y,1), lt_sem(y), {'Color', 'b'}, 1);
    end
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    
end

linkaxes(hsplots, 'xy');
ylim([-0.8 0.8]);
