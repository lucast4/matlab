function lt_neural_POPLEARN_Xcov_Sumcorr(OUTSTRUCT_XCOV, SummaryStruct, ...
    PARAMS, SwitchStruct)
%% lt 2/16/19 - summarize xcov at baseline (pool by motifID)

XLIM = [-0.02 0.02]; % only for the survey plot

%% get moifID for all cases
motifidall = nan(length(OUTSTRUCT_XCOV.bnum), 1);
for i=1:length(OUTSTRUCT_XCOV.bnum)
    bname = SummaryStruct.birds(OUTSTRUCT_XCOV.bnum(i)).birdname;
    motifname = OUTSTRUCT_XCOV.motifname{i};
    
    motifidall(i) = lt_neural_QUICK_MotifID(bname, motifname);
end

assert(~any(isnan(motifidall)));

OUTSTRUCT_XCOV.motifID = motifidall;


%% ================ SEPARATE BY MOTIF X UNIT...


[indsgrp_units, indsgrp_units_U] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.neurpair}); % one idx for each pair of units
indsgrp_motif = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.motifID}); % one idx for each pair of units

maxbird = max(OUTSTRUCT_XCOV.bnum);
x = PARAMS.Xcov_ccLags;

% === for each bird get list of motifs that exist
motiflist_bird = cell(1, maxbird);
for i=1:maxbird
    motiflist = unique(OUTSTRUCT_XCOV.motifID(OUTSTRUCT_XCOV.bnum==i));
    motiflist_bird{i}=motiflist;
end
for i=1:maxbird
    
    figcount=1;
    subplotrows=1;
    subplotcols=10;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    motiflistthis = motiflist_bird{i};
    bname = SummaryStruct.birds(i).birdname;
    [~, motiflist_out, ~] = ...
        lt_neural_QUICK_MotifID(bname);
    % === go thru each neuron....
    for ii=1:length(indsgrp_units_U)
        
        % === go thru each motif, and pliot if data exists
        count = 0;
        for iii=1:length(motiflistthis)
            motifIDthis = motiflistthis(iii);
            
            indsthis = indsgrp_units==indsgrp_units_U(ii) & OUTSTRUCT_XCOV.bnum==i ...
                & OUTSTRUCT_XCOV.motifID==motifIDthis;
            
            if ~any(indsthis)
                continue
            end
            
            % === initiate figure
            if count==0
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                
                % --- label
                npair = OUTSTRUCT_XCOV.neurpair(indsthis, :);
                title([bname '-' num2str(npair)]);
                count=1;
            end
            
            % ################################## PLOT FOR THIS PAIR/MOTIF
            assert(sum(indsthis)==1);
            y = OUTSTRUCT_XCOV.XcovBase(indsthis, :);
            
            % ==== make height indicate the motif
            y = y+iii/3;
            line([x(1) x(end)], [iii/3 iii/3], 'Color', [0.7 0.2 0.2]);
            plot(x, y, '-k', 'LineWidth', 2);
            
            % === text annotate motif
            lt_plot_text(0, max(y), motiflist_out{motifIDthis}, 'r', 8)
        end
        %         set(gca, 'XTick', motiflistthis./3, 'XTickLabel', motiflist_out(motiflistthis));
        axis tight;
            xlim(XLIM);
        lt_plot_zeroline_vert;
    end
    YLIM = ylim;
    ylim([YLIM(1)-0.3 YLIM(2)+0.3]);
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
    xlim(XLIM);
end

%% OVERLAY ALL BASELINE, FOR EACH MOTIF (UNIQUE);

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.motifID});
tlags = PARAMS.Xcov_ccLags;

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

all_xcov = [];
all_bnum = [];

for i=1:length(indsgrpU)
    indsthis = indsgrp==indsgrpU(i);
    
    covmat = OUTSTRUCT_XCOV.XcovBase(indsthis, :);
    %     covmat = OUTSTRUCT_XCOV.XcovWN(indsthis, :);
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis, :));
    mID = unique(OUTSTRUCT_XCOV.motifID(indsthis, :));
    bname = SummaryStruct.birds(bnum).birdname;
    [~, mlist] = lt_neural_QUICK_MotifID(bname);
    motifname = mlist{mID};
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([bname '-' motifname]);
    ylabel('trace = chanpair');
    plot(tlags, covmat', '-', 'Color', [0.7 0.7 0.7]);
    % -- plot mean
    if size(covmat,1)>1
        ymean = mean(covmat,1);
        ysem = lt_sem(covmat);
        shadedErrorBar(tlags, ymean, ysem, {'Color', 'r'}, 1);
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    all_xcov = [all_xcov; mean(covmat,1)];
    all_bnum = [all_bnum; bnum];
    
end

linkaxes(hsplots, 'xy');

% ========== one mean for each bird
for i=1:max(all_bnum)
    indsthis = all_bnum==i;
    if ~any(indsthis)
        continue
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['bird' num2str(i)]);
    %     ylabel();
    covmat = all_xcov(indsthis, :);
    
    plot(tlags, covmat', '-', 'Color', [0.7 0.7 0.7]);
    % -- plot mean
    if size(covmat,1)>1
        ymean = mean(covmat,1);
        ysem = lt_sem(covmat);
        shadedErrorBar(tlags, ymean, ysem, {'Color', 'r'}, 1);
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

% ====== MEAN OF BIRD MEANS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['all']);
ylabel('trace = bird');
allbird_xcov = [];
pcols = lt_make_plot_colors(max(all_bnum), 0,0);
for i=1:max(all_bnum)
    indsthis = all_bnum==i;
    if ~any(indsthis)
        continue
    end
    
    covmat = all_xcov(indsthis, :);
    
    ymean = mean(covmat,1);
    ysem = lt_sem(covmat);
    shadedErrorBar(tlags, ymean, ysem, {'Color', pcols{i}}, 1);
    
    allbird_xcov = [allbird_xcov; ymean];
end
% ------- grand mean
% shadedErrorBar(tlags, mean(allbird_xcov), lt_sem(allbird_xcov), {'Color', 'r'}, 1);

axis tight;
lt_plot_zeroline;
lt_plot_zeroline_vert;



% ==== mean over all motifs
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['all']);
ylabel('trace = motif');
covmat = all_xcov;

plot(tlags, covmat', '-', 'Color', [0.7 0.7 0.7]);

% -- plot mean
ymean = mean(covmat,1);
ysem = lt_sem(covmat);
shadedErrorBar(tlags, ymean, ysem, {'Color', 'r'}, 1);


linkaxes(hsplots, 'xy');
%% ===================== ONE MEAN FOR EACH BIRD (MEAN OF MOTIF X CHAN PAIRS);


[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum});
tlags = PARAMS.Xcov_ccLags;

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(indsgrpU)
    indsthis = indsgrp==indsgrpU(i);
    
    covmat = OUTSTRUCT_XCOV.XcovBase(indsthis, :);
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis, :));
    bname = SummaryStruct.birds(bnum).birdname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([bname]);
    
    plot(tlags, covmat', '-', 'Color', [0.7 0.7 0.7]);
    % -- plot mean
    if size(covmat,1)>1
        ymean = mean(covmat,1);
        ysem = lt_sem(covmat);
        shadedErrorBar(tlags, ymean, ysem, {'Color', 'r'}, 1);
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

linkaxes(hsplots, 'xy');


%% ================ ONE MEAN FOR EACH BIRD (mean of motifs (i.e. 1 trace per motif))


[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.motifID});

all_xcov = [];
all_bnum = [];

for i=1:length(indsgrpU)
    indsthis = indsgrp==indsgrpU(i);
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsthis, :));
    covmat = OUTSTRUCT_XCOV.XcovBase(indsthis, :);
    
    all_xcov = [all_xcov; mean(covmat,1)];
    all_bnum = [all_bnum; bnum];
    
end

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:max(all_bnum)
    indsthis = all_bnum==i;
    if ~any(indsthis)
        continue
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['bird' num2str(i)]);
    ylabel('trace = motif');
    covmat = all_xcov(indsthis, :);
    
    plot(tlags, covmat', '-', 'Color', [0.7 0.7 0.7]);
    % -- plot mean
    if size(covmat,1)>1
        ymean = mean(covmat,1);
        ysem = lt_sem(covmat);
        shadedErrorBar(tlags, ymean, ysem, {'Color', 'r'}, 1);
    end
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

% ==== mean over all motifs
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['all birds']);
ylabel('trace = motif');
covmat = all_xcov;

plot(tlags, covmat', '-', 'Color', [0.7 0.7 0.7]);

% -- plot mean
ymean = mean(covmat,1);
ysem = lt_sem(covmat);
shadedErrorBar(tlags, ymean, ysem, {'Color', 'r'}, 1);


linkaxes(hsplots, 'xy');