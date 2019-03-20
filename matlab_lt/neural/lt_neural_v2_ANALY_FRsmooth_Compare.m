function lt_neural_v2_ANALY_FRsmooth_Compare(OUTDAT, SwitchStruct, ...
    analytype, syltypesneeded, premotorwind, minmotifs, datlevel)
%%

dattype = 'neuron'; % for breaking down shuffle into lower level data.
% dattype = , bird, 'switch', 'expt', 'neuron'
%     nshuffs = 500;

%% ####################### GETTING DEVIATION FROM BASELINE SMOOTHED FR

% prctile_divs = [33 66 100]; % percentiles to divide up data by
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
% prctile_divs = [50 100]; % percentiles to divide up data by
% epochtoplot = 2; % i.e. out of the epochs decided by prctile_divs

% usepercent = 0; % if 1, then gets fr percent deviation from baseline. otherwise uses hz diff
% nbasetime = 60; % 60 minutes bnefore first WN trial is min time for baseline
% nbasetime = []; % 60 minutes bnefore first WN trial is min time for baseline

% analytype = 'AllMinusAll_FRsmooth';
% % AllMinusAll_FRsmooth
% % AllOnlyMinusDiff_FRsmooth
% % AllMinusAllMinusDiff_FRsmooth
% % AllOnlyMinusBase_FRsmooth
%
%
% syltypesneeded = [1 1 1] means needs minimum 1 targ, 1 same, 1 diff.
% [1 0 1] means doesnt care if has same type


%% ==== MA KE SURE DIMENSIONS ARE CORRECT

%% ###################### LIMIT TO NEURONS THAT CONTAIN ALL SYL TYPES?

% =-==== go thru all switches. if bad then throw ou
maxneur = max(OUTDAT.All_neurnum);
numbirds = length(SwitchStruct.bird);

% ======================== SHUFFLE SYL TYPE? % within each neuron
indstokeep = [];
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if isempty(indsthis)
                    continue
                end
                
                % ==== count how many of each syl type there exists
                numtarg = sum(OUTDAT.All_istarg(indsthis)==1);
                numsame = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==1);
                numdiff = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==0);
                
                if all([numtarg numsame numdiff] >= syltypesneeded)
                    % then keep
                    disp([numtarg numsame numdiff]);
                    disp(indsthis)
                    indstokeep = [indstokeep; indsthis];
                end
                
            end
        end
    end
end

disp(['Keeping ' num2str(length(indstokeep)) '/' num2str(length(OUTDAT.All_birdnum)) ' datapoints, passes syltypes required criterion']);

OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);


%% =================== ONLY KEEP SWITCHES THAT HAVE A MINIMUM NUMBER OF MOTIFS

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTDAT.All_birdnum, OUTDAT.All_exptnum, OUTDAT.All_swnum, OUTDAT.All_neurnum});
indstoremove = [];
for i=indsgrpU'
    
    indsthis = indsgrp==i;
    
    nmot = length(unique(OUTDAT.All_motifnum(indsthis)));
    if nmot < minmotifs
        indstoremove = [indstoremove; find(indsthis)];
    end
    
end

indstokeep = ~ismember(1:length(OUTDAT.All_birdnum), indstoremove');
disp(['Keeping ' num2str(sum(indstokeep)) '/' num2str(length(indstokeep)) 'neurons (NOT ENOUGHT MOTIFS)']);
OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);


%% ==================== FOR EACH EXPERIMENT COLLECT BY SYLLABLE TYPE

% OUTDAT.(analytype);
OUTDAT.bnum = OUTDAT.All_birdnum;
OUTDAT.enum= OUTDAT.All_exptnum;
OUTDAT.switch= OUTDAT.All_swnum;
OUTDAT.chanpair= OUTDAT.All_neurnum;
OUTDAT.istarg= OUTDAT.All_istarg;
OUTDAT.issame= OUTDAT.All_issame;

if strcmp(datlevel, 'unit')
    [allbnum, allenum, allswnum, allDat] = ...
        lt_neural_LFP_GrpStats(OUTDAT, analytype);
    [allbnum, allenum, allswnum, allDat_learnZTargDir] = ...
        lt_neural_LFP_GrpStats(OUTDAT, 'AllMinusBase_PitchZ_Ldir');
    
elseif strcmp(datlevel, 'switch')
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat] = ...
        lt_neural_LFP_GrpStats(OUTDAT, analytype);
    
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_learnZTargDir] = ...
        lt_neural_LFP_GrpStats(OUTDAT, 'AllMinusBase_PitchZ_Ldir');
    
end

%% =========== PLOT
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('TARG');
ylabel(analytype);
xlabel(['dat = ' datlevel]);
Y = squeeze(allDat(:,:,1,:));
t = OUTDAT.All_FRsmooth_t{1};
plot(t, Y, '-r');
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(Y(1,:))))]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
ymean = nanmean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
YLIM = ylim;
for j=1:size(Y,1)
    p = signrank(Y(j,:));
    if p<0.05
        plot(t(j), 0.9*YLIM(2), 'xk');
    end
end
% --- one for each bird
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
pcols = lt_make_plot_colors(max(allbnum), 0,0);
for k=unique(allbnum)'
    ymean = nanmean(Y(:, allbnum==k),2);
    ysem = lt_sem(Y(:, allbnum==k)');
    try
        shadedErrorBar(t, ymean, ysem, {'Color', pcols{k}}, 1);
    catch err
        plot(t, ymean, 'Color', pcols{k});
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    for j=1:size(Y,1)
        p = signrank(Y(j,allbnum==k));
        if p<0.05
            plot(t(j), 0.1*k+(0.9)*YLIM(2), 'x', 'Color', pcols{k});
        end
    end
end



% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('SAME');
ylabel(analytype);
xlabel(['dat = ' datlevel]);
Y = squeeze(allDat(:,:,2,:));
t = OUTDAT.All_FRsmooth_t{1};
plot(t, Y, '-r');
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(Y(1,:))))]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
ymean = nanmean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
YLIM = ylim;
for j=1:size(Y,1)
    p = signrank(Y(j,:));
    if p<0.05
        plot(t(j), 0.9*YLIM(2), 'xk');
    end
end
% --- one for each bird
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
pcols = lt_make_plot_colors(max(allbnum), 0,0);
for k=unique(allbnum)'
    ymean = nanmean(Y(:, allbnum==k),2);
    ysem = lt_sem(Y(:, allbnum==k)');
    try
        shadedErrorBar(t, ymean, ysem, {'Color', pcols{k}}, 1);
    catch err
        plot(t, ymean, 'Color', pcols{k});
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    for j=1:size(Y,1)
        p = signrank(Y(j,allbnum==k));
        if p<0.05
            plot(t(j), 0.1*k+(0.9)*YLIM(2), 'x', 'Color', pcols{k});
        end
    end
end

% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('DIFF');
ylabel(analytype);
xlabel(['dat = ' datlevel]);
Y = squeeze(allDat(:,:,3,:));
t = OUTDAT.All_FRsmooth_t{1};
plot(t, Y, '-r');
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(Y(1,:))))]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
ymean = nanmean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
YLIM = ylim;
for j=1:size(Y,1)
    p = signrank(Y(j,:));
    if p<0.05
        plot(t(j), 0.9*YLIM(2), 'xk');
    end
end
% --- one for each bird
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
pcols = lt_make_plot_colors(max(allbnum), 0,0);
for k=unique(allbnum)'
    ymean = nanmean(Y(:, allbnum==k),2);
    ysem = lt_sem(Y(:, allbnum==k)');
    try
        shadedErrorBar(t, ymean, ysem, {'Color', pcols{k}}, 1);
    catch err
        plot(t, ymean, 'Color', pcols{k});
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    for j=1:size(Y,1)
        p = signrank(Y(j,allbnum==k));
        if p<0.05
            plot(t(j), 0.1*k+(0.9)*YLIM(2), 'x', 'Color', pcols{k});
        end
    end
end




% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('TARG - DIFF');
ylabel(analytype);
xlabel(['dat = ' datlevel]);
Y = squeeze(allDat(:,:,1,:)) - squeeze(allDat(:,:,3,:));
t = OUTDAT.All_FRsmooth_t{1};
plot(t, Y, '-r');
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(Y(1,:))))]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
ymean = nanmean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
YLIM = ylim;
for j=1:size(Y,1)
    p = signrank(Y(j,:));
    if p<0.05
        plot(t(j), 0.9*YLIM(2), 'xk');
    end
end
% --- one for each bird
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
pcols = lt_make_plot_colors(max(allbnum), 0,0);
for k=unique(allbnum)'
    ymean = nanmean(Y(:, allbnum==k),2);
    ysem = lt_sem(Y(:, allbnum==k)');
    try
        shadedErrorBar(t, ymean, ysem, {'Color', pcols{k}}, 1);
    catch err
        plot(t, ymean, 'Color', pcols{k});
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    for j=1:size(Y,1)
        p = signrank(Y(j,allbnum==k));
        if p<0.05
            plot(t(j), 0.1*k+(0.9)*YLIM(2), 'x', 'Color', pcols{k});
        end
    end
end


% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('TARG - SAME');
ylabel(analytype);
xlabel(['dat = ' datlevel]);
Y = squeeze(allDat(:,:,1,:)) - squeeze(allDat(:,:,2,:));
t = OUTDAT.All_FRsmooth_t{1};
plot(t, Y, '-r');
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(Y(1,:))))]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
ymean = nanmean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
YLIM = ylim;
for j=1:size(Y,1)
    p = signrank(Y(j,:));
    if p<0.05
        plot(t(j), 0.9*YLIM(2), 'xk');
    end
end
% --- one for each bird
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
pcols = lt_make_plot_colors(max(allbnum), 0,0);
for k=unique(allbnum)'
    ymean = nanmean(Y(:, allbnum==k),2);
    ysem = lt_sem(Y(:, allbnum==k)');
    try
        shadedErrorBar(t, ymean, ysem, {'Color', pcols{k}}, 1);
    catch err
        plot(t, ymean, 'Color', pcols{k});
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    for j=1:size(Y,1)
        p = signrank(Y(j,allbnum==k));
        if p<0.05
            plot(t(j), 0.1*k+(0.9)*YLIM(2), 'x', 'Color', pcols{k});
        end
    end
end

% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('SAME - DIFF');
ylabel(analytype);
xlabel(['dat = ' datlevel]);
Y = squeeze(allDat(:,:,2,:)) - squeeze(allDat(:,:,3,:));
t = OUTDAT.All_FRsmooth_t{1};
plot(t, Y, '-r');
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(Y(1,:))))]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
ymean = nanmean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(t, ymean, ysem, {'Color', 'k'}, 1);
lt_plot_zeroline;
lt_plot_zeroline_vert;
YLIM = ylim;
for j=1:size(Y,1)
    p = signrank(Y(j,:));
    if p<0.05
        plot(t(j), 0.9*YLIM(2), 'xk');
    end
end
% --- one for each bird
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
pcols = lt_make_plot_colors(max(allbnum), 0,0);
for k=unique(allbnum)'
    ymean = nanmean(Y(:, allbnum==k),2);
    ysem = lt_sem(Y(:, allbnum==k)');
    try
        shadedErrorBar(t, ymean, ysem, {'Color', pcols{k}}, 1);
    catch err
        plot(t, ymean, 'Color', pcols{k});
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    for j=1:size(Y,1)
        p = signrank(Y(j,allbnum==k));
        if p<0.05
            plot(t(j), 0.1*k+(0.9)*YLIM(2), 'x', 'Color', pcols{k});
        end
    end
end

linkaxes(hsplots, 'xy');



%% ================ [EXTRACT]  SCALARS

t = OUTDAT.All_FRsmooth_t{1};
indst = t>premotorwind(1) & t<premotorwind(2);
allDat_scal = squeeze(mean(allDat(:, indst, :,:),2));

%% ============== [PLOT] SCALARS
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% =========================
indslearn = 1:size(allDat_scal,2);
% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('TARG -- SAME -- DIFF');
ylabel([analytype]);
indsrow = [1:3];

indscol = ~any(isnan(allDat_scal(indsrow, :))) & indslearn;
Y = allDat_scal(indsrow, indscol);
x = indsrow;
plot(x, Y, '-ok');

ymean = mean(Y,2);
ysem = lt_sem(Y');
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;

YLIM = ylim;
p = signrank(Y(1,:), Y(2, :));
lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);
p = signrank(Y(1,:), Y(3, :));
lt_plot_text(2, 0.9*YLIM(2), ['p(1vs3)=' num2str(p)]);
p = signrank(Y(3,:), Y(2, :));
lt_plot_text(2.5, 0.9*YLIM(2), ['p(3vs2)=' num2str(p)]);


% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('TARG -- DIFF');
ylabel([analytype]);
indsrow = [1 3];

indscol = ~any(isnan(allDat_scal(indsrow, :))) & indslearn;
Y = allDat_scal(indsrow, indscol);
x = indsrow;

plot(x, Y, '-ok');
ymean = mean(Y,2);
ysem = lt_sem(Y');
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;

YLIM = ylim;
p = signrank(Y(1,:), Y(2, :));
lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);



% ############################# SPLIT BY HIGH AND LOW LEARNING EXPERIMENTS
ylearn = squeeze(allDat_learnZTargDir(:,:, 1, :));
ysplit = median(ylearn);



% ========================= HI LEARNING
indslearn = ylearn'>ysplit;

% =====================

% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('HI LEARNING');
xlabel('TARG -- SAME -- DIFF');
ylabel([analytype]);
indsrow = [1:3];

indscol = ~any(isnan(allDat_scal(indsrow, :))) & indslearn;
Y = allDat_scal(indsrow, indscol);
x = indsrow;
plot(x, Y, '-ok');

ymean = mean(Y,2);
ysem = lt_sem(Y');
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;

YLIM = ylim;
p = signrank(Y(1,:), Y(2, :));
lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);
p = signrank(Y(1,:), Y(3, :));
lt_plot_text(2, 0.9*YLIM(2), ['p(1vs3)=' num2str(p)]);
p = signrank(Y(3,:), Y(2, :));
lt_plot_text(2.5, 0.9*YLIM(2), ['p(3vs2)=' num2str(p)]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
Ylearn = squeeze(allDat_learnZTargDir(:,:, indsrow, indscol));
plot(indsrow, Ylearn, '-ok');
ylabel('learn (z, targ dir)');
xlim([0 4]);
lt_plot_zeroline;


% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('TARG -- DIFF');
ylabel([analytype]);
indsrow = [1 3];

indscol = ~any(isnan(allDat_scal(indsrow, :))) & indslearn;

Y = allDat_scal(indsrow, indscol);
x = indsrow;

plot(x, Y, '-ok');
ymean = mean(Y,2);
ysem = lt_sem(Y');
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;

YLIM = ylim;
p = signrank(Y(1,:), Y(2, :));
lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
Ylearn = squeeze(allDat_learnZTargDir(:,:, indsrow, indscol));
plot(indsrow, Ylearn, '-ok');
xlim([0 4]);
ylabel('laern(targ dir, z)');
lt_plot_zeroline;



% ========================= HI LEARNING
indslearn = ylearn'<=ysplit;

% =====================

% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('LO LEARNING');
xlabel('TARG -- SAME -- DIFF');
ylabel([analytype]);
indsrow = [1:3];

indscol = ~any(isnan(allDat_scal(indsrow, :))) & indslearn;
Y = allDat_scal(indsrow, indscol);
x = indsrow;
plot(x, Y, '-ok');

ymean = mean(Y,2);
ysem = lt_sem(Y');
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;

YLIM = ylim;
p = signrank(Y(1,:), Y(2, :));
lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);
p = signrank(Y(1,:), Y(3, :));
lt_plot_text(2, 0.9*YLIM(2), ['p(1vs3)=' num2str(p)]);
p = signrank(Y(3,:), Y(2, :));
lt_plot_text(2.5, 0.9*YLIM(2), ['p(3vs2)=' num2str(p)]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
Ylearn = squeeze(allDat_learnZTargDir(:,:, :, indscol));
plot(indsrow, Ylearn, '-ok');
ylabel('learn (z, targ dir)');
xlim([0 4]);
lt_plot_zeroline;


% ======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('TARG -- DIFF');
ylabel([analytype]);
indsrow = [1 3];

indscol = ~any(isnan(allDat_scal(indsrow, :))) & indslearn;

Y = allDat_scal(indsrow, indscol);
x = indsrow;

plot(x, Y, '-ok');
ymean = mean(Y,2);
ysem = lt_sem(Y');
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
xlim([0 4]);
lt_plot_zeroline;

YLIM = ylim;
p = signrank(Y(1,:), Y(2, :));
lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
Ylearn = squeeze(allDat_learnZTargDir(:,:, indsrow, indscol));
plot(indsrow, Ylearn, '-ok');
xlim([0 4]);
ylabel('laern(targ dir, z)');
lt_plot_zeroline;



%% ================== CORRELATION BETWEEN LEARJING AND NEURAL CAHGNE?
lt_figure; hold on;
xlabel('learning (z, targ dir)');
ylabel([analytype]);
title('TARG (r), DIFF (k) (x jitter added)');

indsrow = [1 3];
indscol = ~any(isnan(allDat_scal(indsrow, :)));

Ytarg = allDat_scal(1, indscol);
Ydiff = allDat_scal(3, indscol);
Ylearn = squeeze(allDat_learnZTargDir(:,:, 1, indscol));

% plot(Ylearn, Ytarg, 'or');
% plot(Ylearn, Ydiff, 'ok');
for j=1:length(Ylearn)
   x = [Ylearn(j) Ylearn(j)];
   y = [Ytarg(j) Ydiff(j)];
   % -- add jitter to x
   x = x+0.02*rand-0.01;
   if y(1)>y(2)
   plot(x,y, '-', 'Color', 'm');
   else
       plot(x,y, '-', 'Color', 'b');
   end
   lt_plot(x(1), y(1), {'Color', 'r'});
%    plot(x(2), y(2), 'ok');
   lt_plot(x(2), y(2), {'Color', 'k'});
end
lt_plot_zeroline;

% ymean = mean(Y,2);
% ysem = lt_sem(Y');
% lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'r'});
% xlim([0 4]);
% lt_plot_zeroline;
% 
% YLIM = ylim;
% p = signrank(Y(1,:), Y(2, :));
% lt_plot_text(1.5, 0.9*YLIM(2), ['p(1vs2)=' num2str(p)]);
% 
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% hsplots = [hsplots; hsplot];
% Ylearn = squeeze(allDat_learnZTargDir(:,:, indsrow, indscol));
% plot(indsrow, Ylearn, '-ok');
% xlim([0 4]);
% ylabel('laern(targ dir, z)');
% lt_plot_zeroline;

