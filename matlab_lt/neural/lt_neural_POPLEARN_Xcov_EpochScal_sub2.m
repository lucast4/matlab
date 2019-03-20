%% all the stuff looking at FF split trials.

%% ========== CONCATENATE FF SPLITS
n1 = size(allDat_ffsplit_base,2) + size(allDat_epochs_ffsplit,1); % base + epochs
n2 = size(allDat_ffsplit_base,1); assert(n2 == size(allDat_epochs_ffsplit,2)); % number of ff splits
allDat_BaseEpochs_FFsplits = nan(n1, n2, 3, size(allDat_epochs_ffsplit,4));

for i=1:size(allDat_epochs_ffsplit,3)
    for ii=1:size(allDat_epochs_ffsplit,4)
   
        tmp1 = allDat_ffsplit_base(:,1, i, ii);
        tmp2 = allDat_epochs_ffsplit(:, :, i, ii);
        tmp3 = [tmp1'; tmp2];
        allDat_BaseEpochs_FFsplits(:,:, i, ii) = tmp3;
    end
end
clear allDat_ffsplit_base
clear allDat_epochs_ffsplit



%% =================== [FFSPLIT] Figure out direction of baseline bias

himinuslo = squeeze(allDat_BaseEpochs_FFsplits(1, 2, syltype, :) - ...
    allDat_BaseEpochs_FFsplits(1, 1, syltype, :)); % hi pitch minus low pitch trials

% === get mean for each expt
[indsgrp, indsgrpU, bes] = lt_tools_grp2idx({allbnum, allenum, allswnum});


lt_figure; hold on;
xlabel('expt ID');
ylabel('hi pitch minus low pitch trials (xcov) [r=extracted means]');
plot(indsgrp, himinuslo, 'ok');
lt_plot_zeroline;
title('afp bias is similar for a given expt across channel pairs');

% ===== get means
tmp = grpstats(himinuslo, indsgrp);
% expand back to origianl datapoints
assert(all(sort(indsgrp)==indsgrp));
assert(all(diff(unique(indsgrp))==1));
allDat_himinuslo_switchmeans = tmp(indsgrp);
plot(indsgrp, allDat_himinuslo_switchmeans, 'rs');


%% ================ [PREPROCESS] IF LEARNDIR DOWN, THEN FLIP DIRECTION OF FF SPLIT DATA

allDat_BaseEpochs_FFsplits_adaptivedir = nan(size(allDat_BaseEpochs_FFsplits,1), ...
    size(allDat_BaseEpochs_FFsplits,2), 1, size(allDat_BaseEpochs_FFsplits,4));
for i=1:size(allDat_BaseEpochs_FFsplits,4)
    
   learndir = allDat_targLearnDir(1,1,1, i);
   
   tmp = squeeze(allDat_BaseEpochs_FFsplits(:,:, syltype, i));
   if learndir==-1
       tmp = fliplr(tmp); % then low FF should be "adaptive" direction
   end
   
   allDat_BaseEpochs_FFsplits_adaptivedir(:,:,1, i) = tmp;
end



%% =================== [PLOT] EACH EXPT [SHOW TIMING]
% === FF SPLITS, HOW CHANGE DURING LEARNING
% === SEPARATE PLOTS SPLIT BY DIRECTION OF LEARNING

[indsgrp, indsgrpU] = lt_tools_grp2idx({allbnum, allenum, allswnum});

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(indsgrpU)
    
    indsthis = indsgrp==indsgrpU(i);
    
    bnum = unique(allbnum(indsthis));
    enum = unique(allenum(indsthis));
    sw = unique(allswnum(indsthis));
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([bname '-' ename '-sw' num2str(sw)]);
    xlabel('time (0=base)');
    ylabel('xcov (r=high ff, k=low ff, bu=mean)');
    
    yxcov = squeeze(allDat_xcov(scalwind, :, syltype, indsthis));
    if size(yxcov,1)==1
        yxcov = yxcov';
    end
    ylearn = squeeze(allDat_learn(:,:, syltype, indsthis));
    yN = squeeze(allDat_Nperbin(:,:, syltype, indsthis));
    ytime = squeeze(allDat_TimeMedian(:,:, syltype, indsthis));
    if size(ytime,1)==1
        ytime = ytime';
    end
    
    x = [0; median(ytime,2)];
    assert(length(unique(ytime(1,:)))==1, 'all neural data should be algned to same time...');
    
    % =================== xcov over bins, for each ff band
    yffsplits = allDat_BaseEpochs_FFsplits(:,:, syltype, indsthis);
    % -- for increasing ff, plot y
    nff = size(yffsplits,2);
    pcols = lt_make_plot_colors(nff, 1, [0.8 0 0]);
    for j=1:nff
       tmp = squeeze(yffsplits(:,j,:,:));
%            x = x(~isnan(tmp));
%            tmp = tmp(~isnan(tmp));
       plot(x, tmp, '-o', 'Color', pcols{j});
       if j==nff
           lt_plot_text(x(end), tmp(end), 'highFF', pcols{j});
       end
    end
    % --- plot mean over ff
    tmp = mean(mean(yffsplits,2),4);
    plot(x, tmp, '-b');
    
%     plot(x, yxcov, '-or');
%     plot(x, ylearn, '-ok');
%     lt_plot_text(x(1), ylearn(end,1), ['N/bin=' num2str(median(yN))], 'm', 8);
%     
%     

% === annotate direction of training
ldir = unique(allDat_targLearnDir(:,:, syltype, indsthis)); assert(length(ldir)==1);
lt_plot_annotation(1, ['l.dir:' num2str(ldir)], 'm');

    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
end

linkaxes(hsplots, 'xy');




% %% =================== [PLOT] EACH EXPT [SHOW TIMING]
% % === FF SPLITS, HOW CHANGE DURING LEARNING
% 
% [indsgrp, indsgrpU] = lt_tools_grp2idx({allbnum, allenum, allswnum});
% 
% figcount=1;
% subplotrows=5;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
% 
% 
% for i=1:length(indsgrpU)
%     
%     indsthis = indsgrp==indsgrpU(i);
%     
%     bnum = unique(allbnum(indsthis));
%     enum = unique(allenum(indsthis));
%     sw = unique(allswnum(indsthis));
%     bname = SwitchStruct.bird(bnum).birdname;
%     ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
%     
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title([bname '-' ename '-sw' num2str(sw)]);
%     xlabel('time (0=base)');
%     ylabel('xcov (r=high ff, k=low ff, bu=mean)');
%     
%     yxcov = squeeze(allDat_xcov(scalwind, :, syltype, indsthis));
%     if size(yxcov,1)==1
%         yxcov = yxcov';
%     end
%     ylearn = squeeze(allDat_learn(:,:, syltype, indsthis));
%     yN = squeeze(allDat_Nperbin(:,:, syltype, indsthis));
%     ytime = squeeze(allDat_TimeMedian(:,:, syltype, indsthis));
%     if size(ytime,1)==1
%         ytime = ytime';
%     end
%     
%     x = [0; median(ytime,2)];
%     assert(length(unique(ytime(1,:)))==1, 'all neural data should be algned to same time...');
%     
%     % =================== xcov over bins, for each ff band
%     yffsplits = allDat_BaseEpochs_FFsplits(:,:, syltype, indsthis);
%     % -- for increasing ff, plot y
%     nff = size(yffsplits,2);
%     pcols = lt_make_plot_colors(nff, 1, [0.8 0 0]);
%     for j=1:nff
%        tmp = squeeze(yffsplits(:,j,:,:));
%        plot(x, tmp, '-o', 'Color', pcols{j});
%        if j==nff
%            lt_plot_text(x(end), tmp(end), 'highFF', pcols{j});
%        end
%     end
%     % --- plot mean over ff
%     tmp = mean(mean(yffsplits,2),4);
%     plot(x, tmp, '-b');
%     
% %     plot(x, yxcov, '-or');
% %     plot(x, ylearn, '-ok');
% %     lt_plot_text(x(1), ylearn(end,1), ['N/bin=' num2str(median(yN))], 'm', 8);
% %     
% %     
% 
% % === annotate direction of training
% ldir = unique(allDat_targLearnDir(:,:, syltype, indsthis)); assert(length(ldir)==1);
% lt_plot_annotation(1, ['l.dir:' num2str(ldir)], 'm');
% 
%     axis tight;
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
%     
% end
% 
% linkaxes(hsplots, 'xy');


%% ================= [PLOT] SEPARATE EXPT BY LEARNING DIRECTION

figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];




% ############################## TRAIN UP
indsthis = squeeze(allDat_targLearnDir(1,1,1,:))==1;

% ========== HI AND LOW FF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('train up');
xlabel('time bin');
ylabel('xcov (hi(r), low(k) FF)');


% -- low ff
ffind = 1;
pcol = 'k';

y = squeeze(allDat_BaseEpochs_FFsplits(:, ffind, syltype, indsthis));
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);

% -- hi ff
ffind = 2;
pcol = 'r';

y = squeeze(allDat_BaseEpochs_FFsplits(:, ffind, syltype, indsthis));
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);

lt_plot_zeroline;

% -- hi ff minus low ff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('train up');
xlabel('time bin');
ylabel('xcov (hi minus lo FF)');

pcol = 'b';
y = squeeze(allDat_BaseEpochs_FFsplits(:, 2, syltype, indsthis) - allDat_BaseEpochs_FFsplits(:, 1, syltype, indsthis));
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
lt_plot_zeroline;


% ############################## TRAIN DOWN
indsthis = squeeze(allDat_targLearnDir(1,1,1,:))==-1;

% ========== HI AND LOW FF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('train down');
xlabel('time bin');
ylabel('xcov (hi(r), low(k) FF)');


% -- low ff
ffind = 1;
pcol = 'k';

y = squeeze(allDat_BaseEpochs_FFsplits(:, ffind, syltype, indsthis));
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);

% -- hi ff
ffind = 2;
pcol = 'r';

y = squeeze(allDat_BaseEpochs_FFsplits(:, ffind, syltype, indsthis));
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
lt_plot_zeroline;

% -- hi ff minus low ff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('train down');
xlabel('time bin');
ylabel('xcov (hi minus lo FF)');

pcol = 'b';
y = squeeze(allDat_BaseEpochs_FFsplits(:, 2, syltype, indsthis) - allDat_BaseEpochs_FFsplits(:, 1, syltype, indsthis));
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
lt_plot_zeroline;



%% ############################# SAME THING, BUT TAKE ALL TRIALS DURING
% TRAINING
% ############################## TRAIN UP
indsthis = squeeze(allDat_targLearnDir(1,1,1,:))==1;
ptit = 'UP';

% ========== HI AND LOW FF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title(['train ' ptit]);
xlabel('base -- WN');
ylabel('xcov (hi(r), low(k) FF)');


% -- low ff
ffind = 1;
pcol = 'k';

y1 = squeeze(allDat_BaseEpochs_FFsplits(1, ffind, syltype, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], ffind, syltype, indsthis),1));
y = [y1'; y2'];
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
xlim([0 x(end)+1]);

% -- hi ff
ffind = 2;
pcol = 'r';

y1 = squeeze(allDat_BaseEpochs_FFsplits(1, ffind, syltype, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], ffind, syltype, indsthis),1));
y = [y1'; y2'];
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
xlim([0 x(end)+1]);
lt_plot_zeroline;


% -- hi ff minus low ff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - WN');
ylabel('xcov (hi minus lo FF)');

pcol = 'b';
y1 = squeeze(allDat_BaseEpochs_FFsplits(1, 2, syltype, indsthis)) - squeeze(allDat_BaseEpochs_FFsplits(1, 1, syltype, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], 2, syltype, indsthis),1)) ...
    - squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], 1, syltype, indsthis),1));

y = [y1'; y2'];
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
xlim([0 x(end)+1]);
lt_plot_zeroline;


% ############################## TRAIN UP
indsthis = squeeze(allDat_targLearnDir(1,1,1,:))==-1;
ptit = 'DOWN';

% ========== HI AND LOW FF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title(['train ' ptit]);
xlabel('base -- WN');
ylabel('xcov (hi(r), low(k) FF)');


% -- low ff
ffind = 1;
pcol = 'k';

y1 = squeeze(allDat_BaseEpochs_FFsplits(1, ffind, syltype, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], ffind, syltype, indsthis),1));
y = [y1'; y2'];
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
xlim([0 x(end)+1]);

% -- hi ff
ffind = 2;
pcol = 'r';

y1 = squeeze(allDat_BaseEpochs_FFsplits(1, ffind, syltype, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], ffind, syltype, indsthis),1));
y = [y1'; y2'];
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
xlim([0 x(end)+1]);
lt_plot_zeroline;


% -- hi ff minus low ff
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - WN');
ylabel('xcov (hi minus lo FF)');

pcol = 'b';
y1 = squeeze(allDat_BaseEpochs_FFsplits(1, 2, syltype, indsthis)) - squeeze(allDat_BaseEpochs_FFsplits(1, 1, syltype, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], 2, syltype, indsthis),1)) ...
    - squeeze(nanmean(allDat_BaseEpochs_FFsplits([2:end], 1, syltype, indsthis),1));

y = [y1'; y2'];
x = 1:size(y,1);
plot(x, y', 'Color', pcol);
shadedErrorBar(x, nanmean(y,2), lt_sem(y'), {'Color', pcol},1);
xlim([0 x(end)+1]);
lt_plot_zeroline;




%% =============== [PLOT] CHANGE IN ADAPTIVE VS. UNADAPTIVE STATE
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

wnbins = 3;
plotraw = 0; % overlay individual pairs.
% =============== WHICH ONES TO PLOT?
% indsthis = allbnum==4;
indsthis = 1:size(allDat_BaseEpochs_FFsplits_adaptivedir,4);
% indsthis = allDat_himinuslo_switchmeans<0;

% --- experiments driving in adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) == (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

% --- experiments driving opposite adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) ~= (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

% =======================
wnbins = wnbins+1; % since bin 1 is always baseline.

if plotraw==1
    YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

% ##############################
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch')
title('r(trials with pitch in adaptive dir)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);

for i=1:size(allDat_BaseEpochs_FFsplits_adaptivedir,2)
    y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, i, 1, indsthis));
    x = 1:size(y,1);
    
    ymean = nanmean(y,2);
    ysem = lt_sem(y');
    
    if plotraw==1
       plot(x, y', '-', 'Color', [0.7 0.7 0.7]); 
    end
        
    shadedErrorBar(x, ymean, ysem, {'Color', pcolors{i}},1);
end
ylim(YLIM);

% ##############################
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch')
title('r(trials with pitch in adaptive dir)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);

for i=1:size(allDat_BaseEpochs_FFsplits_adaptivedir,2)
    y =  [squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, i, 1, indsthis))'
        squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, i, 1, indsthis),1))'];
    x = 1:size(y,1);
    
    ymean = nanmean(y,2);
    ysem = lt_sem(y');
        if plotraw==1
       plot(x, y', '-', 'Color', [0.7 0.7 0.7]); 
    end

    shadedErrorBar(x, ymean, ysem, {'Color', pcolors{i}},1);
    
    % --- p val
    p = signrank(y(2,:), y(1,:));
    lt_plot_text(x(end), ymean(end), ['p=' num2str(p) '[wn vs base]'], pcolors{i});
end
ylim(YLIM);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch) (adaptive minus nonadaptive)')
title('b(actual val); m(minsu base)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);

y = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:,2,1,indsthis) ...
    - allDat_BaseEpochs_FFsplits_adaptivedir(:,1,1,indsthis));
x = 1:size(y,1);

% ---- actual value
ymean = nanmean(y,2);
ysem = lt_sem(y');
shadedErrorBar(x, ymean, ysem, {'Color', 'b'},1);

% ---- subtract base
y = y-y(1,:);
ymean = nanmean(y,2);
ysem = lt_sem(y');
shadedErrorBar(x, ymean, ysem, {'Color', 'm'},1);
    if plotraw==1
       plot(x, y', '-', 'Color', 'm'); 
    end
ylim(YLIM);


% ================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch) (adaptive minus nonadaptive)')
title('b(actual val); m(minsu base)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);

y = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:,2,1,indsthis) ...
    - allDat_BaseEpochs_FFsplits_adaptivedir(:,1,1,indsthis));

% -- mean across bins
y = [y(1,:); nanmean(y(wnbins, :),1)];
x = 1:size(y,1);

% ---- actual value
ymean = nanmean(y,2);
ysem = lt_sem(y');
shadedErrorBar(x, ymean, ysem, {'Color', 'b'},1);
for j=1:size(y,1)
    p = signrank(y(j,:));
    lt_plot_text(j, ymean(j), ['p=' num2str(p) '[vs.0]']);
end
    p = signrank(y(2,:), y(1,:));
    lt_plot_text(x(end), ymean(end), ['p=' num2str(p) '[wn vs base]'], pcolors{i});


% ---- subtract base
y = y-y(1,:);
ymean = nanmean(y,2);
ysem = lt_sem(y');
shadedErrorBar(x, ymean, ysem, {'Color', 'm'},1);
    if plotraw==1
       plot(x, y', '-', 'Color', 'm'); 
    end
ylim(YLIM);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch)(minus base)')
title('r(trials with pitch in adaptive dir)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);

for i=1:size(allDat_BaseEpochs_FFsplits_adaptivedir,2)
    
    y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, i, 1, indsthis));
    y = y-y(1,:);
    x = 1:size(y,1);
%         if plotraw==1
%        plot(x, y', '-', 'Color', [0.7 0.7 0.7]); 
%     end

    ymean = nanmean(y,2);
    ysem = lt_sem(y');
    
    shadedErrorBar(x, ymean, ysem, {'Color', pcolors{i}},1);
end
ylim(YLIM);

% recoded below.
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% hsplots = [hsplots; hsplot];
% xlabel('base - epoch bins');
% ylabel('xcov (trials split by pitch)(minus base)')
% title('r(trials with pitch in adaptive dir)');
% nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
% pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);
% 
% for i=1:size(allDat_BaseEpochs_FFsplits_adaptivedir,2)
%     
%     y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, i, 1, indsthis));
%     y = [y(1,:); mean(y(wnbins,:),1)];
%     
%     y = y-y(1,:);
%     x = 1:size(y,1);
%     
%     ymean = mean(y,2);
%     ysem = lt_sem(y');
%     
%     shadedErrorBar(x, ymean, ysem, {'Color', pcolors{i}},1);
%     
%     
%     p = signrank(y(2,:));
%     lt_plot_text(2, ymean(2), ['p=' num2str(p) '[vs.0]']);
% end


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch)(minus base)')
title('r(trials with pitch in adaptive dir)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);


y = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, :, 1, indsthis),1) - ...
    nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(1, :, 1, indsthis),1));

for i=1:size(allDat_BaseEpochs_FFsplits_adaptivedir,2)
    
%     ymean = 
    x = 1:size(y,1);
    ymean = [0 nanmean(y(i,:),2)];
    ysem = [0 lt_sem(y(i,:)')];
    if plotraw==1
       plot(2+i*0.1, y(i,:)', 'x', 'Color', pcolors{i}); 
    end

    shadedErrorBar(x, ymean, ysem, {'Color', pcolors{i}},1);
    
    p = signrank(y(i,:));
    lt_plot_text(2, ymean(2), ['p=' num2str(p) '[vs.0]']);
end
    p = signrank(y(2,:), y(1,:));
    lt_plot_text(x(end), ymean(end), ['p=' num2str(p) '[adaptive vs. non]'], pcolors{i});
ylim(YLIM);

    
    
%% ========= CLEANER PLOTS [SUMMARY, BAR PLOTS]
figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

wnbins = 3;
plotraw = 1; % overlay individual pairs.

% =============== WHICH ONES TO PLOT?
% indsthis = allbnum==4
indsthis = 1:size(allDat_BaseEpochs_FFsplits_adaptivedir,4);
% indsthis = allDat_himinuslo_switchmeans<0;

% --- experiments driving in adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) == (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

% --- experiments driving opposite adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) ~= (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

% =======================
wnbins = wnbins+1; % since bin 1 is always baseline.
if plotraw==1
YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

% =========================================== COLLECT DATA
y1 = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, 1, 1, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, 1, 1, indsthis),1));
Ynonadapt = [y1 y2];

y1 = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, 2, 1, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, 2, 1, indsthis),1));
Yadapt = [y1 y2];

Ybnum = allbnum(indsthis);

% ############################## PLOTS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base -- WN');
title('nonadaptive trials');
pcol = 'k';
x = [1 2];
Y = Ynonadapt;
if plotraw==1
plot(x, Y', '-', 'Color', pcol);
end
ymean = nanmean(Y,1);
ysem = lt_sem(Y);
lt_plot_bar(x, ymean, {'Errors', ysem});
xlim([0 3]);
ylim(YLIM);


% ############################## PLOTS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base -- WN');
title('adaptive trials');
pcol = 'r';
x = [1 2];
Y = Yadapt;
if plotraw==1
plot(x, Y', '-', 'Color', pcol);
end
ymean = nanmean(Y,1);
ysem = lt_sem(Y);
lt_plot_bar(x, ymean, {'Errors', ysem});
xlim([0 3]);
ylim(YLIM);



% ############################## PLOTS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('BASE(nonad-adapt) -- WN(nonad-adapt)');

% -- base
Y = [Ynonadapt(:,1) Yadapt(:,1)];
x = [1 2];
if plotraw==1
plot(x, Y', '-', 'Color', pcol);
end
ymean = nanmean(Y,1);
ysem = lt_sem(Y);
lt_plot_bar(x, ymean, {'Errors', ysem});

% -- wn
Y = [Ynonadapt(:,2) Yadapt(:,2)];
x = [3 4];
if plotraw==1
plot(x, Y', '-', 'Color', pcol);
end
ymean = nanmean(Y,1);
ysem = lt_sem(Y);
lt_plot_bar(x, ymean, {'Errors', ysem});

% --
xlim([0 5]);
ylim(YLIM);



% ############################## PLOTS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('nonad == adapt)');
ylabel('change from baseline');
pcol = 'b';

% -- wn
Y = [Ynonadapt(:,2) Yadapt(:,2)] - [Ynonadapt(:,1) Yadapt(:,1)];
x = [1 2];
if plotraw==1
plot(x, Y', '-', 'Color', pcol);
end
ymean = nanmean(Y,1);
ysem = lt_sem(Y);
lt_plot_bar(x, ymean, {'Errors', ysem});

p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'vs', 1);
for i=1:size(Y,2)
    p = signrank(Y(:,i));
    if p<0.15
        lt_plot_text(i, max(Y(:,i)), ['p=' num2str(p)]);
    end
end

% -- overlay each bird
for i=1:max(Ybnum)
   y = Y(Ybnum==i,:);
   if isempty(y)
       continue
   end
      
    ymean = nanmean(y,1);
    ysem = lt_sem(y);
    lt_plot(x, ymean, {'Errors', ysem, 'LineStyle', '-'});

    p1 = signrank(y(:,1), y(:,2));
%     lt_plot_pvalue(p, 'vs', 1);
    for ii=2
        p = signrank(y(:,ii));
            lt_plot_text(x(end), ymean(end), ['p(vs)= ' num2str(p1), 'p(2)=' num2str(p) '[bnum' num2str(i)]);
    end
end
% --
xlim([0 3]);
ylim(YLIM);



%% ================== PLOT ALL EXPERIMENTS 
% =========== SCATTER, PLOT FF AND BIAS
indsthis = 1:size(allDat_BaseEpochs_FFsplits_adaptivedir,4);

ybase = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, :, 1, indsthis));
ywn = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir([2:end], :, 1, indsthis),1));
% ywn = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir([end], :, 1, indsthis),1));
ywn = ywn'-ybase';

learnAway = squeeze(nanmean(allDat_learn_Away(:,1,syltype, indsthis),1));
learnRevert = squeeze(nanmean(allDat_learn_Revert(:,1,syltype, indsthis),1));
% learnAway = squeeze(nanmean(allDat_learn_Away(end,1,syltype, indsthis),1));
% learnRevert = squeeze(nanmean(allDat_learn_Revert(end,1,syltype, indsthis),1));
learnwn = [learnRevert learnAway];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('away vs. revert trials [mean over WN bins]');
xlabel('learn(z, away from base');
ylabel('xcov, minus base');
for j=1:size(ywn,1)
    x = learnwn(j,:);
    y = ywn(j,:);
    plot(x, y, '-ok');
end


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('base - epoch bins');
ylabel('xcov (trials split by pitch)(minus base)')
title('r(trials with pitch in adaptive dir)');
nff = size(allDat_BaseEpochs_FFsplits_adaptivedir,2);
pcolors = lt_make_plot_colors(nff, 1, [0.9 0.1 0.1]);

for i=1:size(allDat_BaseEpochs_FFsplits_adaptivedir,2)
    
    y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, i, 1, indsthis));
    y = [y(1,:); nanmean(y(wnbins,:),1)];
    
    y = y-y(1,:);
    x = 1:size(y,1);
    
    ymean = nanmean(y,2);
    ysem = lt_sem(y');
    
    shadedErrorBar(x, ymean, ysem, {'Color', pcolors{i}},1);
    
    
    p = signrank(y(2,:));
    lt_plot_text(2, ymean(2), ['p=' num2str(p) '[vs.0]']);
end