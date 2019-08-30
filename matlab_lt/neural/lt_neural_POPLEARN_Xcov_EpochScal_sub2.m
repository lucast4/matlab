%% all the stuff looking at FF split trials.

%% ========== CONCATENATE FF SPLITS
n1 = size(allDat_ffsplit_base,2) + size(allDat_epochs_ffsplit,1); % base + epochs
n2 = size(allDat_ffsplit_base,1); assert(n2 == size(allDat_epochs_ffsplit,2)); % number of ff splits
allDat_BaseEpochs_FFsplits = nan(n1, n2, 3, size(allDat_epochs_ffsplit,4));

for i=1:size(allDat_epochs_ffsplit,3) % syl types
    for ii=1:size(allDat_epochs_ffsplit,4) % cases
   
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
% himinuslo = squeeze(allDat_BaseEpochs_FFsplits(1, 2, 2, :) - ...
%     allDat_BaseEpochs_FFsplits(1, 1, 2, :)); % hi pitch minus low pitch trials

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

%% ===================== 
if useAd_Nonad_Average_forBaseline ==1
    
    tmp = mean(allDat_BaseEpochs_FFsplits_adaptivedir(1, :, 1, :), 2);
    
    allDat_BaseEpochs_FFsplits_adaptivedir(1, 1, 1, :) = tmp;
    allDat_BaseEpochs_FFsplits_adaptivedir(1, 2, 1, :) = tmp;
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
    if all(isnan(yxcov(:)))
        continue
    end
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
    if all(isnan(yxcov(:)))
        continue
    end
    yN = squeeze(allDat_Nperbin(:,:, syltype, indsthis));
    ytime = squeeze(allDat_TimeMedian(:,:, syltype, indsthis));
    if size(ytime,1)==1
        ytime = ytime';
    end
    
    x = [0; median(ytime,2)];
    assert(length(unique(ytime(1,:)))==1, 'all neural data should be algned to same time...');
    
    % =================== xcov over bins, for each ff band
    yffsplits = allDat_BaseEpochs_FFsplits_adaptivedir(:,:, 1, indsthis);
    % -- for increasing ff, plot y
    nff = size(yffsplits,2);
    pcols = lt_make_plot_colors(nff, 1, [0.8 0 0]);
    for j=1:nff
       tmp = squeeze(yffsplits(:,j,:,:));
%            x = x(~isnan(tmp));
%            tmp = tmp(~isnan(tmp));
       plot(x, tmp, '-o', 'Color', pcols{j});
       if j==nff
           lt_plot_text(x(end), tmp(end), 'adaptiveFFdir', pcols{j});
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




%% [PUBLICATION] =============== [PLOT] CHANGE IN ADAPTIVE VS. UNADAPTIVE STATE

% wnbins = size(allDat_BaseEpochs_FFsplits_adaptivedir,1)-1; % use the last bin.
wnbins = [3 4];
plotraw = 1; % overlay individual pairs.
% =======================
wnbins = wnbins+1; % since bin 1 is always baseline.



% ####################################### WHICH ONES TO PLOT?
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% indsthis = allbnum==4;
indsthis = 1:size(allDat_BaseEpochs_FFsplits_adaptivedir,4);
% indsthis = allDat_himinuslo_switchmeans<0;

% --- experiments driving in adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) == (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

% --- experiments driving opposite adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) ~= (squeeze(allDat_targLearnDir(1,1,1,:))==1); 


if plotraw==1
    YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

lt_neural_POPLEARN_XCov_EpochScal_sub21;


% ####################################### WHICH ONES TO PLOT?
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indsthis = (allDat_himinuslo_switchmeans>0) == (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

% --- experiments driving opposite adaptive direction
% indsthis = (allDat_himinuslo_switchmeans>0) ~= (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

if plotraw==1
    YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

lt_neural_POPLEARN_XCov_EpochScal_sub21;



% ####################################### WHICH ONES TO PLOT?
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% --- experiments driving opposite adaptive direction
indsthis = (allDat_himinuslo_switchmeans>0) ~= (squeeze(allDat_targLearnDir(1,1,1,:))==1); 

if plotraw==1
    YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

lt_neural_POPLEARN_XCov_EpochScal_sub21;



% ####################################### WHICH ONES TO PLOT?
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


indsthis = squeeze(allDat_targLearnDir(1,1,1,:))==1;

if plotraw==1
    YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

lt_neural_POPLEARN_XCov_EpochScal_sub21;
lt_subtitle('driving pitch up');


% ####################################### WHICH ONES TO PLOT?
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

indsthis = squeeze(allDat_targLearnDir(1,1,1,:))==-1;

if plotraw==1
    YLIM = [-0.2 0.5];
else
    YLIM = [-0.05 0.15];
end

lt_neural_POPLEARN_XCov_EpochScal_sub21;
lt_subtitle('driving pitch down');


%% ================== PLOT ALL EXPERIMENTS 
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


wnbins = wnbins; % which bin to plot? (i.e. wn bins, 1,2, ... , +1)

% =========== SCATTER, PLOT FF AND BIAS
indsthis = 1:size(allDat_BaseEpochs_FFsplits_adaptivedir,4);


% #################### DEVIATION FRO BASE, LEARNING AND NEURAL
ybase = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, :, 1, indsthis));
% ywn = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir([2:end], :, 1, indsthis),1));
ywn = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, :, 1, indsthis),1));
ywn = ywn'-ybase';
% ywn = ywn';

% learnAway = squeeze(nanmean(allDat_learn_Away(:,1,syltype, indsthis),1));
% learnRevert = squeeze(nanmean(allDat_learn_Revert(:,1,syltype, indsthis),1));
learnAway = squeeze(nanmean(allDat_learn_Away(wnbins-1,1,syltype, indsthis),1));
learnRevert = squeeze(nanmean(allDat_learn_Revert(wnbins-1,1,syltype, indsthis),1));
learnwn = [learnRevert learnAway];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title(['away vs. revert trials, wnbins: ' num2str(wnbins-1)]);
xlabel('learn(z, away from base');
ylabel('xcov, minus base (this set of trials)');
for j=1:size(ywn,1)
    x = learnwn(j,:);
    y = ywn(j,:);
%     plot(x, y, '-ok');
    plot(x(1), y(1), 'xk');
    plot(x(2), y(2), 'xr');
    if y(2)>y(1)
        plot(x, y, '-m');    
    else
        plot(x, y, '-b');
    end
end
lt_plot(nanmean(learnwn(:,1),1), nanmean(ywn(:,1),1), ...
    {'Xerrors', lt_sem(learnwn(:,1)), 'Errors', lt_sem(ywn(:,1)), 'Color', 'k'});
lt_plot(nanmean(learnwn(:,2),1), nanmean(ywn(:,2),1), ...
    {'Xerrors', lt_sem(learnwn(:,1)), 'Errors', lt_sem(ywn(:,1)), 'Color', 'r'});
% lt_plot(mean(learnwn,1), mean(ywn,1), {'Xerrors', lt_sem(learnwn), 'Errors', lt_sem(ywn)});
p = signrank(ywn(:,2), ywn(:,1));
lt_plot_pvalue(p, 'xcov, away vs. toward base');
xlim([-1.5 2.5]); ylim([-0.35 0.6]);
lt_plot_zeroline;
lt_plot_zeroline_vert;


% #################### DEVIATION FRO BASE, LEARNING AND NEURAL
ybase = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, :, 1, indsthis));
% ywn = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir([2:end], :, 1, indsthis),1));
ywn = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, :, 1, indsthis),1));
% ywn = ywn'-ybase';
ywn = ywn';

% learnAway = squeeze(nanmean(allDat_learn_Away(:,1,syltype, indsthis),1));
% learnRevert = squeeze(nanmean(allDat_learn_Revert(:,1,syltype, indsthis),1));
learnAway = squeeze(nanmean(allDat_learn_Away(wnbins-1,1,syltype, indsthis),1));
learnRevert = squeeze(nanmean(allDat_learn_Revert(wnbins-1,1,syltype, indsthis),1));
learnwn = [learnRevert learnAway];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title(['away vs. revert trials, wnbins: ' num2str(wnbins-1)]);
xlabel('learn(z, away from base');
ylabel('xcov, not minus base (this set of trials)');
for j=1:size(ywn,1)
    x = learnwn(j,:);
    y = ywn(j,:);
%     plot(x, y, '-ok');
    plot(x(1), y(1), 'xk');
    plot(x(2), y(2), 'xr');
    if y(2)>y(1)
        plot(x, y, '-m');    
    else
        plot(x, y, '-b');
    end
end
lt_plot(nanmean(learnwn(:,1),1), nanmean(ywn(:,1),1), ...
    {'Xerrors', lt_sem(learnwn(:,1)), 'Errors', lt_sem(ywn(:,1)), 'Color', 'k'});
lt_plot(nanmean(learnwn(:,2),1), nanmean(ywn(:,2),1), ...
    {'Xerrors', lt_sem(learnwn(:,1)), 'Errors', lt_sem(ywn(:,1)), 'Color', 'r'});
% lt_plot(mean(learnwn,1), mean(ywn,1), {'Xerrors', lt_sem(learnwn), 'Errors', lt_sem(ywn)});
p = signrank(ywn(:,2), ywn(:,1));
lt_plot_pvalue(p, 'xcov, away vs. toward base');
xlim([-1.5 2.5]); ylim([-0.35 0.6]);
lt_plot_zeroline;
lt_plot_zeroline_vert;



% #############
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



% ================ mixed effects model [CHANGE FOR ADAPTIVE DIR
% SIGNIFICANT?]
y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, 2, 1, indsthis));
y = [y(1,:); nanmean(y(wnbins,:),1)];
y = y-y(1,:);
y = y(2, :)';

allcase = lt_tools_grp2idx({allbnum, allenum, allswnum});
allcase = categorical(allcase);
allbnum_cat = categorical(allbnum);
allenum_cat = categorical(allenum);
allswnum_cat = categorical(allswnum);

dattmp = table(y, allcase, allbnum_cat, allenum_cat, allswnum_cat);
% formula = 'y ~ 1';
formula = 'y ~ 1 + (1|allcase)';
% formula = 'y ~ 1 + (1|allbnum_cat) + (1|allenum_cat)';
lme = fitlme(dattmp, formula)



% ================ mixed effects model [CHANGE FOR NONADAPTIVE DIR
% SIGNIFICANT?]
y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, 1, 1, indsthis));
y = [y(1,:); nanmean(y(wnbins,:),1)];
y = y-y(1,:);
y = y(2, :)';

allcase = lt_tools_grp2idx({allbnum, allenum, allswnum});
allcase = categorical(allcase);
allbnum_cat = categorical(allbnum);
allenum_cat = categorical(allenum);
allswnum_cat = categorical(allswnum);

dattmp = table(y, allcase, allbnum_cat, allenum_cat, allswnum_cat);
% formula = 'y ~ 1';
formula = 'y ~ 1 + (1|allcase)';
% formula = 'y ~ 1 + (1|allbnum_cat) + (1|allenum_cat)';
lme = fitlme(dattmp, formula)



% ================ mixed effects model [ADAPTIVE DIR CHANGE MINUS
% NON-ADAPTIVE DIR CHANGE?]
y =  squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, 2, 1, indsthis)) - ...
    squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(:, 1, 1, indsthis));
y = [y(1,:); nanmean(y(wnbins,:),1)];
y = y-y(1,:);
y = y(2, :)';


allcase = lt_tools_grp2idx({allbnum, allenum, allswnum});
allcase = categorical(allcase);
allbnum_cat = categorical(allbnum);
allenum_cat = categorical(allenum);
allswnum_cat = categorical(allswnum);

dattmp = table(y, allcase, allbnum_cat, allenum_cat, allswnum_cat);
% formula = 'y ~ 1';
formula = 'y ~ 1 + (1|allcase)';
% formula = 'y ~ 1 + (1|allbnum_cat) + (1|allenum_cat)';
lme = fitlme(dattmp, formula)


%% ############### SPLIT EXPERIMENTS BY STRONG AND WEAK LEARNING
learnz = squeeze(nanmean(allDat_learn_pu69combined(wnbins-1, 1, syltype, :), 1));
xcovchange = squeeze(nanmean(allDat_xcov_pu69combined(scalwind, wnbins-1, syltype, :), 2));
xcovchange_Adapt = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, 2, 1, :), 1) - ...
    allDat_BaseEpochs_FFsplits_adaptivedir(1, 2, 1, :));
xcovchange_Nonad = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, 1, 1, :), 1) - ...
    allDat_BaseEpochs_FFsplits_adaptivedir(1, 1, 1, :));

% ================
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('learn');
ylabel('xcov change (adaptive)');
x = learnz;
y = xcovchange_Adapt;
plot(x, y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('learn');
ylabel('xcov change (non)');
x = learnz;
y = xcovchange_Nonad;
plot(x, y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('learn');
ylabel('xcov change (ad(red) - non(bk))');
title('slight x jitter added');
x = learnz + 0.012*(rand(size(learnz))-0.5);
y = [xcovchange_Nonad xcovchange_Adapt];
for i=1:length(x)
   plot(x(i), y(i,1), 'ok');
   plot(x(i), y(i,2), 'or');
   if y(i,2)>y(i,1)
       plot([x(i) x(i)], y(i,:), '-r');
   else
       plot([x(i) x(i)], y(i,:), '-', 'Color', [0.7 0.7 0.7]);
   end
   
   lab = [allbnum(i) allenum(i) allswnum(i)];
   lt_plot_text(x(i)+0.01, mean(y(i,:)), num2str(lab), 'm', 8);
   
end
lt_plot_makesquare_plot45line(gca, 'r');

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('learn');
ylabel('xcov change (ad(red) - non(bk))');
title('slight x jitter added');
x = learnz + 0.012*(rand(size(learnz))-0.5);
y = [xcovchange_Nonad xcovchange_Adapt];
for i=1:length(x)
   plot(x(i), y(i,1), 'ok');
   plot(x(i), y(i,2), 'or');
   if y(i,2)>y(i,1)
       plot([x(i) x(i)], y(i,:), '-r');
   else
       plot([x(i) x(i)], y(i,:), '-', 'Color', [0.7 0.7 0.7]);
   end
end
lt_plot_makesquare_plot45line(gca, 'r');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('learn');
ylabel('xcov change (adaptive - nonadap)');
x = learnz;
y = xcovchange_Adapt - xcovchange_Nonad;

plot(x, y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('learn');
ylabel('xcovchange');
x = learnz;
y = xcovchange;
plot(x, y, 'ok');
lt_plot_makesquare_plot45line(gca, 'r');



% ##############################
lt_figure; hold on;

lt_subplot(2,2,1); hold on;
[indsgrp, indsgrpU] = lt_tools_grp2idx({allbnum, allenum, allswnum});
xlabel('learn(z)');
title('each expt mean');
ylabel('xcov, change from base (b=mean, red=adap; k=nonadapt');
[learnz_grp_mean, learnz_grp_sem] = grpstats(learnz, indsgrp, {'mean', 'sem'});
[xcovchange_grp_mean, xcovchange_grp_sem] = grpstats(xcovchange, indsgrp, {'mean', 'sem'});
[xcovchangeAd_grp_mean, xcovchangeAd_grp_sem] = grpstats(xcovchange_Adapt, indsgrp, {'mean', 'sem'});
[xcovchangeNon_grp_mean, xcovchangeNon_grp_sem] = grpstats(xcovchange_Nonad, indsgrp, {'mean', 'sem'});

% ----
lt_plot(learnz_grp_mean-0.01, xcovchange_grp_mean, {'Errors', xcovchange_grp_sem, 'Color', 'b'})

lt_plot(learnz_grp_mean, xcovchangeAd_grp_mean, {'Errors', xcovchangeAd_grp_sem, 'Color', 'r'})
lt_plot(learnz_grp_mean+0.01, xcovchangeNon_grp_mean, {'Errors', xcovchangeNon_grp_sem, 'Color', 'k'})
for i=1:length(learnz_grp_mean)
    plot([learnz_grp_mean(i) learnz_grp_mean(i)+0.01], ...
    [xcovchangeAd_grp_mean(i) xcovchangeNon_grp_mean(i)], '-k');
end

lt_plot_zeroline;
lt_plot_zeroline_vert;

% ================
lt_subplot(2,2,2); hold on;
[indsgrp, indsgrpU] = lt_tools_grp2idx({allbnum, allenum, allswnum});
xlabel('learn(z)');
title('each expt mean');
ylabel('xcov, change from base (ADAPT - NON)');

[xcovchangeAdMinNon_grp_mean, xcovchangeAdMinNon_grp_sem] = grpstats(...
    xcovchange_Adapt - xcovchange_Nonad, indsgrp, {'mean', 'sem'});

% ----
lt_plot(learnz_grp_mean-0.01, xcovchangeAdMinNon_grp_mean, {'Errors', ...
    xcovchangeAdMinNon_grp_sem, 'Color', 'b'})

% lt_plot(learnz_grp_mean, xcovchangeAd_grp_mean, {'Errors', xcovchangeAd_grp_sem, 'Color', 'r'})
% lt_plot(learnz_grp_mean+0.01, xcovchangeNon_grp_mean, {'Errors', xcovchangeNon_grp_sem, 'Color', 'k'})
% for i=1:length(learnz_grp_mean)
%     plot([learnz_grp_mean(i) learnz_grp_mean(i)+0.01], ...
%     [xcovchangeAd_grp_mean(i) xcovchangeNon_grp_mean(i)], '-k');
% end

lt_plot_zeroline;
lt_plot_zeroline_vert;



% ########################### SPLIT EXPERIEMNTS BY HI AND LOW LEARNING
learnthresh = 0.2;

% --
lt_subplot(2,2,3:4); hold on;
xlabel('WEAK LEARN -- STRONG LEARNING -- WEAK/STRONG - TRIALS NOT SPLIT');
title(['learn threshold = ' num2str(learnthresh)]);
ylabel('xcov change (r=adapt; k=nonadapt)');

Ytmp = cell(1,2);

% --------- 
indsthis = learnz<learnthresh;
y = [xcovchange_Nonad(indsthis) xcovchange_Adapt(indsthis)];
Ytmp{1} = y;

indsthis = learnz>learnthresh;
y = [xcovchange_Nonad(indsthis) xcovchange_Adapt(indsthis)];
Ytmp{2} = y;

x = [1 2];
plot(x, Ytmp{1}, '-k');
plot(x(1), Ytmp{1}(:,1), 'ok');
plot(x(2), Ytmp{1}(:,2), 'or');
lt_plot(x+0.2, mean(Ytmp{1},1), {'Errors', lt_sem(Ytmp{1})});

x = [4 5];
plot(x, Ytmp{2}, '-k');
plot(x(1), Ytmp{2}(:,1), 'ok');
plot(x(2), Ytmp{2}(:,2), 'or');
lt_plot(x+0.2, mean(Ytmp{2},1), {'Errors', lt_sem(Ytmp{2})});


% ==========
x = [7 8];
y = cellfun(@(x)diff(x, 1, 2), Ytmp, 'UniformOutput', 0);
lt_plot_MultDist(y, x, 1, 'k', 0)

p = signrank(y{1});
lt_plot_text(x(1), 0.7, ['p=' num2str(p)], 'm', 9)
p = signrank(y{2});
lt_plot_text(x(2), 0.7, ['p=' num2str(p)], 'm', 9)
[p] = ranksum(y{1}, y{2});
lt_plot_pvalue(p);


% ========================== ALL TRIALS, WITHOUT FF SPLIT
Ytmp = cell(1,2);

% --------- 
indsthis = learnz<learnthresh;
y = xcovchange(indsthis);
Ytmp{1} = y;

indsthis = learnz>learnthresh;
y = xcovchange(indsthis);
Ytmp{2} = y;

% ----
x = [10 11];
lt_plot_MultDist(Ytmp, x, 1, 'k', 0)

p = signrank(Ytmp{1});
lt_plot_text(x(1), 0.7, ['p=' num2str(p)], 'm', 9)
p = signrank(Ytmp{2});
lt_plot_text(x(2), 0.7, ['p=' num2str(p)], 'm', 9)
[p] = ranksum(Ytmp{1}, Ytmp{2});
lt_plot_pvalue(p);


% ---
xlim([0 12]);
lt_plot_zeroline;



