function lt_neural_POPLEARN_Xcov_EpochScal(OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, dattype, scalwind, syltype, casestokeep, ...
    mintraindur, mintotaltrain, adhoc_replacelearnwithWind2)
%% lt 3/5/19 - divides up training into epochs -- here PLOTS

%% ====== [PREPROCESS] only plot good experiments

if onlygoodexpt==1
    
    % ===== filter outstruct
    [OUTSTRUCT_XCOV] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
    %
    %     % ===== filter outstruct
    %     [OUTSTRUCT indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, ...
    %         SwitchStruct, 'xcov_spikes');
    
end


%% ====== [PREPROCESS] average, for each experiment

%% group based on syl type

if strcmp(dattype, 'switch')
    fieldtoget = 'LearnTargDir_Z_Epochs';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_learn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'Xcovscal_window_Epochs';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_xcov] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'Xcovscal_window_BaseWN';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_xcov_basewn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'NperBin_Epochs';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_Nperbin] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'TimeBins_Base_Epochs';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_TimeMedian] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
elseif strcmp(dattype, 'chan')
    fieldtoget = 'LearnTargDir_Z_Epochs';
    [allbnum, allenum, allswnum, allDat_learn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'Xcovscal_window_Epochs';
    [allbnum, allenum, allswnum, allDat_xcov] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'Xcovscal_window_BaseWN';
    [allbnum, allenum, allswnum, allDat_xcov_basewn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'NperBin_Epochs';
    [allbnum, allenum, allswnum, allDat_Nperbin] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'TimeBins_Base_Epochs';
    [allbnum, allenum, allswnum, allDat_TimeMedian] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
end



% ============= GET DIFFERENCE FROM BASELINE FOR ALL XCOV SCALARS
Nepoch = size(allDat_xcov,2);
for j=1:size(allDat_xcov,1)
    allDat_xcov(j, :,:,:) = allDat_xcov(j,:,:,:) - repmat(allDat_xcov_basewn(j,1,:,:), 1, Nepoch,1,1);
end

% =========== GET TIME AS DIFFERENCE FROM BASELINE
allDat_TimeMedian = allDat_TimeMedian(:,2:end, :,:) - ...
    repmat(allDat_TimeMedian(1,1,:,:), 1, size(allDat_TimeMedian,2)-1, 1, 1);



%% ==================== [AD HOC] MODIFY SO THAT IS COMPARING WINDOW 1 VS. WINDOW 2
% DO THIS BY REPLACING LEARNIG MATRIX WITH THE SECOND PEAK



% neural1 = squeeze(allDat_xcov(1, :, syltype,:));
% neural2 = squeeze(allDat_xcov(2, :, syltype,:));

if adhoc_replacelearnwithWind2==1
    allDat_learn = reshape(allDat_xcov(2, :, :, :), size(allDat_learn));
end
%% =================== [PLOT] EACH EXPT,

[indsgrp, indsgrpU] = lt_tools_grp2idx({allbnum, allenum, allswnum});

figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


pcol = 'k';
for i=1:length(indsgrpU)
    
    indsthis = indsgrp==indsgrpU(i);
    
    bnum = unique(allbnum(indsthis));
    enum = unique(allbnum(indsthis));
    sw = unique(allbnum(indsthis));
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([bname '-' ename '-sw' num2str(sw)]);
    ylabel('learning (z rel base)')
    xlabel('xcov minus base');
    
    yxcov = squeeze(allDat_xcov(scalwind, :, syltype, indsthis));
    if size(yxcov,1)==1
        yxcov = yxcov';
    end
    ylearn = squeeze(allDat_learn(:,:, syltype, indsthis));
    
    %     disp(size(ylearn,1));
    for j=1:size(ylearn,2)
        plot(yxcov(:, j), ylearn(:,j), '-o', 'Color', pcol);
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

linkaxes(hsplots, 'xy');

%% =================== [PLOT] EACH EXPT,

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
    xlabel('bin');
    ylabel('dev from base (rd=neural)(k=learn)');
    
    yxcov = squeeze(allDat_xcov(scalwind, :, syltype, indsthis));
    if size(yxcov,1)==1
        yxcov = yxcov';
    end
    ylearn = squeeze(allDat_learn(:,:, syltype, indsthis));
    yN = squeeze(allDat_Nperbin(:,:, syltype, indsthis));
    
    
    x = 1:size(yxcov,1);
    plot(x, yxcov, '-or');
    plot(x, ylearn, '-ok');
    lt_plot_text(x(1), ylearn(end,1), ['N/bin=' num2str(median(yN))], 'm', 8);
    
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
end

linkaxes(hsplots, 'xy');

%% =================== [PLOT] EACH EXPT [SHOW TIMING]

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
    xlabel('time (bin median)');
    ylabel('dev from base (rd=neural)(k=learn)');
    
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
    
    x = median(ytime,2);
    assert(length(unique(ytime(1,:)))==1, 'all neural data should be algned to same time...');
    plot(x, yxcov, '-or');
    plot(x, ylearn, '-ok');
    lt_plot_text(x(1), ylearn(end,1), ['N/bin=' num2str(median(yN))], 'm', 8);
    
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
end

linkaxes(hsplots, 'xy');


%% ================ plot distribution of expt durations

lt_figure; hold on;
exptdur = squeeze(allDat_TimeMedian(1,end, 1, :) - allDat_TimeMedian(1,1, 1, :));
plot(exptdur, 'ok');
xlabel('case num');
ylabel('expt dur (hours, last bin minus first');
title('red = min train dur');
line(xlim, [mintraindur mintraindur], 'Color', 'r');

% =============== ONLY KEEP IF DURATION TRAINING LONG ENOUGH.
% LAST BIN MIN US FIRST BIN (OF TRAINING) - THEREFORE NEE DTO HAVE DATA
% THROUGHOUT TRAINING.

indstokeep = exptdur>mintraindur;
allDat_learn = allDat_learn(:,:,:, indstokeep);
allDat_xcov = allDat_xcov(:,:,:, indstokeep);
allDat_xcov_basewn = allDat_xcov_basewn(:,:,:, indstokeep);
allDat_Nperbin = allDat_Nperbin(:,:,:, indstokeep);
allDat_TimeMedian = allDat_TimeMedian(:,:,:, indstokeep);

allbnum = allbnum(indstokeep);
allenum = allenum(indstokeep);
allswnum = allswnum(indstokeep);

%% ============= REMOVE CASES THAT HAVE TOO SHORT DURATION OF EXPEIRMENT
% LAST BIN M INUS BASE - DOESN'T CARE WHETHER HAS DATA IN BETWEEM

recdur = squeeze(allDat_TimeMedian(1,end, 1, :));
indstokeep = recdur>mintotaltrain;
allDat_learn = allDat_learn(:,:,:, indstokeep);
allDat_xcov = allDat_xcov(:,:,:, indstokeep);
allDat_xcov_basewn = allDat_xcov_basewn(:,:,:, indstokeep);
allDat_Nperbin = allDat_Nperbin(:,:,:, indstokeep);
allDat_TimeMedian = allDat_TimeMedian(:,:,:, indstokeep);

allbnum = allbnum(indstokeep);
allenum = allenum(indstokeep);
allswnum = allswnum(indstokeep);



%% ============= [RESTRICT ANALYSIS TO ONLY SPECIFIC CASES]

indsgoodlearn = squeeze(allDat_learn(end, 1, 1, :))>0.1;
indsgoodneural = squeeze(allDat_xcov(scalwind, end, 1, :))>0;

if strcmp(casestokeep, 'goodlearn')
    indstokeep =   indsgoodlearn;
elseif strcmp(casestokeep, 'goodneural')
    indstokeep =   indsgoodneural;
elseif strcmp(casestokeep, 'goodlearnneural')
    indstokeep =   indsgoodlearn & indsgoodneural;
elseif strcmp(casestokeep, 'badlearn')
    indstokeep =  ~indsgoodlearn;
elseif strcmp(casestokeep, 'all')
    % keep all
    indstokeep = logical(ones(size(allDat_learn,4),1));
else
    adfasdfasdff;
end

allDat_learn = allDat_learn(:,:,:, indstokeep);
allDat_xcov = allDat_xcov(:,:,:, indstokeep);
allDat_xcov_basewn = allDat_xcov_basewn(:,:,:, indstokeep);
allDat_Nperbin = allDat_Nperbin(:,:,:, indstokeep);
allDat_TimeMedian = allDat_TimeMedian(:,:,:, indstokeep);

allbnum = allbnum(indstokeep);
allenum = allenum(indstokeep);
allswnum = allswnum(indstokeep);


%% ============ COMPARE LEARNING AND XCOV (ENDPOINTS)
lt_figure; hold on;
Nbins = size(allDat_learn,1);
pcols = lt_make_plot_colors(max(allbnum),0,0);
for i=1:Nbins
    
    learn = squeeze(allDat_learn(i, 1, syltype, :));
    neural = squeeze(allDat_xcov(scalwind, i, syltype, :));
    
    lt_subplot(4,2,i); hold on;
    title(['bin num: ' num2str(i)]);
    xlabel('neural change (rel base');
    ylabel('learn (z, rel base)');
    lt_regress(learn, neural, 1, 0, 1, 1);
    % -- overlay with bird color
    for j=1:length(learn)
        pcol = pcols{allbnum(j)};
        lt_plot(neural(j), learn(j), {'Color', pcol'});
    end
end

% --- average over bins
lt_subplot(4,2,Nbins+1); hold on;
title(['average over bins']);
xlabel('neural change (rel base');
ylabel('learn (z, rel base)');

learn = squeeze(mean(allDat_learn(:, 1, syltype, :),1));
neural = squeeze(mean(allDat_xcov(scalwind, :, syltype, :),2));

lt_regress(learn, neural, 1, 0, 1, 1);
for j=1:length(learn)
    pcol = pcols{allbnum(j)};
    lt_plot(neural(j), learn(j), {'Color', pcol'});
end
% plot(neural, learn, 'ok');
%% ============= [LEARNING RATE VS. XCOV CHANGE]

% ============= 1) compute learning rate
allDat_learnrate = nan(size(allDat_learn,4), 1);

for i=1:size(allDat_learn,4)
    
    % fit line
    learn = squeeze(allDat_learn(:,1,syltype, i));
    learn = [0; learn]; % get time rel baseline
    
    time = squeeze(allDat_TimeMedian(1,:, syltype, i));
    time = [0 time]';
    
    [~,~,~,~,~,SummaryStats] = lt_regress(learn, time, 0);
    
    slope = SummaryStats.slope;
    allDat_learnrate(i) = slope;
end


% =============== 2) xcov change correlate with learning rate?
lt_figure; hold on;

nbins = size(allDat_xcov, 2);
for i=1:nbins
    
    lt_subplot(4,2,i); hold on;
    title(['epoch ' num2str(i)]);
    ylabel('learn rate (z/hr)');
    xlabel('xcov (minus base)');
    
    neural = squeeze(allDat_xcov(scalwind, i, syltype, :));
    learn = allDat_learnrate;
    
    lt_regress(learn, neural, 1, 0, 1, 1);
    
    % -- overlay with bird color
    for j=1:length(learn)
        pcol = pcols{allbnum(j)};
        lt_plot(neural(j), learn(j), {'Color', pcol'});
    end
    
end

% ========== average xcov over bins
lt_subplot(4,2,nbins+1); hold on;
title(['average xcov over epochs']);
ylabel('learn rate (z/hr)');
xlabel('xcov (minus base)');

neural = squeeze(mean(allDat_xcov(scalwind, :, syltype, :), 2));
learn = allDat_learnrate;

lt_regress(learn, neural, 1, 0, 1, 1);
% -- overlay with bird color
for j=1:length(learn)
    pcol = pcols{allbnum(j)};
    lt_plot(neural(j), learn(j), {'Color', pcol'});
end


%% ============= [PLOT]

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% =================== LEARNING TRAJECTORY
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('learning (z rel base)')
xlabel('epoch (bin)');
title('one dat per neural')
Y = squeeze(allDat_learn(:,:, syltype, :));

Y = [zeros(1,size(Y,2)); Y];
x = 1:size(Y,1);
plot(x, Y', '-ok');
ymean = mean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(x, ymean, ysem, {'Color', 'r'},1);
xlim([0 size(Y,1)+1]);
lt_plot_zeroline;
% --- p value for each bin
YLIM = ylim;
tmp = [];
for i=2:size(Y,1)
    y = Y(i,:);
    y = unique(y);
    if isempty(tmp)
        tmp = length(y);
    else
        assert(length(y)==tmp);
    end
    p = signrank(y);
    if p<0.1
        lt_plot_text(i, YLIM(2), ['p=' num2str(p)], 'm', 8);
    end
end


% ==================== PLOT LEARNINGM, BUT ONE DATAPOINT PER TRAJECTORY
Yuni = [];
for i=2:size(Y,1)
    y = unique(Y(i,:), 'stable');
    Yuni = [Yuni; y];
end
Yuni = [zeros(1,size(Yuni,2)); Yuni];
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('learning (z rel base)')
xlabel('epoch (bin)');
title('one dat per learn traj');
Y = Yuni;
x = 1:size(Y,1);
plot(x, Y', '-ok');
ymean = mean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(x, ymean, ysem, {'Color', 'r'},1);
xlim([0 size(Y,1)+1]);
lt_plot_zeroline;
% --- p value for each bin
YLIM = ylim;
tmp = [];
for i=2:size(Y,1)
    y = Y(i,:);
    y = unique(y);
    if isempty(tmp)
        tmp = length(y);
    else
        assert(length(y)==tmp);
    end
    p = signrank(y);
    if p<0.1
        lt_plot_text(i, YLIM(2), ['p=' num2str(p)], 'm', 8);
    end
end
p = signrank(Y(end,:), Y(2,:));
% [~, p] = ttest(Y(end,:), Y(2,:));
lt_plot_pvalue(p, 'last vs 1st train bin');
p = signrank(mean(Y([end-1 end],:),1), mean(Y([2 3],:),1));
lt_plot_pvalue(p, 'last2 vs 1st2 train bin', 2);


% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
ylabel('xcov (minus base, z))')
xlabel('epoch (bin)');
Y = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = [zeros(1,size(Y,2)); Y];
x = 1:size(Y,1);
plot(x, Y', '-ok');
ymean = mean(Y,2);
ysem = lt_sem(Y');
shadedErrorBar(x, ymean, ysem, {'Color', 'r'},1);
xlim([0 size(Y,1)+1]);
lt_plot_zeroline;
% --- p value for each bin
YLIM = ylim;
tmp = [];
for i=2:size(Y,1)
    y = Y(i,:);
    y = unique(y);
    if isempty(tmp)
        tmp = length(y);
    else
        assert(length(y)==tmp);
    end
    p = signrank(y);
    if p<0.1
        lt_plot_text(i, 1.1*max(y), ['p=' num2str(p)], 'm', 8);
    end
end
p = signrank(Y(end,:), Y(2,:));
lt_plot_pvalue(p, 'last vs 1st train bin');
p = signrank(mean(Y([end-1 end],:),1), mean(Y([2 3],:),1));
lt_plot_pvalue(p, 'last2 vs 1st2 train bin', 2);
%% relative timing of WN and neural change

learn = squeeze(allDat_learn(:,:, syltype, :));
neural = squeeze(allDat_xcov(scalwind, :, syltype, :));

% --- nromalize to endpoint (mean endpoint over all cases)
learn = learn./(mean(learn(end,:)));
neural = neural./(mean(neural(end,:)));

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('xcov (minus base)')
ylabel('learn (z rel base)');
title('normalized to mean endpoint');

pcols = lt_make_plot_colors(size(learn,1), 1, [1 0 0]);
for i=1:size(learn,1)
    y = learn(i,:);
    x = neural(i,:);
    
    plot(x,y, 'x', 'Color', pcols{i});
end
for i=1:size(learn,1)
    y = learn(i,:);
    x = neural(i,:);
    lt_plot(mean(x), mean(y), {'Errors', lt_sem(y), 'Xerrors', lt_sem(x), 'Color', pcols{i}});
end
lt_plot_zeroline;
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('binnum')
ylabel('neural - learn (first norm the mean last bin)');

% ----- neural minus learn
y = neural-learn;
x = 1:size(y,1);
plot(x, y, '-ok');
shadedErrorBar(x, mean(y,2), lt_sem(y'), {'Color', 'r'},1);
xlim([0 size(y,1)+1]);
lt_plot_zeroline;
for i=1:size(y,1)
    p = signrank(y(i,:));
    if p<0.2
        lt_plot_text(i, 1.1*max(y(i,:)), ['p=' num2str(p)], 'm', 9);
    end
end
p = signrank(y(1,:), y(end,:));
% p = signrank(mean(y([1 2],:),1), mean(y([end-1 end],:),1));
lt_plot_pvalue(p, 'last bin vs first', 1);

tmp = 1:size(y,1);
xbin = repmat(tmp', 1, size(y,2));
lt_regress(y(:), xbin(:), 1, 0, 1, 1);
% lt_regress(y, x, 1, 0, 1, 1)
%% WITHIN expt correlation between bin num and neural change (no 0) - good, I just did not plot.
if (0)
    Y = squeeze(allDat_xcov(scalwind, :, syltype, :));
    rho_all = [];
    for i=1:size(Y,2)
        y = Y(:,i);
        x = 1:length(y);
        rho = corr(x',y);
        rho_all = [rho_all; rho];
    end
    p = signrank(rho_all);
end
%% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('xcov (minus base)')
ylabel('learn (z rel base)');

X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));

pcols = lt_make_plot_colors(size(X,2), 0,0);
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    
    plot(x, y, '-', 'Color', pcols{i});
end



%%

Yall = {};
% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('rho (xcov vs. learn, over epochs)')
title('no bin shift');
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
rho_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    
    x = [0; x];
    y = [0; y];
    
    %     x = x(2:end);
    %     y = y(1:end-1);
    %     x = x(1:end-1);
    %     y = y(2:end);
    rho = corr(x,y);
    rho_all = [rho_all; rho];
end
lt_plot_histogram(rho_all);
[~, p]= ttest(rho_all);
% p= signrank(rho_all);
lt_plot_pvalue(p, 'ttest',1);
Yall = [Yall, rho_all];

% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('rho (xcov vs. learn, over epochs)')
title('neural(bin n) vs. learn(bin n+1) [cumulative values]');
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
rho_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    x = [0; x];
    y = [0; y];
    
    %     x = x(2:end);
    %     y = y(1:end-1);
    x = x(1:end-1);
    y = y(2:end);
    rho = corr(x,y);
    rho_all = [rho_all; rho];
end
lt_plot_histogram(rho_all);
[~, p]= ttest(rho_all);
% p= signrank(rho_all);
lt_plot_pvalue(p, 'ttest',1);
Yall = [Yall, rho_all];


% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('rho (xcov vs. learn, over epochs)')
title('neural(bin n+1) vs. learn(bin n) [cumulative values]');
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
rho_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    x = [0; x];
    y = [0; y];
    
    %     x = x(2:end);
    %     y = y(1:end-1);
    x = x(2:end);
    y = y(1:end-1);
    rho = corr(x,y);
    rho_all = [rho_all; rho];
end
lt_plot_histogram(rho_all);
[~, p]= ttest(rho_all);
% p= signrank(rho_all);
lt_plot_pvalue(p, 'ttest',1);
Yall = [Yall, rho_all];


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('noshift -- neura(n)vslearn(n+1) -- near(n+1)vslearn(n)');
title('within expt rho');
ylabel('rho');
x = 1:length(Yall);
lt_plot_MultDist(Yall, x, 1);

p = signrank(Yall{1}, Yall{2});
lt_plot_text(1.5, 0.8, ['p=' num2str(p)]);
p = signrank(Yall{3}, Yall{2});
lt_plot_text(2.5, 0.8, ['p=' num2str(p)]);

%%
% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('rho (xcov vs. learn, over epochs)')
title('no bin shift [bin to bin changes]');
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
rho_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    x = [0; x];
    y = [0; y];
    
    x = diff(x);
    y = diff(y);
    rho = corr(x,y);
    rho_all = [rho_all; rho];
end
lt_plot_histogram(rho_all);
[~, p]= ttest(rho_all);
% p= signrank(rho_all);
lt_plot_pvalue(p, 'ttest',1);



% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('rho (xcov vs. learn, over epochs)')
title('neural(bin n) vs. learn(bin n+1) [bin to bin changes]');
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
rho_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    x = [0; x];
    y = [0; y];
    
    x = x(1:end-1);
    y = y(2:end);
    
    x = diff(x);
    y = diff(y);
    rho = corr(x,y);
    rho_all = [rho_all; rho];
end
lt_plot_histogram(rho_all);
[~, p]= ttest(rho_all);
% p= signrank(rho_all);
lt_plot_pvalue(p, 'ttest',1);







% ===================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('rho (xcov vs. learn, over epochs)')
title('neural(bin n+1) vs. learn(bin n) [bin to bin changes]');
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
rho_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    x = [0; x];
    y = [0; y];
    
    x = x(2:end);
    y = y(1:end-1);
    
    x = diff(x);
    y = diff(y);
    rho = corr(x,y);
    rho_all = [rho_all; rho];
end
lt_plot_histogram(rho_all);
[~, p]= ttest(rho_all);
% p= signrank(rho_all);
lt_plot_pvalue(p, 'ttest',1);








% ====================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
title('epoch-epoch differences (scaled from 0-1)');
xlabel('xcov');
ylabel('learn');

Xdiff_all = [];
Ydiff_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    
    x = [0; x];
    y = [0; y];
    
    % -- scale by min and max
    x = (x - min(x))./(max(x)-min(x));
    y = (y - min(y))./(max(y)-min(y));
    
    %     x = diff(x);
    %     y = diff(y);
    %     x = x(1:end-1);
    %     y = y(2:end);
    %     x = x(2:end);
    %     y = y(1:end-1);
    
    % -- collect differences
    Xdiff_all = [Xdiff_all; diff(x)];
    Ydiff_all = [Ydiff_all; diff(y)];
end
lt_regress(Ydiff_all, Xdiff_all, 1, 0, 1, 1);


% ====================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
title('epoch-epoch differences (scaled from 0-1)');
xlabel('xcov(bin N)');
ylabel('learn(bin N+1)');

Xdiff_all = [];
Ydiff_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    
    x = [0; x];
    y = [0; y];
    
    % -- scale by min and max
    x = (x - min(x))./(max(x)-min(x));
    y = (y - min(y))./(max(y)-min(y));
    
    %     x = diff(x);
    %     y = diff(y);
    x = x(1:end-1);
    y = y(2:end);
    %     x = x(2:end);
    %     y = y(1:end-1);
    
    % -- collect differences
    Xdiff_all = [Xdiff_all; diff(x)];
    Ydiff_all = [Ydiff_all; diff(y)];
end
lt_regress(Ydiff_all, Xdiff_all, 1, 0, 1, 1);

% ====================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
title('epoch-epoch differences (scaled from 0-1)');
xlabel('xcov(bin N+1)');
ylabel('learn(bin N)');

Xdiff_all = [];
Ydiff_all = [];
for i=1:size(X,2)
    x = X(:,i);
    y = Y(:,i);
    
    x = [0; x];
    y = [0; y];
    
    % -- scale by min and max
    x = (x - min(x))./(max(x)-min(x));
    y = (y - min(y))./(max(y)-min(y));
    
    %     x = diff(x);
    %     y = diff(y);
    %     x = x(1:end-1);
    %     y = y(2:end);
    x = x(2:end);
    y = y(1:end-1);
    
    % -- collect differences
    Xdiff_all = [Xdiff_all; diff(x)];
    Ydiff_all = [Ydiff_all; diff(y)];
end
lt_regress(Ydiff_all, Xdiff_all, 1, 0, 1, 1);


%% ================== [SHUFFLE ANALYSIS] LEARN AND NEURAL correlate?
lt_figure; hold on;

nshuff = 1000;

% ===============
lt_subplot(4,2,1); hold on;
title('no shift, cumul changes');
xlabel('rho (neural vs. learn)')
doshift = 0;
bintobinchanges =0; % then this, instead of cumulative values.
[rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff);

xperm = mean(rhoperm_all,2);
lt_plot_histogram(xperm);
xdat = mean(rhodat);
line([xdat xdat], ylim, 'Color', 'r');
p = (1+sum(mean(rhoperm_all,2)>mean(rhodat)))/(1+size(rhoperm_all,1));
lt_plot_pvalue(p, 'vs shuff', 1);


% ===============
lt_subplot(4,2,2); hold on;
title('with shift (neural N vs. learn N+1), cumul changes');
doshift = 1;
bintobinchanges =0; % then this, instead of cumulative values.
[rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff);

xperm = mean(rhoperm_all,2);
lt_plot_histogram(xperm);
xdat = mean(rhodat);
line([xdat xdat], ylim, 'Color', 'r');
p = (1+sum(mean(rhoperm_all,2)>mean(rhodat)))/(1+size(rhoperm_all,1));
lt_plot_pvalue(p, 'vs shuff', 1);


% ===============
lt_subplot(4,2,3); hold on;
title('no shift, bin to bin changes');
doshift = 0;
bintobinchanges =1; % then this, instead of cumulative values.
[rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff);

xperm = mean(rhoperm_all,2);
lt_plot_histogram(xperm);
xdat = mean(rhodat);
line([xdat xdat], ylim, 'Color', 'r');
p = (1+sum(mean(rhoperm_all,2)>mean(rhodat)))/(1+size(rhoperm_all,1));
lt_plot_pvalue(p, 'vs shuff', 1);


% ==============
lt_subplot(4,2,4); hold on;
title('with shift (neural N vs. learn N+1), bintobin changes');
doshift = 1;
bintobinchanges =1; % then this, instead of cumulative values.
[rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff);

xperm = mean(rhoperm_all,2);
lt_plot_histogram(xperm);
xdat = mean(rhodat);
line([xdat xdat], ylim, 'Color', 'r');
p = (1+sum(mean(rhoperm_all,2)>mean(rhodat)))/(1+size(rhoperm_all,1));
lt_plot_pvalue(p, 'vs shuff', 1);




% ===============
lt_subplot(4,2,5); hold on;
title('with shift (neural N+1 vs. learn N), cumul changes');
doshift = 2;
bintobinchanges =0; % then this, instead of cumulative values.
[rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff);

xperm = mean(rhoperm_all,2);
lt_plot_histogram(xperm);
xdat = mean(rhodat);
line([xdat xdat], ylim, 'Color', 'r');
p = (1+sum(mean(rhoperm_all,2)>mean(rhodat)))/(1+size(rhoperm_all,1));
lt_plot_pvalue(p, 'vs shuff', 1);


% ==============
lt_subplot(4,2,6); hold on;
title('with shift (neural N+1 vs. learn N), bintobin changes');
doshift = 2;
bintobinchanges =1; % then this, instead of cumulative values.
[rhodat, rhoperm_all] = lt_neural_POPLEARN_Xcov_EpochScal_sub(allDat_xcov, allDat_learn, ...
    allbnum, allenum, allswnum, doshift, bintobinchanges, scalwind, syltype, ...
    nshuff);

xperm = mean(rhoperm_all,2);
lt_plot_histogram(xperm);
xdat = mean(rhodat);
line([xdat xdat], ylim, 'Color', 'r');
p = (1+sum(mean(rhoperm_all,2)>mean(rhodat)))/(1+size(rhoperm_all,1));
lt_plot_pvalue(p, 'vs shuff', 1);




%% =============== MIXED EFFECTS MODEL - CORRELATION BETWEEN LEARNING AND NEURAL?
% ========================= cumulative learning/change, from baseline.

doshift = 0; % then shift by one bin (xcov predict leanring?)
takeincrements = 0; % then is change from bin to bin.
scalewithinexpt = 0;
add0 = 0;

X = squeeze(allDat_xcov(scalwind, :, syltype, :));
Y = squeeze(allDat_learn(:,:, syltype, :));
times = squeeze(allDat_TimeMedian(1,:, syltype, :)); % --- time of bin


% --- add 0, since they are all deviation from baseline
if add0 == 1
    Xneural = [zeros(1, size(X,2)); X];
    Ylearn = [zeros(1, size(Y,2)); Y];
    times = [zeros(1,size(times,2)); times];
else
    Xneural = X;
    Ylearn = Y;
    times = times;
    
end
if scalewithinexpt==1
    Ylearn = (Ylearn - min(Ylearn, [], 1))./(max(Ylearn, [], 1)-min(Ylearn, [], 1));
    Xneural = (Xneural - min(Xneural, [], 1))./(max(Xneural, [], 1)-min(Xneural, [], 1));
    % Xneural = [zeros(1, size(X,2)); X];
    % Ylearn = [zeros(1, size(Y,2)); Y];
end

if doshift==1
    Xneural = Xneural(1:end-1,:);
    Ylearn = Ylearn(2:end,:);
    times = times(1:end-1,:);
elseif doshift==2
    Xneural = Xneural(2:end,:);
    Ylearn = Ylearn(1:end-1,:);
    times = times(1:end-1,:);
end


if takeincrements==1
    % ---- take differences (i.e. increments from bin to bin)
    Xneural = diff(Xneural, 1, 1);
    Ylearn = diff(Ylearn,1,1);
    times = diff(times, 1,1);
end

tmp = 1:size(Xneural,1);
binnum = repmat(tmp', 1, size(Xneural,2));

% --- one ind for each expt
exptID = lt_tools_grp2idx({allbnum, allenum, allswnum});
exptID = repmat(exptID', size(Xneural,1), 1);


% -- birdnum
bnum = repmat(allbnum', size(Xneural,1), 1);

% ---- expand
Xneural = Xneural(:);
Ylearn = Ylearn(:);
binnum = binnum(:);
exptID = categorical(exptID(:));
bnum = categorical(bnum(:));
times = times(:);


dat = table(Xneural, Ylearn, binnum, exptID, bnum, times);
formula = 'Ylearn ~ Xneural';
% formula = 'Ylearn ~ Xneural + (Xneural|exptID)';
% formula = 'Ylearn ~ Xneural + (Xneural|bnum)';
% formula = 'Ylearn ~ binnum + Xneural';
lme = fitlme(dat, formula)


% ========================= predict xcov vhange, learning or time?
% formula = 'Xneural ~ times';
% formula = 'Xneural ~ times + Ylearn';
formula = 'Xneural ~ times + Ylearn + (times|exptID) + (Ylearn|exptID)';
% formula = 'Xneural ~ times + (times|exptID)';
lme = fitlme(dat, formula)


% ======= plot results from lme
xnames = lme.Coefficients.Name;
beta = lme.Coefficients.Estimate;
pval = lme.Coefficients.pValue;
beta_SE = lme.Coefficients.SE;

lt_figure; hold on;
x = 1:length(beta);
lt_plot_bar(x, beta, {'Errors', beta_SE})


