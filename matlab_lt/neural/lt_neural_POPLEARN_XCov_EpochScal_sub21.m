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

% =============== DO LME
indsexpt = categorical(lt_tools_grp2idx({allbnum, allenum}));
xcov_ad_minus_nonad_without_subtr_base = y(2,:)';
dat = table(xcov_ad_minus_nonad_without_subtr_base, indsexpt);
formula = 'xcov_ad_minus_nonad_without_subtr_base ~ 1 + (1|indsexpt)';
lme = fitlme(dat, formula)

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


%%
plotraw = 1;
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
    if all(isnan(ymean))
        continue
    end
    if size(y,1)<2
        continue
    end
    lt_plot(x, ymean, {'Errors', ysem, 'LineStyle', '-'});
    
    if all(isnan(y(:,1)))
        continue
    end
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

%% ========
figcount=1;
subplotrows=1;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

%% ================== PLOT HISTOGRAM OF ALL SYL/CHAN COMBINATIONS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('TARG -- DIFF)');
ylabel('xcov (ad - non, mnus base)');
title('hist is chan x syl,');
lt_plot_annotation(1, 'if empty then conditions not correct', 'r');
binsize = 0.05;

if useAd_Nonad_Average_forBaseline==0 & keptallinds==1
    %  otherwise the extracted data is probably different from original OUTSTRUCT.
    assert(size(OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i},2)==1, 'this code only works if only have one scalar window defined... [an issue with squeeze..]');
    
    % === 1) extract all change from baseline
    tmp = nan(size(OUTSTRUCT_XCOV.bnum));
    for i=1:length(OUTSTRUCT_XCOV.bnum)
        
        ywn_AdMinNon = nanmean(OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit{i}{scalwind}(wnbins-1,2),1) ...
            - nanmean(OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit{i}{scalwind}(wnbins-1,1),1);
        
        ybase_AdMinNon = OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i}(2) ...
            - OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i}(1);
        
        
        
        y = ywn_AdMinNon - ybase_AdMinNon;
        y = y*OUTSTRUCT_XCOV.learndirTarg(i); % so positive is in adaptive direction
        
        tmp(i) = y;
    end
    
    Y = {};
    Y{1} = tmp(OUTSTRUCT_XCOV.istarg==1);
    Y{2} = tmp(OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==0);
    
    % ============ WHAT IS GOOD BINSIZE
    y = cell2mat(cellfun(@(x)x', Y, 'UniformOutput', 0));
    tmp = round(100*max(abs(y)))/100;
    tmp = tmp+(binsize-mod(tmp, binsize));
    % tmp = tmp+binsize/2;
    xcenter = -tmp:binsize:tmp;
    
    distributionPlot(Y, 'xValues', [1 3], 'showMM', 4, 'addSpread', 0, 'color', ...
        'b', 'histOpt', 0, 'divFactor', xcenter);
    
    
    
    
end

%% ================== PLOT HISTOGRAM OF ALL SYL/CHAN COMBINATIONS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('nonadapt -- adapt)');
ylabel('xcov (ad - non, mnus base)');
title('hist is chan x syl, line dat is each neuron 1 point');
lt_plot_annotation(1, 'if empty then conditions not correct', 'r');
binsize = 0.05;



if useAd_Nonad_Average_forBaseline==0 & keptallinds==1
    %  otherwise the extracted data is probably different from original OUTSTRUCT.
    assert(size(OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i},2)==1, 'this code only works if only have one scalar window defined... [an issue with squeeze..]');
    
    % === 1) extract all change from baseline
    Yall_ad = nan(size(OUTSTRUCT_XCOV.bnum));
    Yall_nonad = nan(size(OUTSTRUCT_XCOV.bnum));
    for i=1:length(OUTSTRUCT_XCOV.bnum)
        
        %         ywn_AdMinNon = nanmean(OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit{i}{scalwind}(wnbins-1,2),1) ...
        %             - nanmean(OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit{i}{scalwind}(wnbins-1,1),1);
        %
        %         ybase_AdMinNon = OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i}(2) ...
        %             - OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i}(1);
        
        yAd = nanmean(OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit{i}{scalwind}(wnbins-1, 2),1) ...
            - OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i}(2); % change from baseline
        yNon = nanmean(OUTSTRUCT_XCOV.xcovscalEpochs_window_FFsplit{i}{scalwind}(wnbins-1, 1),1) ...
            - OUTSTRUCT_XCOV.xcovscalBase_window_FFsplit{i}(1);
        
        if OUTSTRUCT_XCOV.learndirTarg(i)==-1
            % -- flip ad and nonad
            A = yAd;
            yAd = yNon;
            yNon = A;
        end
        Yall_ad(i) = yAd;
        Yall_nonad(i) = yNon;
        
    end
    
    Y = {};
    
    indsthis = OUTSTRUCT_XCOV.istarg==1;
    Y{1} = Yall_nonad(indsthis);
    Y{2} = Yall_ad(indsthis);
    lt_plot_annotation(2, ['N=' num2str(cellfun(@(x)length(x), Y))], 'm');
    
    % ============ WHAT IS GOOD BINSIZE
    y = cell2mat(cellfun(@(x)x', Y, 'UniformOutput', 0));
    tmp = round(100*max(abs(y)))/100;
    tmp = tmp+(binsize-mod(tmp, binsize));
    % tmp = tmp+binsize/2;
    xcenter = -tmp:binsize:tmp;
    
    distributionPlot(Y, 'xValues', [1 3], 'showMM', 4, 'addSpread', 0, 'color', ...
        'b', 'histOpt', 0, 'divFactor', xcenter);
    
    
    
    % ######################## OVERLAY EACH BIRD
        % -- wn
        Y = [Ynonadapt(:,2) Yadapt(:,2)] - [Ynonadapt(:,1) Yadapt(:,1)];
        x = [1 3];
        % if plotraw==1
        % plot(x, Y', '-', 'Color', pcol);
        % end
        ymean = nanmean(Y,1);
        ysem = lt_sem(Y);
        lt_plot_bar(x, ymean, {'Errors', ysem});
    
        % p = signrank(Y(:,1), Y(:,2));
        % lt_plot_pvalue(p, 'vs', 1);
        % for i=1:size(Y,2)
        %     p = signrank(Y(:,i));
        %     if p<0.15
        %         lt_plot_text(i, max(Y(:,i)), ['p=' num2str(p)]);
        %     end
        % end
    
        % -- overlay each bird
        for i=1:max(Ybnum)
            y = Y(Ybnum==i,:);
            if isempty(y)
                continue
            end
    
            ymean = nanmean(y,1);
            ysem = lt_sem(y);
            if all(isnan(ymean))
                continue
            end
            if size(y,1)<2
                continue
            end
            lt_plot(x, ymean, {'Errors', ysem, 'LineStyle', '-'});
    
            if all(isnan(y(:,1)))
                continue
            end
            p1 = signrank(y(:,1), y(:,2));
            %     lt_plot_pvalue(p, 'vs', 1);
            for ii=2
                p = signrank(y(:,ii));
                lt_plot_text(x(end), ymean(end), ['p(vs)= ' num2str(p1), 'p(2)=' num2str(p) '[bnum' num2str(i)]);
            end
        end
        % --
        xlim([-1 5]);
        ylim([-0.7 0.7]);
    
end
%% ======================= lme modeling, account for shared experiment

% =========================================== COLLECT DATA
y1 = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, 1, 1, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, 1, 1, indsthis),1));
Ynonadapt = y2 - y1;

y1 = squeeze(allDat_BaseEpochs_FFsplits_adaptivedir(1, 2, 1, indsthis));
y2 = squeeze(nanmean(allDat_BaseEpochs_FFsplits_adaptivedir(wnbins, 2, 1, indsthis),1));
Yadapt = y2 - y1;

Yresponse = Yadapt - Ynonadapt; % difference (adapt-nonadapt) of differences (Wn - base)
Ybnum = categorical(allbnum(indsthis));
Yenum = categorical(allenum(indsthis));
YenumUnique = categorical(lt_tools_grp2idx({allbnum(indsthis), allenum(indsthis)}));

lt_figure; hold on;
dat = table(Yresponse, Ynonadapt, Yadapt, Ybnum, YenumUnique);

% ================ adapt - nonadapt
% formula = 'Yresponse ~ 1 + (1|Ybnum) + (1|Ybnum:YenumUnique)';
formula = 'Yresponse ~ 1 + (1|YenumUnique)';
% formula = 'Yresponse ~ 1 + (1|YenumUnique:Ybnum)';
% % formula = 'Yresponse ~ 1 + (1|YenumUnique) + (1|Ybnum)';
% formula = 'Yresponse ~ 1 + (1|Ybnum)';
lme = fitlme(dat, formula, 'StartMethod', 'random')

lt_subplot(3,3,1); hold on;
xlabel('expt id');
ylabel('change in xcov, adapt-nonadap');
plot(YenumUnique, Yresponse, 'o');
lt_plot_zeroline;

lt_subplot(3,3,2); hold on;
title(['lme, ' formula]);
ylabel('effect');
lt_tools_lme_plotEffects(lme, 0);

% -------- plot model predictions
lt_subplot(3,3,3); hold on;
title('model fits');
plot(dat.YenumUnique, lme.fitted, 'ok');


% ================ adapt
formula = 'Yadapt ~ 1 + (1|YenumUnique)';
% formula = 'Yadapt ~ 1 + (1|YenumUnique:Ybnum)';
% % formula = 'Yresponse ~ 1 + (1|YenumUnique) + (1|Ybnum)';
% formula = 'Yresponse ~ 1 + (1|Ybnum)';
lme = fitlme(dat, formula);

lt_subplot(3,3,4); hold on;
xlabel('expt id');
ylabel('change in xcov, adapt');
plot(YenumUnique, Yadapt, 'o');
lt_plot_zeroline;

lt_subplot(3,3,5); hold on;
title(['lme, ' formula]);
ylabel('effect');
lt_tools_lme_plotEffects(lme, 0);

% -------- plot model predictions
lt_subplot(3,3,6); hold on;
title('model fits');
plot(dat.YenumUnique, lme.fitted, 'ok');

% ================ nonadapt
formula = 'Ynonadapt ~ 1 + (1|YenumUnique)';
% % formula = 'Yresponse ~ 1 + (1|YenumUnique) + (1|Ybnum)';
% formula = 'Yresponse ~ 1 + (1|Ybnum)';
lme = fitlme(dat, formula);

lt_subplot(3,3,7); hold on;
xlabel('expt id');
ylabel('change in xcov, nonadapt');
plot(YenumUnique, Ynonadapt, 'o');
lt_plot_zeroline;

lt_subplot(3,3,8); hold on;
title(['lme, ' formula]);
ylabel('effect');
lt_tools_lme_plotEffects(lme, 0);

% -------- plot model predictions
lt_subplot(3,3,9); hold on;
title('model fits');
plot(dat.YenumUnique, lme.fitted, 'ok');
