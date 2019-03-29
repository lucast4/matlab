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



    
    