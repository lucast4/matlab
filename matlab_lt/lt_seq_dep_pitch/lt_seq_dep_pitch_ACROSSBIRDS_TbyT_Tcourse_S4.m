function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S4(NtargCell, NnontargCell, ...
    FFbinnedCell, TbinnedCell, FFsinglebinMat, indstoplot)


%% =================== [PLOT] SHOW SAMPLE SIZE

lt_figure; hold on;

% ############################ NUMBER OF RENDS IN BIN FOR TARG AND NONTARG
lt_subplot(4,2,1); hold on;
xlabel('targ/hi -- targ/lo -- nontarg/hi -- nontarg/lo');
ylabel('log2(N)');
title('number renditions');

Nall = [];

% ---------------- targ, high
N = NtargCell(indstoplot,2);
N = cellfun(@mean, N);
Nall = [Nall N];

% ---------------- targ, lo
N = NtargCell(indstoplot,1);
N = cellfun(@mean, N);
Nall = [Nall N];

% ---------------- Nontarg, high
N = NnontargCell(indstoplot,2);
N = cellfun(@mean, N);
Nall = [Nall N];

% ---------------- Nontarg, low
N = NnontargCell(indstoplot,1);
N = cellfun(@mean, N);
Nall = [Nall N];

boxplot(log2(Nall));

% ####################### SCATTER PLOT
lt_subplot(4,2,3); hold on;
title('targ');
xlabel('lo dens (log2N)');
ylabel('hi dens');

% ---------- TARG
N = NtargCell(indstoplot,:);
N = log2(cellfun(@mean, N));
plot(N(:,1), N(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'k', -1);


lt_subplot(4,2,4); hold on;
title('Nontarg');
xlabel('lo dens');
ylabel('hi dens');

% ---------- TARG
N = NnontargCell(indstoplot,:);
N = log2(cellfun(@mean, N));
plot(N(:,1), N(:,2), 'ok');
lt_plot_makesquare_plot45line(gca, 'k', -1);


%% ================= [PLOT] FFdev
lt_figure; hold on;


% ############################ BINNED FF DEV 
lt_subplot(3,1,1); hold on;

for i=indstoplot
   
    % --- low dens
    indcol = 1;
    plotcol = 'b';
    
    x = TbinnedCell{i,indcol};
    y = FFbinnedCell{i, indcol};
    plot(x+0.2*indcol,y, '-x', 'Color', plotcol);
    
    % --- hi dens
    indcol = 2;
    plotcol = 'r';
    
    x = TbinnedCell{i,indcol};
    y = FFbinnedCell{i, indcol};
    plot(x+0.2*indcol,y, '-x', 'Color', plotcol);
    
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);

% ############################ BINNED FF DEV [mean across syls] 
lt_subplot(3,1,2); hold on;

% ============================ low dens
indcol = 1;
plotcol = 'b';

% run
X = [];
Y = [];
for i=indstoplot
       
    x = TbinnedCell{i,indcol};
    y = FFbinnedCell{i, indcol};
    
    X = [X x];
    Y = [Y y'];
end

[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
xmean = unique(X);
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', plotcol});

% ============================ hi dens
indcol = 2;
plotcol = 'r';

% run
X = [];
Y = [];
for i=indstoplot
       
    x = TbinnedCell{i,indcol};
    y = FFbinnedCell{i, indcol};
    
    X = [X x];
    Y = [Y y'];
end

[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
xmean = unique(X);
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', plotcol});

% =========== formatn
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);


% ###################################### TAKE ALL RENDS IN ONE BIN
lt_subplot(3,2,6); hold on;

Y = FFsinglebinMat(indstoplot,:);
X = [1 2];
plot(X, Y', '-', 'Color', [0.7 0.7 0.7]);
% -- means
lt_plot(X+0.1, mean(Y), {'Errors', lt_sem(Y)});

xlim([0 3]);
lt_plot_zeroline;
% --- signifnicace
[~, p] = ttest(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'ttest', 1);
