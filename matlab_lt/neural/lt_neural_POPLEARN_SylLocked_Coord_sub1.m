figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

[~, indsort] = sort(rhoall);
nplot = 3;

% --- collect all rho and power to plot in end
rhotmp = [];
powtmp = [];
frtmp = [];

% ====== LOW CORREALTION TRIALS
trialstoplot = indsort(1:nplot)';
ptit = 'LOW CORR';
if strcmp(bregionthis, 'LMAN')
    pcol = [0.2 0.6 0.2];
elseif strcmp(bregionthis, 'RA')
    pcol = 'r';
end

for tt=trialstoplot
    
    spk = cellfun(@(x)(x(tt)), spkdat);
    fr = cellfun(@(x)(x(:, tt)),  frmat, 'UniformOutput', 0);
    lfp = cellfun(@(x)x(:, tt), lfpdat,  'UniformOutput', 0);
    rho = rhoall(tt);
    frmn = frmean(tt);
    pow = lfppow_all(tt);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([ptit ', tr#' num2str(tt)]);
    ylabel('0.1*lfp; spikes');
    % -- plot spikes
    for j=1:length(spk)
        lt_neural_PLOT_rasterline(spk{j}, j, pcol);
    end
    % -- plot lfp
    xlfp = PARAMS.THIS.lfpx;
    for j=1:length(lfp)
        y = 0.1*lfp{j};
        plot(xlfp, y, 'Color', pcol);
    end
    % -- annotate rho and power
    tstr = ['rho=' num2str(rho) ',pow=' num2str(pow)];
    lt_plot_annotation(1, tstr, 'b');
    rhotmp = [rhotmp; rho];
    powtmp = [powtmp ; pow];
    frtmp = [frtmp; frmn];
    lt_plot_zeroline_vert;
    lt_plot_zeroline;
end

% ====== MEDIAN CORREALTION TRIALS
trialstoplot = indsort(floor(length(indsort)/2):floor(length(indsort)/2)+nplot-4)';
ptit = 'MID CORR';

for tt=trialstoplot
    
    spk = cellfun(@(x)(x(tt)), spkdat);
    fr = cellfun(@(x)(x(:, tt)),  frmat, 'UniformOutput', 0);
    lfp = cellfun(@(x)x(:, tt), lfpdat,  'UniformOutput', 0);
    frmn = frmean(tt);
    rho = rhoall(tt);
    pow = lfppow_all(tt);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([ptit ', tr#' num2str(tt)]);
    ylabel('0.1*lfp; spikes');
    % -- plot spikes
    for j=1:length(spk)
        lt_neural_PLOT_rasterline(spk{j}, j, pcol);
    end
    % -- plot lfp
    xlfp = PARAMS.THIS.lfpx;
    for j=1:length(lfp)
        y = 0.1*lfp{j};
        plot(xlfp, y, 'Color', pcol);
    end
    % -- annotate rho and power
    tstr = ['rho=' num2str(rho) ',pow=' num2str(pow)];
    lt_plot_annotation(1, tstr, 'b');
    rhotmp = [rhotmp; rho];
    powtmp = [powtmp ; pow];
    frtmp = [frtmp; frmn];
    lt_plot_zeroline_vert;
    lt_plot_zeroline;
end


% ====== HIGH CORREALTION TRIALS
trialstoplot = indsort(end-nplot+1:end)';
ptit = 'HIGH CORR';

for tt=trialstoplot
    
    spk = cellfun(@(x)(x(tt)), spkdat);
    fr = cellfun(@(x)(x(:, tt)),  frmat, 'UniformOutput', 0);
    lfp = cellfun(@(x)x(:, tt), lfpdat,  'UniformOutput', 0);
    rho = rhoall(tt);
    frmn = frmean(tt);
    pow = lfppow_all(tt);
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title([ptit ', tr#' num2str(tt)]);
    ylabel('0.1*lfp; spikes');
    % -- plot spikes
    for j=1:length(spk)
        lt_neural_PLOT_rasterline(spk{j}, j, pcol);
    end
    % -- plot lfp
    xlfp = PARAMS.THIS.lfpx;
    for j=1:length(lfp)
        y = 0.1*lfp{j};
        plot(xlfp, y, 'Color', pcol);
    end
    % -- annotate rho and power
    tstr = ['rho=' num2str(rho) ',pow=' num2str(pow)];
    lt_plot_annotation(1, tstr, 'b');
    rhotmp = [rhotmp; rho];
    powtmp = [powtmp ; pow];
    frtmp = [frtmp; frmn];
    lt_plot_zeroline_vert;
    lt_plot_zeroline;
end


% =================== PLOT BAR SHOWING ALL EXTRACTED PARAMS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('selected trials in order');
x = 1:length(rhotmp);
plot(x-0.2, rhotmp, 'or');
plot(x+0.2, 0.1*powtmp, 'ob');
ylabel('r=corr; b=pow(0.1*)');
lt_plot_zeroline;

% =================== correlation all trials
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all trials');
xlabel('corr');
ylabel('pow(1x)');
lt_regress(lfppow_all, rhoall, 1, 0, 1, 1, pcol);

% =================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all trials');
xlabel('mean fr');
ylabel('pow(1x)');
lt_regress(lfppow_all, frmean, 1, 0, 1, 1, pcol);

% ========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all trials');
xlabel('mean fr');
ylabel('corr');
lt_regress(rhoall, frmean, 1, 0, 1, 1, pcol);

% ================ FORMAT ALL
linkaxes(hsplots, 'xy');