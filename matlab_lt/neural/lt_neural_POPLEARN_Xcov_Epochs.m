function lt_neural_POPLEARN_Xcov_Epochs(OUTSTRUCT, OUTSTRUCT_XCOV, PARAMS, ...
    onlygoodexpt, SwitchStruct, dattype)
%% lt 3/5/19 - divides up training into epochs -- here PLOTS
XLIM = [-0.05 0.05];

%% ====== [PREPROCESS] only plot good experiments

if onlygoodexpt==1
    
    % ===== filter outstruct
    [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
    
    % ===== filter outstruct
    [OUTSTRUCT indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, ...
        SwitchStruct, 'xcov_spikes');
    
end


%% ====== [PREPROCESS] average, for each experiment

%% group based on syl type

if strcmp(dattype, 'switch')
    %     fieldtoget = 'XcovgramBase';
    %     [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_base] = ...
    %         lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    %
    %     fieldtoget = 'XcovgramWN';
    %     [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_wn] = ...
    %         lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
elseif strcmp(dattype, 'chan')
    fieldtoget = 'Xcovslice_epochs';
    [allbnum, allenum, allswnum, allDat_epochs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'XcovBase';
    [allbnum, allenum, allswnum, allDat_base] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'XcovWN';
    [allbnum, allenum, allswnum, allDat_wn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
end


%% ===== PLOT EACH EXPERIMENT, each syl

pcols = lt_make_plot_colors(size(allDat_epochs,1), 1, [1 0 0]);
t = PARAMS.Xcov_ccLags;

figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% =======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG');
xlabel('lag (LMAN-RA)');
ind3 = 1; % target

% ---- 1) plot baseline
dat = squeeze(allDat_base(:,:, ind3, :));
ymean = mean(dat, 2);
plot(t, ymean, ':', 'Color', 'k', 'LineWidth', 2);

% ----- 2) WN epochs
dat = squeeze(allDat_epochs(:, :, ind3, :));
for i=1:size(dat,1)
    y = squeeze(dat(i,:,:));
    ymean = mean(y,2);
    plot(t, ymean, 'Color', pcols{i}, 'LineWidth', 2);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
end

% % --- 3) plot WN
% dat = squeeze(allDat_wn(:,:, ind3, :));
% ymean = mean(dat, 2);
% plot(t, ymean, ':', 'Color', 'm', 'LineWidth', 2);

% -- fomrat
axis tight;
xlim(XLIM);
lt_plot_zeroline_vert;


% =======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('TARG');
xlabel('lag (LMAN-RA)');

% ---- 1) get baseline
datBase = squeeze(allDat_base(:,:, ind3, :));

% ----- 2) WN epochs
dat = squeeze(allDat_epochs(:, :, ind3, :));
for i=1:size(dat,1)
    y = squeeze(dat(i,:,:));
    y = y - datBase;
    ymean = mean(y,2);
    plot(t, ymean, 'Color', pcols{i}, 'LineWidth', 2);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
end

% % --- 3) plot WN
% dat = squeeze(allDat_wn(:,:, ind3, :));
% ymean = mean(dat, 2);
% plot(t, ymean, ':', 'Color', 'm', 'LineWidth', 2);

% -- fomrat
axis tight;
xlim(XLIM);
lt_plot_zeroline_vert;
lt_plot_zeroline;



% =================== PLOT EACH EPOCH SEPARATELY
% ------- RUN
dat = squeeze(allDat_epochs(:, :, ind3, :));
for i=1:size(dat,1)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('TARG');
    xlabel('lag (LMAN-RA)');
    ylabel(['epoch(WN): ' num2str(i)]);
    
    
    % ------------ baseline
    datbase = squeeze(allDat_base(:,:, ind3, :));
    ymean = mean(datbase, 2);
    ysem = lt_sem(datbase');
    shadedErrorBar(t, ymean, ysem, {'Color', 'k', 'LineStyle', '--'}, 1);
    
    
    % -------- epoch (this)
    y = squeeze(dat(i,:,:));
    ymean = mean(y,2);
    ysem = lt_sem(y');
    %     shadedErrorBar(t, ymean, ysem, {'Color', pcols{i}}, 1);
    shadedErrorBar(t, ymean, ysem, {'Color', 'r'}, 1);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
    
    
    % -------- FINAL
    datWN = squeeze(allDat_wn(:,:, ind3, :));
    ymean = mean(datWN, 2);
    ysem = lt_sem(datWN');
    plot(t, ymean, 'Color', 'm');
    
    % -- fomrat
    axis tight;
    xlim(XLIM);
    lt_plot_zeroline_vert;
    
end







% =======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF');
xlabel('lag (LMAN-RA)');
ind3 = 3; % target

% ---- 1) plot baseline
dat = squeeze(allDat_base(:,:, ind3, :));
ymean = mean(dat, 2);
plot(t, ymean, ':', 'Color', 'k', 'LineWidth', 2);

% ----- 2) WN epochs
dat = squeeze(allDat_epochs(:, :, ind3, :));
for i=1:size(dat,1)
    y = squeeze(dat(i,:,:));
    ymean = mean(y,2);
    plot(t, ymean, 'Color', pcols{i}, 'LineWidth', 2);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
end

% % --- 3) plot WN
% dat = squeeze(allDat_wn(:,:, ind3, :));
% ymean = mean(dat, 2);
% plot(t, ymean, ':', 'Color', 'm', 'LineWidth', 2);

% -- fomrat
axis tight;
xlim(XLIM);
lt_plot_zeroline_vert;


% =======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('TARG');
xlabel('lag (LMAN-RA)');

% ---- 1) get baseline
datBase = squeeze(allDat_base(:,:, ind3, :));

% ----- 2) WN epochs
dat = squeeze(allDat_epochs(:, :, ind3, :));
for i=1:size(dat,1)
    y = squeeze(dat(i,:,:));
    y = y - datBase;
    ymean = mean(y,2);
    plot(t, ymean, 'Color', pcols{i}, 'LineWidth', 2);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
end

% % --- 3) plot WN
% dat = squeeze(allDat_wn(:,:, ind3, :));
% ymean = mean(dat, 2);
% plot(t, ymean, ':', 'Color', 'm', 'LineWidth', 2);

% -- fomrat
axis tight;
xlim(XLIM);
lt_plot_zeroline_vert;
lt_plot_zeroline;



% =================== PLOT EACH EPOCH SEPARATELY
% ------- RUN
dat = squeeze(allDat_epochs(:, :, ind3, :));
for i=1:size(dat,1)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     title('TARG');
    xlabel('lag (LMAN-RA)');
    ylabel(['epoch(WN): ' num2str(i)]);
    
    
    % ------------ baseline
    datbase = squeeze(allDat_base(:,:, ind3, :));
    ymean = mean(datbase, 2);
    ysem = lt_sem(datbase');
    shadedErrorBar(t, ymean, ysem, {'Color', 'k', 'LineStyle', '--'}, 1);
    
    
    % -------- epoch (this)
    y = squeeze(dat(i,:,:));
    ymean = mean(y,2);
    ysem = lt_sem(y');
    %     shadedErrorBar(t, ymean, ysem, {'Color', pcols{i}}, 1);
    shadedErrorBar(t, ymean, ysem, {'Color', 'r'}, 1);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
    
    
    % -------- FINAL
    datWN = squeeze(allDat_wn(:,:, ind3, :));
    ymean = mean(datWN, 2);
    ysem = lt_sem(datWN');
    plot(t, ymean, 'Color', 'm');
    
    % -- fomrat
    axis tight;
    xlim(XLIM);
    lt_plot_zeroline_vert;
    
end



%% ===== PLOT EACH EXPERIMENT, each syl

pcols = lt_make_plot_colors(size(allDat_epochs,1), 1, [1 0 0]);
t = PARAMS.Xcov_ccLags;

figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


blist = unique(allbnum);
for jj=1:length(blist)
    birdthis = blist(jj);
    indsthis = allbnum==birdthis;
    
% =======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG');
xlabel('lag (LMAN-RA)');
ind3 = 1; % target

% ---- 1) plot baseline
dat = squeeze(allDat_base(:,:, ind3, indsthis));
ymean = mean(dat, 2);
plot(t, ymean, ':', 'Color', 'k', 'LineWidth', 2);

% ----- 2) WN epochs
dat = squeeze(allDat_epochs(:, :, ind3, indsthis));
for i=1:size(dat,1)
    y = squeeze(dat(i,:,:));
    ymean = mean(y,2);
    plot(t, ymean, 'Color', pcols{i}, 'LineWidth', 2);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
end

% % --- 3) plot WN
% dat = squeeze(allDat_wn(:,:, ind3, :));
% ymean = mean(dat, 2);
% plot(t, ymean, ':', 'Color', 'm', 'LineWidth', 2);

% -- fomrat
axis tight;
xlim(XLIM);
lt_plot_zeroline_vert;


% =======================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
% title('TARG');
xlabel('lag (LMAN-RA)');

% ---- 1) get baseline
datBase = squeeze(allDat_base(:,:, ind3, indsthis));

% ----- 2) WN epochs
dat = squeeze(allDat_epochs(:, :, ind3, indsthis));
for i=1:size(dat,1)
    y = squeeze(dat(i,:,:));
    y = y - datBase;
    ymean = mean(y,2);
    plot(t, ymean, 'Color', pcols{i}, 'LineWidth', 2);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
end

% % --- 3) plot WN
% dat = squeeze(allDat_wn(:,:, ind3, :));
% ymean = mean(dat, 2);
% plot(t, ymean, ':', 'Color', 'm', 'LineWidth', 2);

% -- fomrat
axis tight;
xlim(XLIM);
lt_plot_zeroline_vert;
lt_plot_zeroline;



% =================== PLOT EACH EPOCH SEPARATELY
% ------- RUN
dat = squeeze(allDat_epochs(:, :, ind3, indsthis));
for i=1:size(dat,1)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
title('TARG');
    xlabel('lag (LMAN-RA)');
    ylabel(['epoch(WN): ' num2str(i)]);
    
    
    % ------------ baseline
    datbase = squeeze(allDat_base(:,:, ind3, indsthis));
    ymean = mean(datbase, 2);
    ysem = lt_sem(datbase');
    shadedErrorBar(t, ymean, ysem, {'Color', 'k', 'LineStyle', '--'}, 1);
    
    
    % -------- epoch (this)
    y = squeeze(dat(i,:,:));
    ymean = mean(y,2);
    ysem = lt_sem(y');
    %     shadedErrorBar(t, ymean, ysem, {'Color', pcols{i}}, 1);
    shadedErrorBar(t, ymean, ysem, {'Color', 'r'}, 1);
    lt_plot_text(t(end), ymean(end), ['epoch' num2str(i)], pcols{i}, 9);
    
    
    % -------- FINAL
    datWN = squeeze(allDat_wn(:,:, ind3, indsthis));
    ymean = mean(datWN, 2);
    ysem = lt_sem(datWN');
    plot(t, ymean, 'Color', 'm');
    
    % -- fomrat
    axis tight;
    xlim(XLIM);
    lt_plot_zeroline_vert;
    
end

end

% =========
linkaxes(hsplots, 'xy');
