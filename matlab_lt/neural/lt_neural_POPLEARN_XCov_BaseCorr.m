function lt_neural_POPLEARN_XCov_BaseCorr(OUTSTRUCT_XCOV, PARAMS, onlygoodexpt, SwitchStruct)
%% lt 2/27/19 - plots sliding dot product (baseline)


%% only start of experimetns?
if onlygoodexpt==1
    % ===== filter outstruct
    expttype = 'xcov_spikes';
%     [OUTSTRUCT] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct, expttype);
    [OUTSTRUCT_XCOV] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, SwitchStruct, expttype);
end

%%
% === lag window
% windtoplot = [-0.035 0.035];
windtoplot = [-0.04 0.04];

% =========== get matrix
Yxcorr = OUTSTRUCT_XCOV.XcovBase_NoMinShuff;
txcorr = PARAMS.Xcov_ccLags;

% ---
indx = txcorr>=windtoplot(1) & txcorr<=windtoplot(2);
Yxcorr = Yxcorr(:, indx);
txcorr = txcorr(indx);

% ---
maxbird = max(OUTSTRUCT_XCOV.bnum);
%%
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ==== all data
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all dat');
xlabel('LMAN -- RA');
ylabel('sliding dot prod');
Y = Yxcorr';
x = txcorr;
% --- subtract mean for each case
Y = Y-mean(Y,1);
x = txcorr;

Ystd = std(Y, [], 2);
Ymean = mean(Y, 2);
Ysem = lt_sem(Y');
shadedErrorBar(x, Ymean, Ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, Ymean, Ysem, {'Color', 'k'}, 1);
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('all dat');
xlabel('LMAN -- RA');
ylabel('sliding dot prod (mean centered)');
Y = Yxcorr';
x = txcorr;
% --- subtract mean for each case
Y = Y-mean(Y,1);

% --- sort by max
[~, indmax] = max(Y);
[~, indsort] = sort(indmax);
Y = Y(:, indsort);
clim = [-prctile(abs(Y(:)), 95) prctile(abs(Y(:)), 95)];
lt_neural_Coher_Plot(Y, x, 1:size(Y,2),1, '', clim);


% ==== all data
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('each bird');
xlabel('LMAN -- RA');
ylabel('sliding dot prod');
pcols = lt_make_plot_colors(maxbird, 0, 0);
for i=1:maxbird
    indsthis = OUTSTRUCT_XCOV.bnum==i;
    if ~any(indsthis)
        continue
    end
    Y = Yxcorr(indsthis, :)';
    x = txcorr;
    % --- subtract mean for each case
    Y = Y-mean(Y,1);
    x = txcorr;
    
    Ystd = std(Y, [], 2);
    Ymean = mean(Y, 2);
    Ysem = lt_sem(Y');
    shadedErrorBar(x, Ymean, Ystd, {'Color', 0.3*pcols{i}}, 1);
    shadedErrorBar(x, Ymean, Ysem, {'Color', pcols{i}}, 1);
    lt_plot_zeroline_vert;
end


% ==== all data
for i=1:maxbird
    indsthis = OUTSTRUCT_XCOV.bnum==i;
    if ~any(indsthis)
        continue
    end
    
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title(['bird ' num2str(i)]);
xlabel('LMAN -- RA');
ylabel('sliding dot prod');

Y = Yxcorr(indsthis, :)';
    x = txcorr;
    % --- subtract mean for each case
    Y = Y-mean(Y,1);
    x = txcorr;
    
    % --- sort by max
    [~, indmax] = max(Y);
    [~, indsort] = sort(indmax);
    Y = Y(:, indsort);
    lt_neural_Coher_Plot(Y, x, 1:size(Y,2),1, '', clim);
end