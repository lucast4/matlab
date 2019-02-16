function [datbase, datWN, datWN_notminshuff, datbase_notminshuff, Xq, NanCount] = ...
    lt_neural_POPLEARN_XCov_sub1(datmat_real, datmat_shuff, dosmooth, ...
    dosmooth_sigma, inds_base, inds_WN, xbins, plotraw, xcovver)
if ~exist('plotraw', 'var')
    plotraw = 0;
end
%% ==================== SMOOTH?
if dosmooth==1
    % SMOOTH PARAMS (gaussian)
    N = 2*dosmooth_sigma/0.001;
    win = gausswin(N);
    win = win./sum(win);
    % INTERP PARAMS.
    Xq = xbins(1):0.001:xbins(end);
    
    % 1) INTERPOLATE
    try
        tmp = interp1(xbins, datmat_real', Xq);
    catch err
        tmp = interp1(lags_sec, datmat_real', Xq);
    end
    
    datmat_real = tmp';
    % 2) SMOOTH
    for kkk=1:size(datmat_real,1)
        datmat_real(kkk,:) = conv(datmat_real(kkk, :), win, 'same');
    end
    
    % 1) INTERPOLATE
    tmp = interp1(xbins, datmat_shuff', Xq);
    datmat_shuff = tmp';
    % 2) SMOOTH
    for kkk=1:size(datmat_shuff,1)
        datmat_shuff(kkk,:) = conv(datmat_shuff(kkk, :), win, 'same');
    end
else
    Xq = [];
end

%% ============= NOTE DOWN WHAT FRACTION OF TRIALS ARE NAN

Nnan = sum(any(isnan(datmat_real)'));
Ntot = length(any(isnan(datmat_real)'));

NanCount = [Nnan Ntot];
%% get xcov
% ================ 2versions, either using SHIFT PREDICTOR CACLULATE IN 2
% DIRECTIONS OR IN ONE.
if size(datmat_shuff,1) == 2*(size(datmat_real,1)-1)
    % ============= NEW VERSION, WITH SHUFFLE DONE IN BOTH DIRECTIOSN
    % I.E. N VS N-1 AND N VS N+1. shuffle has more trials (N=2*(Norig-1)). need
    % to match indices of shuffle back to original data. do this in follwoing
    % way: each data trial is matched to the shuffle trials that are n vs. n+1
    % and n vs. n-1.
    % note: alterantive is to be slightly more conservative, only shuflfles for
    % which both trials are within the WN window. i.e. window is: [2(n1-1)+1
    % 2(n2-1)] where n1 and n2 are the first and last data trials. This is
    % difficult when there are skipped trials... easier to just do this.
    
    % --- find the correct shuffle trials
    tmp = sort([1 2:size(datmat_real,1)-1 2:size(datmat_real,1)-1 size(datmat_real,1)]); % tmp(2) tells you the ind (in real dat) that shuff trial 2 corresponds to
    assert(length(tmp)==size(datmat_shuff,1));
    
    inds_base_shuff = find(ismember(tmp, inds_base));
    inds_WN_shuff = find(ismember(tmp, inds_WN));
elseif size(datmat_shuff,1)==size(datmat_real,1)
    % ======== OLD VERSION, SHUFF INDICES ARE EXACTYL SAME AS INDS FOR REAL
    % DAT
    inds_base_shuff = inds_base;
    inds_WN_shuff = inds_WN;
    
else
    disp('doesnt make sense, not sure how shift predictors were cpmputed...')
    pause;
end
% datbase = mean(datmat_real(inds_base,:)) ...
%     - mean(datmat_shuff(inds_base,:));
%
% datWN = mean(datmat_real(inds_WN,:)) ...
%     - mean(datmat_shuff(inds_WN,:));
%
% datWN_notminshuff = mean(datmat_real(inds_WN,:));
% datbase_notminshuff = mean(datmat_real(inds_base,:));

datbase = nanmean(datmat_real(inds_base,:)) ...
    - nanmean(datmat_shuff(inds_base_shuff,:));

datWN = nanmean(datmat_real(inds_WN,:)) ...
    - nanmean(datmat_shuff(inds_WN_shuff,:));

datWN_notminshuff = nanmean(datmat_real(inds_WN,:));
datbase_notminshuff = nanmean(datmat_real(inds_base,:));


% ################## OTHER VERSIONS
% ====== 1) Z-SCORE
% for each lag, get distribution of values over trials. use that to
% z-transform actual data
if strcmp(xcovver, 'zscore')
    % 1) BASE
    yshuff = datmat_shuff(inds_base_shuff,:);
    ydat = datmat_real(inds_base,:);
    
    yz_base = (ydat - nanmean(yshuff,1))./nanstd(yshuff, [], 1);
    
    
    % 2) WN
    yshuff = datmat_shuff(inds_WN_shuff,:);
    ydat = datmat_real(inds_WN,:);
    
    yz_WN = (ydat - nanmean(yshuff,1))./nanstd(yshuff, [], 1);
    
    % ========== REPLACE OUTPUT
    datbase = nanmean(yz_base);
    datWN = nanmean(yz_WN);
elseif strcmp(xcovver, 'scale')
    % ======== 2) MATCH THE MAGNITUDE OF SHUFFLES DURING WN VS. BASELINE
%     multfact = nanmean(nanmean(datmat_shuff(inds_base_shuff,:)))/nanmean(nanmean(datmat_shuff(inds_WN_shuff,:))); % one mult for all lags
    multfact = nanmean(datmat_shuff(inds_base_shuff,:))./nanmean(datmat_shuff(inds_WN_shuff,:)); % one for each lag
    datWN_scaled = multfact.*datWN;
end

if plotraw==1
    %% =========== sanity check, plotting raw xcov.
    
    lt_figure; hold on;
    lt_subplot(3,2,1:2); hold on;
    title('k=base; r=wn [dash = base]');
    
    % === base shuff
    dat = datmat_shuff(inds_base_shuff,:);
    y = mean(dat);
    ysem = lt_sem(dat);
    x = 1:length(y);
    shadedErrorBar(x, y, ysem, {'Color', 'k', 'LineStyle', '--'}, 1);
    
    % === base dat
    dat = datmat_real(inds_base,:);
    y = mean(dat);
    ysem = lt_sem(dat);
    x = 1:length(y);
    shadedErrorBar(x, y, ysem, {'Color', 'k', 'LineStyle', '-'}, 1);
    
    % === wn shuff
    dat = datmat_shuff(inds_WN_shuff,:);
    y = mean(dat);
    ysem = lt_sem(dat);
    x = 1:length(y);
    shadedErrorBar(x, y, ysem, {'Color', 'r', 'LineStyle', '--'}, 1);
    
    
    % === wn dat
    dat = datmat_real(inds_WN,:);
    y = mean(dat);
    ysem = lt_sem(dat);
    x = 1:length(y);
    shadedErrorBar(x, y, ysem, {'Color', 'r', 'LineStyle', '-'}, 1);
    
    lt_plot_zeroline;
    
    lt_subplot(3,2,3); hold on;
    title('z-scored(mean of)');
    plot(mean(yz_base), '-k');
    plot(mean(yz_WN), '-r');
    
    lt_subplot(3,2,4); hold on;
    title('z-scored(median)');
    plot(median(yz_base), '-k');
    plot(median(yz_WN), '-r');
    
    lt_subplot(3,2,5); hold on;
    title('dat - shuff [k=base, r=wn; m=wn(scaled)]');
    plot(datbase, '-k');
    plot(datWN, '-r');
    plot(datWN_scaled, '-m');
    
    % === plot change over time
    lt_subplot(3,2,6); hold on;
    x = 1:size(datmat_shuff(:,12),1);
    plot(x, datmat_shuff(:,12), 'ok');
    plot([x(1) x(2:2:end-1) x(end)], datmat_real(:,12), 'or');
end