function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Moment(DATBYREND, TrialStruct)
%% lt 9/30/18 - moment to moment learning, predicted by catch and FF.

%% ======================== EXTRACT RELEVANT EXPERIMENTS

getsiglearn = 0;
% birdstoget = 13:17;
% birdstoget = [13 17];
birdstoget = 'notSDP';
syltype = 'same';
% syltype = 'same';
% syltype = 'all';
% syltype = 'diff';
% syltype = 'nontarg';
minbaserends = 0; % has at least this many baseline renditions.
minrends = 20;

[Inds_sylcounter, Inds_birdnum, Inds_exptnum, ...
    Inds_IsSame, Inds_IsTarg] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
    syltype, getsiglearn, birdstoget, minbaserends, minrends);


%% ######################## PREDICTING LEARNING BY LOCAL LEARNING, HIT, AND CATCH
DATTMP = struct;
% NOTE: DO NOT CHANGE ORDER!! of following 3 extractions...
use_targ_locallearn = 1;
use_nminus2 = 0; % if 1, then use n-2 as local learning. if 0. thne use n-1.

if ~strcmp(syltype, 'targ')
    % then should use minus 1, so deviations will be agaiunst your own syl
    use_nminus2 = 0;
end

% ------------------ TRAINING, NOTCATCH
istrain = 1;
iscatch = 0;
DATTMP.Train = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
    TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn);

% ------------------ TRAINING, CATCH
istrain = 1;
iscatch = 1;
DATTMP.Catch = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
    TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn);

% ------------------ BASELINE, C/NC
istrain = 0;
iscatch = [0 1];
DATTMP.Base = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
    TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn);

% ------------------ use the local learning for the syllable of interest,
% or for the target syllable in that experiment?



%% ########################### DIFFERENCE BETWEEN TRAIN AND CATCH FOR TRIAL N+1?
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

binsize = 10;
timedevwind = [2 60]; % minuites

% ======= collect for one large analysis
xallall = [];
yallall = [];
timedevallall = [];
iscatchall = [];
exptcounter = [];
count = 1;
for j=1:length(DATTMP.Train)
    
    if length(DATTMP.Catch)<j
        % then this expt doesn't have catch data ...
        continue
    end
    
    % ========== training
    x = DATTMP.Train(j).learnlocal; % (trial n minus n-1)
%     y = DATTMP.Train(j).ffdev_first+DATTMP.Train(j).learnlocal; % (trial n+1 minus n)
    y = DATTMP.Train(j).ffdev_first; % (trial n+1 minus n)
    timedev = DATTMP.Train(j).tdev_first;
    
    % -- remove nans
    indtmp = ~isnan(x);
    timedev = timedev(indtmp);
    y = y(indtmp);
    x = x(indtmp);
    
    indtmp = timedev>timedevwind(1) & timedev<timedevwind(2);
    x = x(indtmp);
    y = y(indtmp);
    timedev = timedev(indtmp);
    
        % ========== catch
    xcatch = DATTMP.Catch(j).learnlocal;
%     ycatch = DATTMP.Catch(j).ffdev_first+DATTMP.Catch(j).learnlocal;
    ycatch = DATTMP.Catch(j).ffdev_first;
    timedevcatch = DATTMP.Catch(j).tdev_first;
    
    indtmp = ~isnan(xcatch);
    timedevcatch = timedevcatch(indtmp);
    ycatch = ycatch(indtmp);
    xcatch = xcatch(indtmp);
    
    indtmp = timedevcatch>timedevwind(1) & timedevcatch<timedevwind(2);
    xcatch = xcatch(indtmp);
    ycatch = ycatch(indtmp);
    timedevcatch = timedevcatch(indtmp);

    % ========= CONCATENATE INTO ONE ARRAY
    xall = [x; xcatch];
    yall = [y; ycatch];
    timedevall = [timedev; timedevcatch];
    iscatch = [zeros(size(x)); ones(size(xcatch))];
    
    
    % ========= PLOT
    % -- 1) overlay both data with linear regression
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    %         title([bname '-' ename ',' fname]);
    xlabel('local learn (n minus n-1)');
    ylabel('ff dev (n+1 minus n-1)');
    
    if (0)
    % ----- plot each all trials with lin reg overlay
    lt_regress(y, x, 1, 0, 1, 1, 'k');
    lt_regress(ycatch, xcatch, 1, 0, 1, 1, 'b');
    else
    % ----- overlay smoothed and datapoints
    plot(x,y, 'kx');
    plot(xcatch, ycatch, 'bx');
    
    [~, indtmp] = sort(x);
    ysort = y(indtmp);
    xsort = x(indtmp);   
    tmpy = lt_running_stats(ysort, binsize);
    tmpx = lt_running_stats(xsort, binsize);
    shadedErrorBar(tmpx.Mean, tmpy.Mean, tmpy.SEM, {'Color', 'k'},1);
    
    [~, indtmp] = sort(xcatch);
    ysort = ycatch(indtmp);
    xsort = xcatch(indtmp);   
    tmpy = lt_running_stats(ysort, binsize);
    tmpx = lt_running_stats(xsort, binsize);
    shadedErrorBar(tmpx.Mean, tmpy.Mean, tmpy.SEM, {'Color', 'b'},1);
    end
    
    % -- 2) overlay distributions for all data
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('local learn (n minus n-1)');
    [~, xtmp] = lt_plot_histogram(x, [], 1, 1, '', 1, 'k');
    lt_plot_histogram(xcatch, xtmp, 1, 1, '', 1, 'b');
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title('ffdev (n+1 minus n-1)');
    [~, xtmp] = lt_plot_histogram(y, [], 1, 1, '', 1, 'k');
    lt_plot_histogram(ycatch, xtmp, 1, 1, '', 1, 'b');
    
    % --- 3) Overlay smoothed curves
    
    
    % ============ COLLECT ACROSS EXPT
    xallall = [xallall; xall];
    yallall = [yallall; yall];
    timedevallall = [timedevallall; timedevall];
    iscatchall = [iscatchall; iscatch];
    exptcounter = [exptcounter; count*ones(size(xall))];
    count = count+1;
end

lt_figure; hold on;
lt_regress(yallall(iscatchall==0), xallall(iscatchall==0), 1, 0, 1, 1, 'k');
lt_regress(yallall(iscatchall==1), xallall(iscatchall==1), 1, 0, 1, 1, 'b');


% =================
inds = timedevallall<2 & iscatchall==0;
mean(yallall(inds))
inds = timedevallall<2 & iscatchall==1;
mean(yallall(inds))
inds = timedevallall>2 & iscatchall==0;
mean(yallall(inds))
inds = timedevallall>2 & iscatchall==1;
mean(yallall(inds))


lt_figure; hold on;
[Ybyexpt_wn Ybyexpt_wn_sem] = grpstats(yallall(iscatchall==0), exptcounter(iscatchall==0), {'mean', 'sem'});
[Ybyexpt_catch Ybyexpt_catch_sem] = grpstats(yallall(iscatchall==1), exptcounter(iscatchall==1), {'mean', 'sem'});
xlabel('WN');
ylabel('catch')
title('trial(n+1)-trial(n)');
% plot(Ybyexpt_wn, Ybyexpt_catch, 'ok');
% lt_plot_makesquare_plot45line(gca, 'b', -10);
lt_plot_45degScatter(Ybyexpt_wn, Ybyexpt_catch, 'k', 1, 1, Ybyexpt_wn_sem, Ybyexpt_catch_sem);


% ========================= LME, account for 
% y1 = AllPairs_Means(:,1);
% y2 = AllPairs_Means(:,2);
tbl = table(yallall, xallall, iscatchall, exptcounter);

mdl = 'yallall ~ xallall + iscatchall + (1|exptcounter) + (-1 + xallall|exptcounter) + (-1 + iscatchall|exptcounter)';  % full model
mdl = 'yallall ~ xallall + iscatchall + (1|exptcounter) + (-1 + xallall|exptcounter)';  % full model

lme = fitlme(tbl, mdl)





%% ########################### PLOT (absolute balue of local earning - i.e.
% more dwivation from mmedian predict more laernign?)

figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

takeabs = 1; % takes absolute value of local learn (i.e. might predict strong increase in next
% trial when current trial is more deviated from local.
fnames = fieldnames(DATTMP);
for j=1:length(DATTMP.Train)
    
    for fname = fnames'
        fname = fname{1};
        % ================
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        %         title([bname '-' ename ',' fname]);
        xlabel('local learn (n minus n-1)');
        ylabel('ff dev (n+1 minus n-1)');
        
        try
            if use_targ_locallearn==1
                x = DATTMP.(fname)(j).learnlocal;
            else
                x = DATTMP.(fname)(j).learnlocal;
            end
            
            if use_nminus2==1
                y = DATTMP.(fname)(j).ffdev_first + DATTMP.(fname)(j).learnlocal;
            else
                y = DATTMP.(fname)(j).ffdev_first;
            end
            
            if takeabs==1
                x = abs(x);
            end
            plot(x, y, 'x');
            lt_plot_makesquare_plot45line(gca, 'b');
            lt_regress(y, x, 0, 0, 1, 1, 'r', 1);
        catch err
        end
    end
end

linkaxes(hsplots, 'xy');



%% ==================================== EXTRACT SLOPES (TRAIN, CATCH, BASE)
% =========== 1) COLLECT
plotON = 1;
plotWNhits = 1;
[SlopesAll, Learn_And_FFdevAll] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_7(DATTMP, ...
    use_nminus2, plotON, [], plotWNhits);


% ########################### PLOT (separately plots positive and negative
% deviations)
if use_nminus2==0
    lt_figure;
    lt_plot_text(0, 0.5, 'NOTE: using (n-2)-(n-1) for ffdev');
end

% ========== SUMMARY PLOT
lt_figure; hold on;
excludeIfMissCatch = 1;

% =========
lt_subplot(3,2,1); hold on;
xlabel('NEG(train) - NEG(catch) | POS(train) - POS(catch)');
ylabel('slope (trial n+2 vs trial n+1)');

Y = SlopesAll(:, [1 2 4 5]);
if excludeIfMissCatch==1
    Y = Y(~any(isnan(Y)'), :);
end
x = 1:size(Y,2);
plot(x,Y, '-ok');
xlim([0 5]);
lt_plot_zeroline;

% ----------
lt_subplot(3,2, 2); hold on;
xlabel('NEG -- POS -- AVERAGE(POS,NEG)');
ylabel('slope (subtract catch, pos = more learning than catch)');
Ythis = [Y(:,1)-Y(:,2) Y(:,3)-Y(:,4)];
Ythis(:,1) = -Ythis(:,1); % flip so that positive is more learning than catch.
Ythis = [Ythis mean(Ythis, 2)]; % take average for 3rd column.
x = 1:size(Ythis,2);
plot(x, Ythis, '-ok');
xlim([0 3]);
lt_plot_zeroline;

% ============ REPLOT ALL DATA, USING MEASURE WHERE UP IS MORE LEARNING
Yadaptive = [1-SlopesAll(:,1:3) SlopesAll(:,4:6)];

% --------------- 1) all slopes
lt_subplot(3,2,3); hold on;
excludeIfMissCatch = 1;
xlabel('NEG(train) - NEG(catch) - NEG(base) | POS(train) - POS(catch) - POS(catch)');
ylabel('slope (trial n+2 vs trial n+1) [adaptive dir]');

Y = Yadaptive;
if excludeIfMissCatch==1
    Y = Y(~any(isnan(Y)'), :);
end
x = 1:size(Y,2);
for j=1:size(Y,1)
    plot(x, Y(j,:), '-ok');
end
xlim([0 7]);

% ---------------- 2) MEAN DEVIATIONS
lt_subplot(3,2,4); hold on;
excludeIfMissCatch = 1;
xlabel('NEG(train) - NEG(catch) - NEG(base) | POS(train) - POS(catch) - POS(catch)');
ylabel('mean ffdev (trial n+2 minus n)');

Y = Learn_And_FFdevAll;
if excludeIfMissCatch==1
    Y = Y(~any(isnan(SlopesAll)'), :);
end

functmp = @(x)(mean(x(:,2)));
FFdev_mean = cell2mat(cellfun(functmp, Y, 'UniformOutput', 0));

x = 1:size(FFdev_mean,2);
plot(x, FFdev_mean, '-ok');
xlim([0 7]);
lt_plot_zeroline;



%% =============== COLLECT SLOPES, SEPARATING BY TIME BIN OF DEVIATION
tbinedges = [2]; % minutes, will fill in the edges

% ========== COLLECT MULTIPLE TIME BINS AND PLOT
plotON = 0;
twind = [2 60];
[SlopesAll, Learn_And_FFdevAll] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_7(DATTMP, ...
    use_nminus2, plotON, twind);


%% ============= CROSS CORRELATION OF DEVIATIONS

maxwind = 50;
fnames = fieldnames(DATTMP);

CCAll = cell(length(DATTMP.Train), 3); % [syls x [train, catch, base]];
for j=1:length(DATTMP.Train)
    rowcount = 1;
    for fname = fnames'
        fname = fname{1};
        
        if length(DATTMP.(fname))>=j
            y1 = DATTMP.(fname)(j).learnlocal;
            y2 = DATTMP.(fname)(j).ffdev_first;
            
            if any(isnan(y1)) | any(isnan(y2))
                assert(sum(isnan(y1))/length(y1) < 0.2, 'too much nans...');
                assert(sum(isnan(y2))/length(y2) < 0.2, 'too much nans...');
                
                indstmp = any(isnan([y1 y2])');
                y1(indstmp) = [];
                y2(indstmp) = [];
            end
            
            [cc, lags] = xcov(y1, y2, 20, 'coeff');
            
            
        else
            cc = [];
            lags = [];
        end
        
        % ============
        CCAll{j, rowcount} = cc;
        rowcount = rowcount+1;
    end
end

% ============== IMAPORTANT
lt_figure;
lt_plot_text(0, 0.5, 'NOTE: do not plot catch, since they are not temporally sequenced...');


% ================ plot
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

onlykeepfull=1;
if onlykeepfull==1
    CCAll(any(cellfun(@isempty, CCAll)'), :) = [];
end


% ====== 1) overlay all xcovs
pcollist = {'r', 'm', 'k'};
for j=1:size(CCAll,1)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots =[hsplots hsplot];
    
    % ---- train
    indtmp = 1;
    pcol = 'r';
    
    cc = CCAll{j, indtmp};
    plot(lags, cc, '-', 'Color', pcol);
    
    %     % ---- catch
    %     indtmp = 2;
    %     pcol = 'm';
    %
    %     cc = CCAll{j, indtmp};
    %     plot(lags, cc, '-', 'Color', pcol);
    
    % ---- base
    indtmp = 3;
    pcol = 'k';
    
    cc = CCAll{j, indtmp};
    plot(lags, cc, '-', 'Color', pcol);
    
    % ======
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

linkaxes(hsplots, 'xy');

lt_figure; hold on;
title('r=train; k=base');
CCAlltmp = cellfun(@transpose, CCAll, 'UniformOutput', 0);
for i=[1 3]
    lt_subplot(2,2,i); hold on;
    xlabel('local laerning (trial n+1 minus n)');
    ylabel('ffdev (usually trial n+2 minus trial n+1)');
    
    pcol = pcollist{i};
    
    ccall = cell2mat(CCAlltmp(:, i));
    %
    %     ccall = reshape(cell2mat(CCAll(:,i)), size(CCAll,1), []);
    
    plot(lags, ccall, '-k');
    ymean = mean(ccall,1);
    ysem = lt_sem(ccall);
    x = 1:length(ymean);
    
    lt_plot(lags, ymean, {'Errors', ysem, 'Color', pcol});
    
    % ======
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

% ======== subtract  base from train
ccall_train = cell2mat(CCAlltmp(:, 1));
ccall_base = cell2mat(CCAlltmp(:, 3));

ccnorm = ccall_train - ccall_base;

lt_subplot(2,2,4); hold on;
xlabel('local laerning (trial n+1 minus n)');
ylabel('ffdev (usually trial n+2 minus trial n+1)');
title('train minus base');
pcol = 'b';

plot(lags, ccnorm, '-k');
ymean = mean(ccnorm,1);
ysem = lt_sem(ccnorm);
x = 1:length(ymean);

lt_plot(lags, ymean, {'Errors', ysem, 'Color', pcol});









