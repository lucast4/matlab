function OUTSTRUCT = lt_opto_ExtrBirdDat(dirtoplot, twind, SwitchTimes, ...
    plotON)
%% lt 10/11/18 - extract dat for a bird. have to first run Opto analy stuff...
% INPUTS:
% --- BU6
% dirtoplot = '/bluejay5/lucas/birds/bu6bu98/Opto_Stim_analy/Reversion1';
% twind = 1;
%
% SwitchTimes = {...
%     '29Sep2018-2021', 'of-up', ...
%     '04Oct2018-2400', 'up-dn', ...
%     '09Oct2018-2400', 'dn-up'...
%     };


if ~exist('plotON', 'var')
    plotON = 1;
end

%% PARAMS

% =========== LOAD DATA
load([dirtoplot '/PARAMS.mat']);
load([dirtoplot '/DATSTRUCT.mat']);

% ========= PLOTCOLS
clear pcolstruct
pcolstruct.All = [0.3 0.3 0.3];
pcolstruct.StimNotCatch = 'b';
pcolstruct.StimCatch = [0.7 0.2 0.2];

%% ========= CONVERT DAT TIME TO DATENUMS
% =========
numdays = length(DATSTRUCT.data);
assert(PARAMS.global.LastDate_num-PARAMS.global.FirstDate_num+1 == numdays, 'then my assumption about inds equaling days is incorrect');

%
firstday_datenum = PARAMS.global.FirstDate_num;
for i=1:numdays
    if isempty(DATSTRUCT.data{i})
        continue
    end
    fieldsthis = fieldnames(DATSTRUCT.data{i}.timewindow{twind});
    for fthis = fieldsthis'
        tvals = DATSTRUCT.data{i}.timewindow{twind}.(fthis{1}).timevals;
        
        tvals_dnum = firstday_datenum+(i-1)+(tvals/24);
        DATSTRUCT.data{i}.timewindow{twind}.(fthis{1}).timevals_dnum = tvals_dnum;
    end
end


%% =================== PLOT RAW DATA [all data for this bird]
if plotON==1;
    lt_figure; hold on;
    for i=1:numdays
        if isempty(DATSTRUCT.data{i})
            continue
        end
        fieldsthis = fieldnames(DATSTRUCT.data{i}.timewindow{twind});
        
        for fthis = fieldsthis'
            if (0)
                tvals = DATSTRUCT.data{i}.timewindow{twind}.(fthis{1}).timevals;
                tvals = (24*(i-1))+tvals;
            else
                tvals = DATSTRUCT.data{i}.timewindow{twind}.(fthis{1}).timevals_dnum;
            end
            ffvals = DATSTRUCT.data{i}.timewindow{twind}.(fthis{1}).ffvals;
            
            % ======= PLOT
            % - raw
            pcol = pcolstruct.(fthis{1});
            plot(tvals, ffvals, 'x', 'Color', pcol);
            % - mean
            x = ceil(tvals(end));
            lt_plot(x, mean(ffvals), {'Errors', lt_sem(ffvals), 'Color', pcol});
        end
    end
    
    %%
    if ~isempty(SwitchTimes)
    % ========= CONVERT SWITCH TIMES
    [SwitchTimesMat, SwitchTimesMat_conting] = ...
        lt_tools_SwitchTimes(SwitchTimes, PARAMS.global.LastDate_num+1);
    
    % ====== put lines for swich times
    for i=1:length(SwitchTimes)/2
        dtime = datenum(SwitchTimes{(i*2-1)}, 'ddmmmyyyy-HHMM');
        conting = SwitchTimes{(i*2)};
        line([dtime dtime], ylim);
        YLIM = ylim;
        lt_plot_text(dtime, YLIM(2), conting, 'k');
    end
    
    for i=1:size(SwitchTimesMat,1)
        sw1 = SwitchTimesMat(i,1);
        if sw1==0
            XLIM = xlim;
            sw1 = XLIM(1);
        end
        sw2 = SwitchTimesMat(i,2);
        cont = SwitchTimesMat_conting(i);
        
        pcol = [rand rand rand];
        line([sw1 sw2], [YLIM(2) YLIM(2)], 'Color', pcol);
        lt_plot_text(mean([sw1 sw2]), YLIM(2), num2str(cont), pcol);
    end
    end
end

% --- date tick
datetick('x', 'mm/dd');
rotateXLabels(gca, 90);

%% ########### SUMMARY FOR AFP BIAS AS FUNCTION OF DIRECTION OF TRAINIG  [ExTRACT]
OUTSTRUCT.All_Dprime = [];
OUTSTRUCT.All_FFmeanStimMinusNostim = [];
OUTSTRUCT.All_FFmedianStimMinusNostim = [];
OUTSTRUCT.All_pvalranksum = [];
OUTSTRUCT.All_traindir = [];
if ~isempty(SwitchTimes)
for i=1:numdays
    
    if isempty(DATSTRUCT.data{i})
        continue
    end
    
    fieldsthis = fieldnames(DATSTRUCT.data{i}.timewindow{twind});
    
    if ~(any(strcmp(fieldsthis, 'StimNotCatch')) & any(strcmp(fieldsthis, 'StimCatch')))
        % then this day doesn't have stim and catch.
        continue
    end
    
    f2 = DATSTRUCT.data{i}.timewindow{twind}.StimNotCatch.ffvals;
    t2 = DATSTRUCT.data{i}.timewindow{twind}.StimNotCatch.timevals_dnum;
    
    f1 = DATSTRUCT.data{i}.timewindow{twind}.StimCatch.ffvals;
    t1 = DATSTRUCT.data{i}.timewindow{twind}.StimCatch.timevals_dnum;
    
    % -------- check that distributions for tvals are similar
    assert(ranksum(t1, t2)>0.01, 'are time values different? then is problem as this assumes interleaves stim/nostim');
    
    % --- calculate dprime
    dprime = (mean(f2)-mean(f1))/(sqrt(0.5*(var(f1)+var(f2))));
    ffmeandiff_stim = mean(f2)-mean(f1);
    ffmedian = median(f2) - median(f1);
    p = ranksum(f2, f1);
    
    % --- what direction learning today?
    epochthis = all(([t1 t2]>SwitchTimesMat(:,1))') & ...
        all(([t1 t2]<SwitchTimesMat(:,2))'); % find the current "epoch"
    assert(sum(epochthis)==1, 'means trials in this dat span 2 epochs?');
    learndir = SwitchTimesMat_conting(epochthis);
    
    % =================== OUTPUT
    OUTSTRUCT.All_Dprime = [OUTSTRUCT.All_Dprime; dprime];
    OUTSTRUCT.All_FFmeanStimMinusNostim = [OUTSTRUCT.All_FFmeanStimMinusNostim; ffmeandiff_stim];
    OUTSTRUCT.All_FFmedianStimMinusNostim = [OUTSTRUCT.All_FFmedianStimMinusNostim; ffmedian];
    OUTSTRUCT.All_pvalranksum = [OUTSTRUCT.All_pvalranksum; p];
    OUTSTRUCT.All_traindir = [OUTSTRUCT.All_traindir; learndir];
end
end

%% ######### GET ORIGINAL DATA FOR EACH DAY

for i=1:numdays
    
    if isempty(DATSTRUCT.data{i})
        continue
    end
    
    
    
end
%% =============== PLOT
if plotON==1
    lt_figure; hold on;
    
    % --------------- SUMMARY
    lt_subplot(3,2,1); hold on;
    xlabel('train dir');
    ylabel('stim minus nostim (ff)');
    y = OUTSTRUCT.All_FFmeanStimMinusNostim;
    
    % - run
    plot(OUTSTRUCT.All_traindir, y, 'ok');
    x = unique(OUTSTRUCT.All_traindir);
    [ymean, ysem] = grpstats(y, OUTSTRUCT.All_traindir, {'mean', 'sem'});
    lt_plot(x+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
    [h, p] = ttest2(y(OUTSTRUCT.All_traindir==-1), y(OUTSTRUCT.All_traindir==1));
    lt_plot_pvalue(p, '-1 vs 1', 1);
    lt_plot_zeroline;
    
    % --------------- SUMMARY
    lt_subplot(3,2,2); hold on;
    xlabel('train dir');
    ylabel('dprime (stim minus base)');
    y = OUTSTRUCT.All_Dprime;
    
    % - run
    plot(OUTSTRUCT.All_traindir, y, 'ok');
    x = unique(OUTSTRUCT.All_traindir);
    [ymean, ysem] = grpstats(y, OUTSTRUCT.All_traindir, {'mean', 'sem'});
    lt_plot(x+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
    [h, p] = ttest2(y(OUTSTRUCT.All_traindir==-1), y(OUTSTRUCT.All_traindir==1));
    lt_plot_pvalue(p, '-1 vs 1', 1);
    lt_plot_zeroline;
    
    % --------------- SUMMARY
    lt_subplot(3,2,3); hold on;
    xlabel('train dir');
    ylabel('median (stim minus base)');
    y = OUTSTRUCT.All_FFmedianStimMinusNostim;
    
    % - run
    plot(OUTSTRUCT.All_traindir, y, 'ok');
    x = unique(OUTSTRUCT.All_traindir);
    [ymean, ysem] = grpstats(y, OUTSTRUCT.All_traindir, {'mean', 'sem'});
    lt_plot(x+0.1, ymean, {'Errors', ysem, 'Color', 'r'});
    [h, p] = ttest2(y(OUTSTRUCT.All_traindir==-1), y(OUTSTRUCT.All_traindir==1));
    lt_plot_pvalue(p, '-1 vs 1', 1);
    lt_plot_zeroline;
end
