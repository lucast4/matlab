function OUTSTRUCT = lt_opto_ExtrBirdDat(dirtoplot, twind, SwitchTimes, ...
    plotON, StartDaySkipTime, onlylongepoch)
% --- only keep if onset of epoch is at least 24 hours before onset of
% stim
if ~exist('onlylongepoch', 'var')
    onlylongepoch=0;
end
% onlylongepoch=0;
lightson = 7; % 7am

%% NOTE:
% - if multiple siwtches on  agiven day, takes summary of only data preceding
% first switch
% - this means will not include stim that occurs right after an epoch
% switch.
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
if plotON==1
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
            x = ceil(tvals(end))-0.1;
            lt_plot(x, mean(ffvals), {'Errors', lt_sem(ffvals), 'Color', pcol});
        end
        % ---------- SKIP EARLY DATA?
        % --- assume 7am wakeup time
        tmp = floor(min(tvals))+lightson/24 + StartDaySkipTime/24;
        line([tmp tmp], ylim, 'Color', [0.7 0.3 0.7]);
        
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



%% =================== PLOT RAW DATA [MEAN VALUES
if plotON==1
    lt_figure; hold on;
    allff = [];
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
            allff = [allff mean(ffvals)];
            % ======= PLOT
            % - raw
            pcol = pcolstruct.(fthis{1});
            %             plot(tvals, ffvals, 'x', 'Color', pcol);
            % - mean
            x = median(tvals);
            if strcmp(fieldsthis, 'All')
                
                lt_plot(x, mean(ffvals), {'Errors', lt_sem(ffvals), 'Color', pcol});
            end
            
            
        end
        
    end
    
    %%
    if ~isempty(SwitchTimes)
        % ========= CONVERT SWITCH TIMES
        [SwitchTimesMat, SwitchTimesMat_conting] = ...
            lt_tools_SwitchTimes(SwitchTimes, PARAMS.global.LastDate_num+1);
        
        % ====== put lines for swich times
        YLIMreal = ylim;
        YLIM = [0 max(allff)];
        for i=1:length(SwitchTimes)/2
            dtime = datenum(SwitchTimes{(i*2-1)}, 'ddmmmyyyy-HHMM');
            conting = SwitchTimes{(i*2)};
            line([dtime dtime], [YLIMreal(1) YLIM(2)]);
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
OUTSTRUCT.All_FFdaymean = [];
OUTSTRUCT.All_Daynum = [];

% ===== indivbidual trials
OUTSTRUCT.AllTrialsStimEpoch_ff = {};
OUTSTRUCT.AllTrialsStimEpoch_tval = {};
OUTSTRUCT.AllTrialsStimEpoch_tstim_relonset = {};
OUTSTRUCT.AllTrialsStimEpoch_isStim = {};

[SwitchTimesMat, SwitchTimesMat_conting] = ...
    lt_tools_SwitchTimes(SwitchTimes, PARAMS.global.LastDate_num+1);

assert(length(PARAMS.indiv) == length(DATSTRUCT.data), 'I assumed that these indices are "days"');

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
        
        
        if isfield(DATSTRUCT.data{i}.timewindow{twind}, 'All')
            f3 = DATSTRUCT.data{i}.timewindow{twind}.All.ffvals;
        else
            f3 = [];
        end
        f_all = [f1 f2 f3];
        
        % -------- check that distributions for tvals are similar
        if (0) % ignore this since should be obiuos by inspection of trial data whether this is the case
            assert(ranksum(t1, t2)>0.01, 'are time values different? then is problem as this assumes interleaves stim/nostim');
        end
        
        %% ======== EXTRACT THIGNS THAT REQUIRE GOING BACK TO EARLIER PROCESSED DATA
        if length(PARAMS.indiv{i}.Params)>1
            %  then ignore the one for "All", i.e. only take during stim
            %  epoch
            indParam = [];
            for ii=1:length(PARAMS.indiv{i}.Params)
                if length(PARAMS.indiv{i}.Params{ii}.FieldsToCheck)==2 ...
                        & all(ismember({'NotStim_StimCatch', 'StimNotCatch'}, PARAMS.indiv{i}.Params{ii}.FieldsToCheck))
                    % === then this is for stim epoch...
                    indParam = ii;
                elseif length(PARAMS.indiv{i}.Params{ii}.FieldsToCheck)==2 ...
                        & all(ismember({'StimCatch', 'StimNotCatch'}, PARAMS.indiv{i}.Params{ii}.FieldsToCheck))
                    indParam = ii;
                end
            end
        else
            indParam = 1;
        end
        
        
        % ========== Load data
        tmp = strfind(PARAMS.indiv{i}.Params{indParam}.savefolder , 'bluejay');
        if PARAMS.indiv{i}.Params{indParam}.savefolder(tmp+7)=='3'
            PARAMS.indiv{i}.Params{indParam}.savefolder(tmp+7)='2'; % resasign to bluejay 2...
        end
        try
            if PARAMS.indiv{i}.Params{indParam}.PLOT_TimeWindow_savedir(tmp+7)=='3'
                PARAMS.indiv{i}.Params{indParam}.PLOT_TimeWindow_savedir(tmp+7)='2'; % resasign to bluejay 2...
            end
        catch err
        end
        
        
        if exist([PARAMS.indiv{i}.Params{indParam}.savefolder '/StatsStruct.mat'], 'file')
            ststruct = load([PARAMS.indiv{i}.Params{indParam}.savefolder '/StatsStruct.mat']);
        else
            ststruct = load([PARAMS.indiv{i}.Params{indParam}.PLOT_TimeWindow_savedir '/StatsStruct.mat']);
        end
        fnamesthis = fieldnames(ststruct.StatsStruct);
        
        tval_all = []; % datenum
        tSinceStim_all = []; % time from syl onset to stim (pos means stim occured first)
        isStim_all = []; % 1 for stim, 0 for stim catch or no stim
        pitch_all = [];
        twindname = PARAMS.indiv{i}.Params{indParam}.TimeField{twind};
        
        for ff=1:length(fnamesthis)
            fthis = fnamesthis{ff};
            tval = ststruct.StatsStruct.(fthis).datenum;
            tsincetrig = ststruct.StatsStruct.(fthis).TimeSinceLastTrig; % THIS IS LASER STIM.
            pitch = ststruct.StatsStruct.(fthis).WINDOWED.(twindname).Pitch.vals';
            tval2 = ststruct.StatsStruct.(fthis).WINDOWED.(twindname).Time.datenum;
            %            tval2 = DATSTRUCT.data{i}.timewindow{twind}.(fthis).timevals_dnum;
            %                 outlierinds = ststruct.StatsStruct.(fthis).WINDOWED.(twindname).OutlierInds;
            
            % === sanity check to make sure both dartasets are aligned
            if ~(all(tval==tval2))
                
                disp('why arent they aligned? likely becuase of outliers:...');
                disp(outlierinds);
                keyboard
            end
            
            % ============= SAVE TO OUTPUT
            tval_all = [tval_all; tval'];
            pitch_all = [pitch_all; pitch];
            tSinceStim_all = [tSinceStim_all; tsincetrig'];
            
            if strcmp(fthis, 'StimNotCatch')
                % --- STIM ON
                isstim = ones(size(tval))';
            elseif strcmp(fthis, 'NotStim_StimCatch') | strcmp(fthis, 'StimCatch')
                % ---- STIM OFF or CATCH
                isstim = zeros(size(tval))';
            else
                disp('WHAT IS THIS?')
                keyboard;
            end
            isStim_all = [isStim_all; isstim];
        end
        
        %% ========================
        % NOTE: two codes. fgirst is old version. then second is updated to
        % only keep if is before the switch on a given day. if there is no
        % switch on a day then the output is identical (i.e. will not
        % subsample t and f). run both is fine, as the second does sanimity
        % check to make sure they are aligned.
        
        % --- what direction learning today?
        epochthis = find(all(([t1 t2]>SwitchTimesMat(:,1))') & ...
            all(([t1 t2]<SwitchTimesMat(:,2))')); % find the current "epoch"
        if length(epochthis)~=1
            % then these trials are not in just one well-defined learning epoch
            % learndir is not well defined
            learndir = nan;
        else
            learndir = SwitchTimesMat_conting(epochthis);
        end
        
        % =================== IF THERE IS SWITCH WITHIN TRIALS, THEN TAKE
        % ALL TRIALS B EFORE SWIETCH AND THROW OUT REST OF DATA
        %         if isnan(learndir) % then there is siwtch within trials ...
        assert(length(unique(floor([t1 t2])))==1, 'dat should all be from one dat');
        tmp = mean([t1 t2]>SwitchTimesMat(:,1),2); % any epoch with onset within trials will give non-(0 or 1) output
        %             assert(~all(ismember(tmp, [0 1])), 'no onset within trials ...?');
        
        % ----
        if ~isempty(epochthis)
            assert(epochthis==find(tmp==1, 1, 'last'), 'this makes sure this code ("=== IF THER ER IS...) is identical to above (what direction learnig.., qwhich it replaces');
        end
        
        % ----
        epochthis = find(tmp==1, 1, 'last');
        learndir = SwitchTimesMat_conting(epochthis);
        offsetthis = SwitchTimesMat(epochthis,2);
        
        % ------------ SUBSAMPLE TRIALS
        indstmp = t1<offsetthis;
        t1 = t1(indstmp);
        f1 = f1(indstmp);
        
        indstmp = t2<offsetthis;
        t2 = t2(indstmp);
        f2 = f2(indstmp);
        
        % ======================= if switch starts today, then skip
        onsetthis = SwitchTimesMat(epochthis,1);
        if floor(onsetthis)==unique(floor(t1))
            continue
        elseif floor(onsetthis)>unique(floor(t1))
            % not posible...
            fasdfsdafasfasdf;
        end
        
        
        if onlylongepoch==1
            t_firsttrial = min([t1 t2]);
            epochdur = t_firsttrial - onsetthis;
            if epochdur<1
                continue
            end
        end
        
        
        % ========================== throw out data at start of day
        tmp = floor(min(t1))+lightson/24 + StartDaySkipTime/24; % assume lgihts on at 7am
        
        indstmp = t1>tmp;
        t1 = t1(indstmp);
        f1 = f1(indstmp);
        
        indstmp = t2>tmp;
        t2 = t2(indstmp);
        f2 = f2(indstmp);
        
        
        %% ======================= MAP BACK TO FIGURE OUT TRIG TIMES
        idunique = pitch_all-tval_all; % each trial marked by unique combination of pitch and time.
        assert(length(unique(idunique))==length(idunique), 'the follwoing assumes each trial assigned to unqiue pitch...');
        %         assert(length(unique(pitch_all))==length(pitch_all), 'the follwoing assumes each trial assigned to unqiue pitch...');
        
        %         [tmp, ind1] = intersect(pitch_all, f1); assert(all(ismember(tmp, f1)), 'did not matchall of f, why?');
        %         tstim1 = tSinceStim_all(ind1);
        %         isstim1 = isStim_all(ind1);
        %
        %         [tmp, ind1] = intersect(pitch_all, f2); assert(all(ismember(tmp, f2)), 'did not matchall of f, why?');
        %         tstim2= tSinceStim_all(ind1);
        %         isstim2 = isStim_all(ind1);
        
        [tmp, ind1] = intersect(idunique, f1-t1); assert(all(ismember(tmp, f1-t1)), 'did not matchall of f, why?');
        tstim1 = tSinceStim_all(ind1);
        isstim1 = isStim_all(ind1);
        
        [tmp, ind1] = intersect(idunique, f2-t2); assert(all(ismember(tmp, f2-t2)), 'did not matchall of f, why?');
        tstim2= tSinceStim_all(ind1);
        isstim2 = isStim_all(ind1);
        
        % ==================== SAVE INDIVIDUAL TRIALS AND STIM TIMES
        ffall = [f1 f2];
        timeall = [t1 t2];
        tstimall = [tstim1' tstim2'];
        isstimall = [isstim1' isstim2'];
        
        if isempty(tstimall)
            % ========== not sure why... just ignore this day I guess
            tstimall = nan;
            isstimall = nan;
        end
        if length(isstimall)~=length(ffall)
            % NOT EXACTLY SURE WHY TO BE HONEST...
            tstimall = nan;
            isstimall = nan;
        end
        %% calculate stats
        dprime = (mean(f2)-mean(f1))/(sqrt(0.5*(var(f1)+var(f2))));
        ffmeandiff_stim = mean(f2)-mean(f1);
        ffmedian = median(f2) - median(f1);
        p = ranksum(f2, f1);
        ffdaymean = mean([f2 f1]);
        
        % --- overlay on plot
        if plotON==1
            % - mean
            %             x = ceil(t1(end))+0.1;
            x = ceil(t1(end))-0.2;
            % -- stim
            lt_plot(x, mean(f2), {'Errors', lt_sem(f2), 'Color', pcolstruct.StimNotCatch, 'MarkerFaceColor', 'none'});
            %                 plot(x, mean(f2), 'o', 'Color', pcolstruct.StimNotCatch);
            % -- not stim
            lt_plot(x, mean(f1), {'Errors', lt_sem(f1), 'Color', pcolstruct.StimCatch, 'MarkerFaceColor', 'none'});
            %                 plot(x, mean(f1), 'o', 'Color', pcolstruct.Stimatch);
        end
        
        % =================== OUTPUT
        OUTSTRUCT.All_Daynum = [OUTSTRUCT.All_Daynum; i];
        OUTSTRUCT.All_Dprime = [OUTSTRUCT.All_Dprime; dprime];
        OUTSTRUCT.All_FFdaymean = [OUTSTRUCT.All_FFdaymean; mean(f_all)];
        OUTSTRUCT.All_FFmeanStimMinusNostim = [OUTSTRUCT.All_FFmeanStimMinusNostim; ffmeandiff_stim];
        OUTSTRUCT.All_FFmedianStimMinusNostim = [OUTSTRUCT.All_FFmedianStimMinusNostim; ffmedian];
        OUTSTRUCT.All_pvalranksum = [OUTSTRUCT.All_pvalranksum; p];
        OUTSTRUCT.All_traindir = [OUTSTRUCT.All_traindir; learndir];
        
        OUTSTRUCT.AllTrialsStimEpoch_ff = [OUTSTRUCT.AllTrialsStimEpoch_ff; ffall];
        OUTSTRUCT.AllTrialsStimEpoch_tval = [OUTSTRUCT.AllTrialsStimEpoch_tval; timeall];
        OUTSTRUCT.AllTrialsStimEpoch_tstim_relonset = [OUTSTRUCT.AllTrialsStimEpoch_tstim_relonset; tstimall];
        if isempty(isstimall)
            keyboard
        end
        OUTSTRUCT.AllTrialsStimEpoch_isStim = [OUTSTRUCT.AllTrialsStimEpoch_isStim; isstimall];
%         if length(isstimall)~=length(ffall)
%             keyboard
%         end
    end
end


%% ========= REMOVE ANY DAYS THAT HAVE NAN FOR LEARN DIR

indstokeep = ~isnan(OUTSTRUCT.All_traindir);
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);



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
