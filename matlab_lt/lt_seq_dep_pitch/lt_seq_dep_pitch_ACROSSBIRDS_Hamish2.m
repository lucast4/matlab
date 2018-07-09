%% lt 5/24/18 - does direction of AFP bias change during learning?
% ACROSS EXPERIMENTS FOR THE SAME BIRD/SYLLABLE?
% NOTE: sumary data for PBS and MUSC (summary across expts) uses mean of
% day means, with perfectly matched days across PBS and MUSC.

%% LT 6/20/16
function lt_seq_dep_pitch_ACROSSBIRDS_Hamish2(SeqDepPitch_AcrossBirds, Params, ...
    plotsametype)

%% 1) SORT OUT ONLY THE THE EXPEIRMENTS THAT HAVE LMAN INACTIVATION DATA
% copy strcuture, save backup.
filter = 'LMAN';
[SeqDepPitch_AcrossBirds, NumBirds]=lt_seq_dep_pitch_ACROSSBIRDS_ExtractStruct(SeqDepPitch_AcrossBirds, filter);

%%
useLearningRelLastBlineDay=0; % otherwise will use rel entire baseline
useWNday2=0; % otherwise will use 1. using 2 allows to get some expt without songs on day 1
takeDay2forExptWithNoDay1=1; % otherwise will throw out those experiments.


%% ========================= PLOT EACH EXPERIMENTS

       All_learndir = [];
       All_basebiaspval = [];
       All_FF_BASE_PBS = [];
       All_FF_WN_PBS = [];
       All_FF_BASE_MUSC = [];
       All_FF_WN_MUSC = [];
        All_Birdnum = [];
        
        All_WNdayrange = [];
        
       All_CV_BASE_PBS = [];
       All_CV_WN_PBS = [];
       All_CV_BASE_MUSC = [];
       All_CV_WN_MUSC = [];

       
count = 0;
for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    %     figcount=1;
    %     subplotrows=4;
    %     subplotcols=2;
    %     fignums_alreadyused=[];
    %     hfigs=[];

    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        sametypesyls = [SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STDS, ...
            SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique_STSS];
        
        if plotsametype==0
            SylsToPlot={targsyl};
        elseif plotsametype==1
            SylsToPlot = [{targsyl} sametypesyls];
        end
        
        % ====== PLOT RAW DAT FOR THIS DAY TO COMPARE TO EXTRACTED STATS
        plotLarge=1;
        BirdToPlot=birdname;
        ExptToPlot=exptname;
        overlayMeans=1;
        UseSylColors=0; % 0 is black;
        flipsign=1; % plot pos, if neg
        use_std=0; % applies to mean plots only!! (std, instead of sem)
        plotRawFF=1; % 1=raw FF, 0= baseline subtracted (MUSC-MUSC)
        OverlayLMANStats=1; % plots mean shift in single target learning window (defined below)
        OverlayMUSC_days=[];
        plotRunningCV = 0;
        lt_seq_dep_pitch_ACROSSBIRDS_PlotRawLMAN(SeqDepPitch_AcrossBirds,...
            Params, BirdToPlot, ExptToPlot, SylsToPlot, overlayMeans, ...
            plotRawFF, UseSylColors, flipsign, use_std, OverlayLMANStats, ...
            OverlayMUSC_days, plotLarge, plotRunningCV)
       
        count = count+1;
        
       % ################################ COLLECT STATS
%        ffminusbase_PBS = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl).meanFF_pbs;
%        ffminusbase_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl).meanFF_musc;
%        
%        ffbase_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(targsyl).meanFF_WithinTimeWindow;
%        ffbase_PBS = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(targsyl).meanFF_WithinTimeWindow;
%        
%        
%        fflearn_MUSC = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl)
       
       
       % ############################## COLLECT STATS
       DatThis = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning;
       
       if isfield(DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC, 'final_extracted_window')
       
       %% ============================================= MUSC
       % ============================= BASE
       tvals = DatThis.AllDays_PlotLearning.EpochData_MUSC.Baseline.(targsyl).Tvals_WithinTimeWindow;
       ffvals = DatThis.AllDays_PlotLearning.EpochData_MUSC.Baseline.(targsyl).rawFF_WithinTimeWindow;
       
       [ffmean, ffstd] = grpstats(ffvals, floor(tvals), {'mean', 'std'}); % mean for each day
       ffmean_BASE_MUSC = mean(ffmean); % mean of means
       basedays = unique(floor(tvals));
       % --- save for comparison
       ffvals_base_MUSC = ffvals;
       
       % ############## CV
       ffCV_BASE_MUSC = mean(ffstd./ffmean);
       
       
       % ============================ DUR WN
       % days to extract, during learing.
       WNdays = DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.dayInds;
       
       % tvals, confirm that is same as previous extraction
       tvals_OLD = DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl).TvalsWithinWindow_MUSC;
       tvals = cell2mat(DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).Tvals_WithinTimeWindow(WNdays));
       assert(all(tvals == tvals_OLD), 'not aligned with previous data?');
       
       % mean of day means of FF
       ffmean = cellfun(@mean, DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).FFvals_WithinTimeWindow(WNdays));
       ffmean_WN_MUSC = mean(ffmean);
       
       ffstd = cellfun(@std, DatThis.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).FFvals_WithinTimeWindow(WNdays));
       ffCV_WN_MUSC = mean(ffstd./ffmean);
       
       
       %% ============================================= PBS
       % ============================ BASE
       tvals = DatThis.AllDays_PlotLearning.EpochData.Baseline.(targsyl).Tvals_WithinTimeWindow;
       ffvals = DatThis.AllDays_PlotLearning.EpochData.Baseline.(targsyl).rawFF_WithinTimeWindow;
       
       % -- only keep days that overlap with MUSC baseline days
       indstmp = ismember(floor(tvals), basedays);
       tvals = tvals(indstmp);
       ffvals = ffvals(indstmp);
       
       [ffmean, ffstd] = grpstats(ffvals, floor(tvals), {'mean', 'std'}); % mean for each day
       ffmean_BASE_PBS = mean(ffmean); % mean of means
       assert(all(unique(floor(tvals)) == basedays), 'PBS and MUSC not aligned...');
       % --- save for comparison
       ffvals_base_PBS = ffvals;
       
       % -------- CV
       ffCV_BASE_PBS = mean(ffstd./ffmean);

       

       % =========================== DUR WN
       % days to extract, during learing.
       
       % tvals, confirm that is same as previous extraction
       tvals_OLD = DatThis.AllDays_PlotLearning.EpochData.PBS_and_MUSC.final_extracted_window.(targsyl).TvalsWithinWindow_PBS;
       tvals = cell2mat(DatThis.AllDays_PlotLearning.DataMatrix.(targsyl).Tvals_WithinTimeWindow(WNdays));
       assert(all(tvals == tvals_OLD), 'not aligned with previous data?');
       
       % mean of day means of FF
       ffmean = cellfun(@mean, DatThis.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow(WNdays));
       ffmean_WN_PBS = mean(ffmean);
       
       % ------ CV
       ffstd = cellfun(@std, DatThis.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals_WithinTimeWindow(WNdays));
       ffCV_WN_PBS = mean(ffstd./ffmean);
       
       %% ====================== PLOT EXTRACTED VALUES ON LEARNING
       % TRAJECTORY
       firstday = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.SeqFilter.FirstDay;
       basedayinds = lt_convert_EventTimes_to_RelTimes(firstday, basedays);
       basedayinds = basedayinds.JustDays_rel;

       % --- PBS
       line([basedayinds(1)-0.5 basedayinds(end)+0.5], [ffmean_BASE_PBS ffmean_BASE_PBS], 'Color', 'b', 'LineWidth', 3);
       line([WNdays(1)-0.5 WNdays(end)+0.5], [ffmean_WN_PBS ffmean_WN_PBS], 'Color', 'b', 'LineWidth', 3);
       
       % --- MUSC
       line([basedayinds(1)-0.5 basedayinds(end)+0.5], [ffmean_BASE_MUSC ffmean_BASE_MUSC], 'Color', 'm', 'LineWidth', 3);
       line([WNdays(1)-0.5 WNdays(end)+0.5], [ffmean_WN_MUSC ffmean_WN_MUSC], 'Color', 'm', 'LineWidth', 3);
       
       
       % ############################# COLLECT STATS
       learndir = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
       All_learndir = [All_learndir; learndir];
       
       [~, p] = ttest2(ffvals_base_MUSC, ffvals_base_PBS);
       if p<0.05
           lt_plot_pvalue(p, 'base bias (ttest)', 1);
       end
       All_basebiaspval = [All_basebiaspval; p];
       
       All_FF_BASE_PBS = [All_FF_BASE_PBS; ffmean_BASE_PBS];
       All_FF_WN_PBS = [All_FF_WN_PBS; ffmean_WN_PBS];
       All_FF_BASE_MUSC = [All_FF_BASE_MUSC; ffmean_BASE_MUSC];
       All_FF_WN_MUSC = [All_FF_WN_MUSC; ffmean_WN_MUSC];
       
       All_CV_BASE_PBS = [All_CV_BASE_PBS; ffCV_BASE_PBS];
       All_CV_WN_PBS = [All_CV_WN_PBS; ffCV_WN_PBS];
       All_CV_BASE_MUSC = [All_CV_BASE_MUSC; ffCV_BASE_MUSC];
       All_CV_WN_MUSC = [All_CV_WN_MUSC; ffCV_WN_MUSC];
       
       
       All_Birdnum = [All_Birdnum; i];
       
       WNday1 = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd + ...
           SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumEmptyDays_StartWN_FromZscoreCode;
       All_WNdayrange = [All_WNdayrange; [WNdays(1) WNdays(end)]-WNday1+1];

       
       else
          lt_plot_annotation(1, 'NOT COLLECTING DATA! no dat in final window...', 'm'); 
       end
       
       
    end
end
disp(['Num expts plotted: ' num2str(count)]);


%% ============= PLOT SUMMARY
lt_figure; hold on;

% --------
lt_subplot(3,2,1); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('FF minus base dur learn (dir of learn)');

bias = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
learn_PBS = All_learndir.*(All_FF_WN_PBS - All_FF_BASE_PBS);
learn_MUSC = All_learndir.*(All_FF_WN_MUSC - All_FF_BASE_PBS);
for j=1:length(bias)
    line([bias(j) bias(j)], [learn_PBS(j) learn_MUSC(j)], 'Color', [0.7 0.7 0.7]);
end 
plot(bias, learn_PBS, 'ko');
plot(bias, learn_MUSC, 'ro');


% ---------
lt_subplot(3,2,2); hold on;
xlabel('learn (targ dir)');
ylabel('consolidated learning (hz)');

allbias = All_FF_BASE_PBS - All_FF_BASE_MUSC;
x = All_learndir.*(All_FF_WN_PBS - All_FF_BASE_PBS);
y = All_learndir.*(All_FF_WN_MUSC - All_FF_BASE_PBS);

% -- bias in direciton fo learing
indstokeep = sign(allbias) == sign(All_learndir)

xtmp = x(indstokeep);
ytmp = y(indstokeep);
plot(xtmp, ytmp, 'bo');


% -- bias in direciton opposite to learing
indstokeep = sign(allbias) ~= sign(All_learndir);

xtmp = x(indstokeep);
ytmp = y(indstokeep);
plot(xtmp, ytmp, 'mo');

% ---
lt_plot_makesquare_plot45line(gca, 'k');

% ---------
lt_subplot(3,2,3); hold on;
ylabel('fraction consolidation');
xlabel('afp bias (in dir of learnig');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = (All_FF_WN_MUSC - All_FF_BASE_PBS)./(All_FF_WN_PBS - All_FF_BASE_PBS);
lt_regress(y, x, 1);
% ---- connect same bird with lines
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
    
    [~, indsort] = sort(xthis);
    xthis = xthis(indsort);
    ythis = ythis(indsort);
    
    for k=1:length(xthis)-1
       line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]); 
    end
end
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ------
lt_subplot(3,2,4); hold on;
title('no difference in duration from learn day1');
xlabel('afp bias, dir of learn');
ylabel('days (b=1st; r=last)');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
yblue = All_WNdayrange(:,1);
yred = All_WNdayrange(:,2);

% lt_regress(yblue, x, 1, 0, 1, 1, 'b');
% lt_regress(yred, x, 1, 0, 1, 1, 'r');
plot(x, yred+0.1, 'or');
plot(x, yblue-0.1, 'ob');
for i=1:length(x)
line([x(i) x(i)], [yblue(i)-0.1 yred(i)+0.1], 'Color', [0.6 0.6 0.6]);
end


% ===================== CV CHANGE FROM BASELINE
lt_subplot(3,2,5); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('CV (bu=base; rd=Train)');
title('PBS data');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y_blue = All_CV_BASE_PBS;
y_mag = All_CV_WN_PBS;

plot(x, y_blue, 'ob');
plot(x, y_mag, 'or');
for i=1:length(x)
    if y_blue(i)>y_mag(i)
        line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'b'); 
    else
       line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'r');
    end
end
lt_plot_zeroline;

% -----------------
lt_subplot(3,2,6); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('CV (bu=base; rd=Train)');
title('MUSC data');

x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y_blue = All_CV_BASE_MUSC;
y_mag = All_CV_WN_MUSC;

plot(x, y_blue, 'ob');
plot(x, y_mag, 'or');
for i=1:length(x)
    if y_blue(i)>y_mag(i)
        line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'b'); 
    else
       line([x(i) x(i)], [y_blue(i) y_mag(i)], 'Color', 'r');
    end
end
lt_plot_zeroline;



%% =========================== CV THINGS
lt_figure; hold on;

% ============== 1) 
lt_subplot(3,2,1); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('change in CV (difference)');
title('PBS');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = All_CV_WN_PBS - All_CV_BASE_PBS;
lt_regress(y, x, 1);
lt_plot_zeroline;

% ============== 
lt_subplot(3,2,2); hold on;
xlabel('AGAINST -- TOWARDS');
ylabel('Change in CV');
title('PBS');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
x = x>0;
y = All_CV_WN_PBS - All_CV_BASE_PBS;

% plot(x, y, 'ok');

[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});

lt_plot([0 1]+0.2, ymean, {'Errors', ysem, 'Color', 'r'});

% --- lines between paired experiments
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
%     
%     [~, indsort] = sort(xthis);
%     xthis = xthis(indsort);
%     ythis = ythis(indsort);
%     
%     for k=1:length(xthis)-1
%        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]); 
%     end
    % plot with own jitter and color
    xthis = xthis + 0.3*(rand-0.5);
    pcol = 0.9*[rand rand rand];
    lt_plot(xthis, ythis, {'Color', pcol});    
end


xlim([-1 2]);
lt_plot_zeroline;
p = ranksum(y(x==0), y(x==1));
lt_plot_pvalue(p, 'ranksum',1);


% ============== 
lt_subplot(3,2,3); hold on;
xlabel('AGAINST -- TOWARDS');
ylabel('Change in CV');
title('MUSC');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
x = x>0;
y = All_CV_WN_MUSC - All_CV_BASE_MUSC;

% plot(x, y, 'ok');

[ymean, ysem] = grpstats(y, x, {'mean', 'sem'});

lt_plot([0 1]+0.2, ymean, {'Errors', ysem, 'Color', 'r'});

% --- lines between paired experiments
for j=1:NumBirds
    indsthis = find(All_Birdnum==j);
    
    xthis = x(indsthis);
    ythis = y(indsthis);
%     
%     [~, indsort] = sort(xthis);
%     xthis = xthis(indsort);
%     ythis = ythis(indsort);
%     
%     for k=1:length(xthis)-1
%        line([xthis(k) xthis(k+1)], [ythis(k) ythis(k+1)], 'Color', [0.5 0.5 0.5]); 
%     end
    % plot with own jitter and color
    xthis = xthis + 0.3*(rand-0.5);
    pcol = 0.9*[rand rand rand];
    lt_plot(xthis, ythis, {'Color', pcol});    
end


xlim([-1 2]);
lt_plot_zeroline;
p = ranksum(y(x==0), y(x==1));
lt_plot_pvalue(p, 'ranksum',1);




% ============== 1) 
lt_subplot(3,2,4); hold on;
xlabel('baseline AFP bias (pos = dir of learning)');
ylabel('change in CV (difference)');
title('MUSC');
x = All_learndir.*(All_FF_BASE_PBS - All_FF_BASE_MUSC);
y = All_CV_WN_MUSC - All_CV_BASE_MUSC;
lt_regress(y, x, 1);
lt_plot_zeroline;


% ============== 1) 
lt_subplot(3,2,5); hold on;
xlabel('magnitude of AFP bias (during train)');
ylabel('change in CV (train - base)');
title('PBS');
x = All_learndir.*(All_FF_WN_PBS - All_FF_WN_MUSC);
y = All_CV_WN_PBS - All_CV_BASE_PBS;
lt_regress(y, x, 1);
lt_plot_zeroline;

% ===============
lt_subplot(3,2,6); hold on;
xlabel('magnitude of AFP bias, baseline normalized (during train)');
ylabel('change in CV (train - base)');
title('PBS');
x =  All_learndir.* ((All_FF_WN_PBS - All_FF_BASE_PBS) - (All_FF_WN_MUSC - All_FF_BASE_MUSC));
y = All_CV_WN_PBS - All_CV_BASE_PBS;
lt_regress(y, x, 1);
lt_plot_zeroline;


%% ========== BASELINE BIAS A FUNCTION OF BASELINE PITCH (PBS)

lt_figure; hold on;

xlabel('baseline pitch (PBS)');
ylabel('baseline pitch (MUSC)');

x = All_FF_BASE_PBS;
y = All_FF_BASE_MUSC;

plot(x,y, 'ok');
lt_plot_makesquare_plot45line(gca, 'b');















