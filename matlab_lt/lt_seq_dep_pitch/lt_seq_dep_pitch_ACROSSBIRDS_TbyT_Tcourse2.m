function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2(TrialStruct, ParamsTrial, ...
    ignoreLMANexpt, songasrend)

%% 7/31/18 - lt diverged from lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse, here 
% focusing on key extractions and analyses. 


%% ##################################################################
%% ############################# PREPROCESSING

Numbirds = length(TrialStruct.birds);

%% ======= CONVERT TO USING SONG AS REND?
if songasrend==1
    disp('CONVERTING FROM USING RENDS TO USING SONGS!');
    for i=1:Numbirds
        Numexpt = length(TrialStruct.birds(i).exptnum);
        
        for ii=1:Numexpt
            
            % ---------- SKIP IF NO DATA
            if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
                disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                continue
            end
            
            Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
            
            % ================== go thru all syls
            for ss =1:Numsyls
                
                t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
                t_dnum = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals_datenum;
                ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
                
                % ----------------- GO THRU EACH SONG AND COLLECT
                ff_songs = grpstats(ff, t, {'mean'});
                
                [~, indtmp] = unique(t);
                t_songs = t(indtmp);
                t_dnum = t_dnum(indtmp);
                
                assert(length(t_songs) == length(ff_songs));
                
                
                % ----------------- PUT BACK INTO STRUCTURE                
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals = ff_songs;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals = t_songs;
                TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals_datenum = t_dnum;
                
                % ------------- for all renditions calculate deviation from
                % recent trials
                % === method1 - fit regression line to one hour of data
                % (directly preceding this rendition...) record deviation from
                % that hour's prediction
                
            end
        end
    end
end

%% ====== CONVERT ALL FF TO FF DEVIATION FROM BASELINE (IN LEARN DIR)
% AND DECIDE IF USE SONG AS REND

disp('NOTE TO SELF: ffvals will be automatically flipped to be dir of learning (so no need to flip again later)');
pause;

for i=1:Numbirds
    Numexpt = length(TrialStruct.birds(i).exptnum);
    
    for ii=1:Numexpt
        
        % ---------- SKIP IF NO DATA
        if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
            continue
        end
        
        % --------- ignore if lMAN?
        if ignoreLMANexpt==1
            if isfield(TrialStruct.birds(i).exptnum(ii), 'LMANinactivated')
                % if no field then is not LMAN inact experiemts...
                isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
                if isLMAN==1
                    disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                    continue
                end
            end
        end
        
        Numsyls = length(TrialStruct.birds(i).exptnum(ii).sylnum);
        %         birdname =  TrialStruct.birds(i).birdname;
        %         exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        targlearndir = TrialStruct.birds(i).exptnum(ii).targlearndir;
        
        % =========== collect syls for this experiment
        %         ffedges_allsyls =[];
        %         tedges_allsyls = [];
        %         istarg_allsyls = [];
        %         issame_allsyls =[];
        %         sylnames_allsyls = {};
        
        % ===================== FOR ALL SYLS, each rendition collect time and
        % deviation from running avg.
        
        for ss =1:Numsyls
            
            
            
            % ============== subplot for this syl
            t = TrialStruct.birds(i).exptnum(ii).sylnum(ss).Tvals;
            ff = TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals;
            
            % --------------- lines for base mean and 1std
            WNontime = TrialStruct.birds(i).exptnum(ii).WNontime;
            
            %             basedays = TrialStruct.birds(i).exptnum(ii).BaseDays;
            %             wndays = TrialStruct.birds(i).exptnum(ii).WNDays;
            indsbase = t<WNontime;
            %             indsbase = t<basedays(end)+1;
            ffmean_base = nanmean(ff(indsbase));

            
            % --------------- subtract mean and flip if negative larning
            ff = (ff-ffmean_base).*targlearndir;
            TrialStruct.birds(i).exptnum(ii).sylnum(ss).FFvals = ff;
            
            % ------------- for all renditions calculate deviation from
            % recent trials
            % === method1 - fit regression line to one hour of data
            % (directly preceding this rendition...) record deviation from
            % that hour's prediction
                        
        end
    end
end

%% ========== FOR ALL RENDS, COMPUTE FLANKING RENDITIONS (FF AND T DEV)

% =============== BINS
% xedges = [-60 -31:5:-1 1:5:31 60]; % minutes
% xedges = [-30:10:-10 -0.1 0.1 10:10:30];
xedges = [-30 -0.5 0.5 30]; % GOOD, for comparison to baseline
% xedges = [-30:2.5:-2.5 -0.1 0.1 2.5:2.5:30];
% xedges = [-60 -30:5:30 60];
% xedges = [-60 -30:5:-5 -0.1 0.1 5:5:30 60];
% xedges = [-40 -20:5:20 40];
% xedges = [-60 -15 -3 0 3 15 60];
% xedges = [-1440 0 1440];
% % NOTE: the edge of xedges will be used as flank. so will only keep rends
% with at least that much time flanking.

xcenters = xedges(1:end-1)+diff(xedges)/2;

flipTargNontarg=0; % if 1, then for each targ/nontarg pair, flips the dataset
% so that is target dev relative to target pitch, as function of nontarget
% density [default =0];

useMeanFFRef = 0; % 1: uses mean in window; 0; uses value for each rend
flanktime_targ = 2.5; % minutes (will do +/- this number). This only used if
% note: this is used to decide how many rends of targ/nontarg are present
% during reference period.

tflankplot = min(abs([xedges(1) xedges(end)]))./60;
tflankplot = 0.5;
% tflankplot = 0;

twind_plot = [xedges(1) xedges(end)]./60;

TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
    ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg, ...
    useMeanFFRef);




%% ====== TO COLLECT DATA TO PLOT (NEW VERSION, COLLECTS FIRST, THEN PLOTS)
% Instead of collecitng and plotting at same time

% ======== what data to take?
% dotrain = 1; % during training?
singleRendOnly=1; % 1: only takes one datapt (first nontarg after the referen)
% if 0, then takes all rends in window


% ====== what data to use to compute density?
densitymethod = 'refperiod';
% refperiod: will use reference period
% entirebin: will use time from rendtion to maxtime_fromref (see maxtime_fromref below)
% beforefirstnontarg = will get number that is >= ref time(i.e.0) and < time of first nontarget post-reference.


% ============ method for decideing hi and lo density trials
cutoffmethod = 'medianslice';
% medianslice: for each nontarg bin, finds median for targ
% medianoverall: overall median of targ
% medianslice_rand: gets median by first addaing random jitter. useful for
% i) when even num bins. if only one y bin, then forces that to be in low
% density bin. if odd num bins, then will split favoring putting
% more data into low density category. 
% targffdev: then splits based on mean FF dev of target syllables 


% ======== summary plot: what time period to plot (locked to ref period)
% NOTE: This is also used to define period for estimating density if the
% choice above for densitymethod is entirebin.
% NOTE ONLY INPORTANT FOR DENSITY MEASURE.
mintime_fromref = 5;
maxtime_fromref = 30; % minutes


% =============== DEFAULTS
% singleRendOnly=0; % only take data up to first nontarg after the referen
% densitymethod = 'refperiod';
% cutoffmethod = 'medianslice';
% mintime_fromref = 5;
% maxtime_fromref = 30; % minutes
% minrends_inbin = 20;


% ================ RUN
DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt);


assert(all(isnan(DATBYREND.IsDurTrain) == isnan(DATBYREND.Density_isHigh)), 'then nan likely matches up');
assert(all(isnan(DATBYREND.Density_targ) == isnan(DATBYREND.Density_nontarg)), 'then nan likely matches up');



%% ################################################
%% ######################## PLOTS

%% TO DO:
if (0)
%% =========== [PLOT DIAGNOSTIC] FOR INDIVIDUAL EXPERIMENTS

if (1)
    birdtoplot = 13;
    expttotplot = 1;
    
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S7(DATBYREND, TrialStruct, ...
        birdtoplot, expttotplot, mintime_fromref, maxtime_fromref, xedges);
end

% ========== PLOT MULTIPLE
birdtoplot_list = [13:17];
% birdtoplot_list = 1:3:17;
expttoplot_list = [1:10];
for bb=birdtoplot_list
   for ee=expttoplot_list
       lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S7(DATBYREND, TrialStruct, ...
        bb, ee, mintime_fromref, maxtime_fromref, xedges);
   end
end
end



%% ================== SUMMARY - for each experiment plot all renditions, or subset
%% ====================== extract binned data
dosametype = 1; % 1=same; 0: diff type;
onlyifSigLearn = 0; % 1=yes, 0=don't care

% ----------------- % method for equalizing density
densitymethod = 'default'; % then whatever was used to classify has high and low density
% FFwindow; % then from 0(inclusive) to whatever window used to get FFdev

% ----------- min num trials in this bin...
minrends_inbin = 1;
minsongs_inbin = 1;

% ============================= TRAINING
dotrain = 1; % 1 = train, 0 = base
[~, ~, FFbinnedCell, TbinnedCell, ~, ...
    ~, NumDatPerRend, BirdNumMat] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_noDens(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin);


% ============================= BASELINE
dotrain = 0;
[~, ~, FFbinnedCell_BASE, TbinnedCell_BASE, ~, ...
    ~, NumDatPerRend_BASE, BirdNumMat_BASE] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_noDens(DATBYREND, ...
    dotrain, dosametype, onlyifSigLearn, densitymethod, xedges, mintime_fromref, ...
    maxtime_fromref, minrends_inbin, minsongs_inbin);


%% =================== distrubution of labeling densityl;
lt_figure; hold on;

% === distribution
lt_subplot(2,2,1); hold on;
indstmp = ~isnan(NumDatPerRend);
lt_plot_histogram(NumDatPerRend(indstmp));


lt_subplot(2,2,2); hold on;
xlabel('birdnum');
ylabel('renditions per locked syl');

indstmp = ~isnan(NumDatPerRend);
x = BirdNumMat(indstmp);
y = NumDatPerRend(indstmp);
plot(x,y, 'ok');



%% ===================== [PLOT] WHICH INDS TO PLOT?
onlyPlotIfDenseLabeling = 0; % decides if dense by looking at number of datapoints
% per rendition.
labelthresh = 0; % in log2 units, only matters if onlyPlotIfDenseLabeling=1
birdstokeep = [13:17]; % leave empty to take all birds

% ================== 1) ONLY PLOT IF HAVE DATA FOR BOTH LO AND HIGH DENS
% -------- 1) only plot if have data
indstoplot = ~cellfun('isempty', FFbinnedCell);

% -------- 2) only keep if dense labeling?
if onlyPlotIfDenseLabeling==1
     indstoplot = indstoplot & log2(NumDatPerRend)>labelthresh;
end

% --------- 3) birds to keep
if ~isempty(birdstokeep)
    indstoplot = indstoplot & ismember(BirdNumMat, birdstokeep);
end


% --------- FINAL, convert to numerical indeices
indstoplot = find(indstoplot)';



%% =================== [PLOT] SHOW SAMPLE SIZE

lt_figure; hold on;

% ############################ BINNED FF DEV 
lt_subplot(3,2,1); hold on;

for i=indstoplot

    x = TbinnedCell{i};
    y = FFbinnedCell{i};
    plot(x,y, '-x', 'Color', 'b');
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);

% ############################ BINNED FF DEV (BASELINE)
lt_subplot(3,2,2); hold on;

for i=indstoplot

    x = TbinnedCell_BASE{i};
    y = FFbinnedCell_BASE{i};
    plot(x,y, '-x', 'Color', 'b');
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);


% ############################## IN EACH BIN, COMPARE DAT VS. BASE
lt_subplot(3,2,3); hold on;

for i=indstoplot

    x = TbinnedCell{i};
    y = FFbinnedCell{i};
    
    xbase = TbinnedCell_BASE{i};
    ybase = FFbinnedCell_BASE{i};
    
%     X = [xbase'-2 x'+2];
%     Y = [ybase y];

    for cc = xcenters
       
        X = [cc-2 cc+2];
        Y = [ybase(xbase==cc) y(x==cc)];
        
        if length(Y)==2
            plot(X, Y, '-');
        end
        
    end
%     for ii=1:size(X,1)
%         plot(X(ii,:), Y(ii,:), '-');
%     end
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);



% ############################ BINNED FF DEV [mean across syls] 
lt_subplot(3,2,4); hold on;
title('mean');
X = [TbinnedCell{indstoplot}];
Y = cell2mat(FFbinnedCell(indstoplot))';

[ymean, ysem] = grpstats(Y, X, {'mean', 'sem'});
xmean = unique(X);
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', 'k'});


% % ###################################### TAKE ALL RENDS IN ONE BIN
% lt_subplot(3,2,6); hold on;
% 
% Y = FFsinglebinMat(indstoplot,:);
% X = [1 2];
% plot(X, Y', '-', 'Color', [0.7 0.7 0.7]);
% % -- means
% lt_plot(X+0.1, mean(Y), {'Errors', lt_sem(Y)});
% 
% xlim([0 3]);
% lt_plot_zeroline;
% % --- signifnicace
% [~, p] = ttest(Y(:,1), Y(:,2));
% lt_plot_pvalue(p, 'ttest', 1);

% ######################################## TRAINING AND BASELINE IN SAME
% PLOT
% -------- 1) only keep syls with both train and baseline
lt_subplot(3,2,5); hold on;
title('train vs base');
for i=indstoplot

    x = TbinnedCell{i};
    y = FFbinnedCell{i};
    
    xbase = TbinnedCell_BASE{i};
    ybase = FFbinnedCell_BASE{i};
    
    % only keep bins that both base and WN have
    [xboth, ind1, ind2] = intersect(x, xbase);
    
    y = y(ind1);
    ybase = ybase(ind2);
    
    if ~isempty(xboth)
    plot(xboth, y, 'r');
    plot(xboth, ybase, 'k-');
    end
end
lt_plot_zeroline;
lt_plot_zeroline_vert;
xlim([-30 30]);





%% ================= [SUMAMRY PLOTS] - GOOD
% -------------- INPUTPARAMS
getsame = 1;
gettrain = 1;
getsiglearn = 0;
birdstoget = 13:17;
nbin = 10; % for smothing;

% ===================================== 1) RAW DATA FOR ALL SYLS
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


sylmax = max(DATBYREND.Sylcounter);
% --- collect all times and plot distribution
t_all = [];
for ss = 1:sylmax
   
    if getsiglearn==1
        sigthis = 1;
    else
        sigthis = [0 1];
    end
       
    indsthis = DATBYREND.IsSame==getsame & DATBYREND.IsDurTrain==gettrain ...
        & DATBYREND.IsTarg==0 & ismember(DATBYREND.SigLearn, sigthis) ...
        & DATBYREND.Sylcounter==ss;
    
    if ~any(indsthis)
        continue
    end
   
    % ---- want this syl?
    bnum = unique(DATBYREND.Birdnum(indsthis));
    enum = unique(DATBYREND.Exptnum(indsthis));
    snum = unique(DATBYREND.Sylnum(indsthis));
    siglearn = unique(DATBYREND.SigLearn(indsthis));
    
%     if bnum==13 & enum==1 & snum==3
%         keyboard
%     end
%     
    if ~isempty(birdstoget)
        if ismember(bnum, birdstoget)==0
            continue
        end
    end
    
    
    % ============================ PLOT RAW DAT, AND SMOOTHED RUNNING
    tthis = cell2mat(DATBYREND.Time_dev(indsthis));
    tthis = tthis*(24*60);
    
    ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
    
    
    % =================== PLOT
    birdname = TrialStruct.birds(bnum).birdname;
    exptname = TrialStruct.birds(bnum).exptnum(enum).exptname;
    sylthis = TrialStruct.birds(bnum).exptnum(enum).sylnum(snum).syl;
    if siglearn==1
        pcol = 'b';
    elseif siglearn==0
        pcol = 'r';
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([birdname '-' exptname '-' sylthis]);
    
    plot(tthis, ffthis, 'x', 'Color',pcol);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % --------------- COLLECT FOR PLOT
    t_all = [t_all; tthis];
end

linkaxes(hsplots, 'xy');

% ======================= PLOT HISTO OF ALL TIMEDEVS
lt_figure; hold on;
title('all devs, all rends');
lt_plot_histogram(t_all);



% ################################ 1) SMOOTHED
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


sylmax = max(DATBYREND.Sylcounter);
% --- collect all times and plot distribution
t_all = [];
for ss = 1:sylmax
   
    if getsiglearn==1
        sigthis = 1;
    else
        sigthis = [0 1];
    end
       
    indsthis = DATBYREND.IsSame==getsame & DATBYREND.IsDurTrain==gettrain ...
        & DATBYREND.IsTarg==0 & ismember(DATBYREND.SigLearn, sigthis) ...
        & DATBYREND.Sylcounter==ss;
    
    if ~any(indsthis)
        continue
    end
   
    % ---- want this syl?
    bnum = unique(DATBYREND.Birdnum(indsthis));
    enum = unique(DATBYREND.Exptnum(indsthis));
    snum = unique(DATBYREND.Sylnum(indsthis));
    siglearn = unique(DATBYREND.SigLearn(indsthis));
    
%     if bnum==13 & enum==1 & snum==3
%         keyboard
%     end
%     
    if ~isempty(birdstoget)
        if ismember(bnum, birdstoget)==0
            continue
        end
    end
    
    
    % ============================ PLOT RAW DAT, AND SMOOTHED RUNNING
    tthis = cell2mat(DATBYREND.Time_dev(indsthis));
    tthis = tthis*(24*60);
    
    ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
    
    % ================ running smoother
    % - first sort
    [~, indsort] = sort(tthis);
    tthis = tthis(indsort);
    ffthis = ffthis(indsort);
    if length(tthis)<1.5*nbin
        continue
    end
    tthis = lt_running_stats(tthis, nbin);
    ffthis = lt_running_stats(ffthis, nbin);
    
    % =================== PLOT
    birdname = TrialStruct.birds(bnum).birdname;
    exptname = TrialStruct.birds(bnum).exptnum(enum).exptname;
    sylthis = TrialStruct.birds(bnum).exptnum(enum).sylnum(snum).syl;
    if siglearn==1
        pcol = 'b';
    elseif siglearn==0
        pcol = 'r';
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([birdname '-' exptname '-' sylthis]);
    
%     shadedErrorBar(tthis.Mean, ffthis.Mean, ffthis.SEM, {'Color', pcol});
    lt_plot(tthis.Median, ffthis.Mean, {'Errors', ffthis.SEM, 'Color', pcol});
    
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % --------------- COLLECT FOR PLOT
    t_all = [t_all; tthis];
end

linkaxes(hsplots, 'xy');


% ===================================== 1) DISTRIBUTION OF DEVIATIONS
% SANITY CHECK, SHOULD SUM TO THE TOTAL LEARNING...

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


sylmax = max(DATBYREND.Sylcounter);
LearnExpectAll = [];
LearnActualAll = [];
for ss = 1:sylmax
   
    if getsiglearn==1
        sigthis = 1;
    else
        sigthis = [0 1];
    end
       
    indsthis = DATBYREND.IsSame==getsame & DATBYREND.IsDurTrain==gettrain ...
        & DATBYREND.IsTarg==0 & ismember(DATBYREND.SigLearn, sigthis) ...
        & DATBYREND.Sylcounter==ss;
    
    if ~any(indsthis)
        continue
    end
   
    % ---- want this syl?
    bnum = unique(DATBYREND.Birdnum(indsthis));
    enum = unique(DATBYREND.Exptnum(indsthis));
    snum = unique(DATBYREND.Sylnum(indsthis));
    siglearn = unique(DATBYREND.SigLearn(indsthis));
    learnthis = unique(DATBYREND.LearnMag_regr(indsthis));
    
%     if bnum==13 & enum==1 & snum==3
%         keyboard
%     end
%     
    if ~isempty(birdstoget)
        if ismember(bnum, birdstoget)==0
            continue
        end
    end
    
    
    % ============================ PLOT RAW DAT, AND SMOOTHED RUNNING
    tthis = cell2mat(DATBYREND.Time_dev(indsthis));
    tthis = tthis*(24*60);
    
    ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
    
    
    % =================== PLOT
    birdname = TrialStruct.birds(bnum).birdname;
    exptname = TrialStruct.birds(bnum).exptnum(enum).exptname;
    sylthis = TrialStruct.birds(bnum).exptnum(enum).sylnum(snum).syl;
    if siglearn==1
        pcol = 'b';
    elseif siglearn==0
        pcol = 'r';
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([birdname '-' exptname '-' sylthis]);
    
    lt_plot_histogram(ffthis);
    lt_plot_zeroline_vert;
    
    % ================================== HOW MUCH TOTAL LEARNING FOR THIS
    % SYL?
    learn_expected = sum(ffthis);
    
    LearnExpectAll = [LearnExpectAll; learn_expected];
    LearnActualAll = [LearnActualAll; learnthis];
end

    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['all syls']);
    xlabel('expected learnign (sum of all ff devs)');
    ylabel('actual learning (eends of regression line)');
    plot(LearnExpectAll, LearnActualAll, 'ok');
    lt_plot_makesquare_plot45line(gca, 'r');
    



    %% ==================== [SUMMARY PLOT, ADAPTIVE BINS]
    % I.E. for each syl, get divide into first and second half, matching sample
    % sizes
    % 3 bins: 0-4min, the rest divided into 2
    
    plotraw = 1;
    edgelist = [2 3 4 5]; % list of edges to use 
    edgelist = [3]; % list of edges to use 
    minrendsinbin = 4;
    
    getsame = 1;
    gettrain = 1;
    getsiglearn = 0;
    birdstoget = 13:17;
    colorscheme = 'learnsig';
    % choices: bird; learnsig
    
    %     plotraw = 0;
%     edgelist = [2 3 4 5]; % list of edges to use 

    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_1(DATBYREND, TrialStruct, getsame, gettrain, ...
    getsiglearn, birdstoget, plotraw, edgelist, minrendsinbin, colorscheme)
       
        
    %% ==================== [SUMMARY PLOT, INPUT NEW BINEDGES]
    % I.E. for each syl, get divide into first and second half, matching sample
    % sizes
    % 3 bins: 0-4min, the rest divided into 2
    
    % =====================
    xedgethis = [0:4:32];
    
    figcount=1;
    subplotrows=5;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    sylmax = max(DATBYREND.Sylcounter);
    
    % --- collect
    Xall = [];
    Yall = [];
    
    for ss = 1:sylmax
        
        if getsiglearn==1
            sigthis = 1;
        else
            sigthis = [0 1];
        end
        
        indsthis = DATBYREND.IsSame==getsame & DATBYREND.IsDurTrain==gettrain ...
            & DATBYREND.IsTarg==0 & ismember(DATBYREND.SigLearn, sigthis) ...
            & DATBYREND.Sylcounter==ss;
        
        if ~any(indsthis)
            continue
        end
        
        % ---- want this syl?
        bnum = unique(DATBYREND.Birdnum(indsthis));
        enum = unique(DATBYREND.Exptnum(indsthis));
        snum = unique(DATBYREND.Sylnum(indsthis));
        siglearn = unique(DATBYREND.SigLearn(indsthis));
        
        
        if ~isempty(birdstoget)
            if ismember(bnum, birdstoget)==0
                continue
            end
        end
        
        
        % ============================ PLOT RAW DAT, AND SMOOTHED RUNNING
        tthis = cell2mat(DATBYREND.Time_dev(indsthis));
        tthis = tthis*(24*60);
        
        ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
        
        % --------- get 3 bins
        tbins = discretize(tthis, xedgethis);
        
        % ------ collect binned balues
        tbinned = grpstats(tthis, tbins, {'mean'});
        ffbinned = grpstats(ffthis, tbins, {'mean'});
        ffbinned_sem = grpstats(ffthis, tbins, {'sem'});
        
        
        % =================== PLOT
        birdname = TrialStruct.birds(bnum).birdname;
        exptname = TrialStruct.birds(bnum).exptnum(enum).exptname;
        sylthis = TrialStruct.birds(bnum).exptnum(enum).sylnum(snum).syl;
        if siglearn==1
            pcol = 'b';
        elseif siglearn==0
            pcol = 'r';
        end
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title([birdname '-' exptname '-' sylthis]);
        
        lt_plot(tbinned, ffbinned, {'Errors', ffbinned_sem, 'Color', pcol});
        
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
        
        
        % ---------- put lines for bin edges
        for j=1:length(xedgethis)
            line([xedgethis(j) xedgethis(j)], ylim, 'Color', 'r');
        end
        
        % --------- put individual datapoints
        plot(tthis, ffthis, '.', 'Color', pcol);
        
        % ------------- COLLECT
%         Xall = [Xall; tbinned'];
%         Yall = [Yall; ffbinned'];
        
    end
    
    linkaxes(hsplots, 'xy');
    
%     
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots hsplot];
%     title('summary, each syl');
%     
%     for j=1:size(Xall,1)
%         plot(Xall(j, :), Yall(j,:), '-o');
%     end
%     Ymean = mean(Yall,1);
%     Ysem = lt_sem(Yall);
%     Xmean = mean(Xall, 1);
%     lt_plot_bar(Xmean, Ymean, {'Errors', Ysem, 'Color', 'k'});
%     
% 
% 

























