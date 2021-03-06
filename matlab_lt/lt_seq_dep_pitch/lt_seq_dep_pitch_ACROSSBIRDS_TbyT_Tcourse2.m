function [DATBYREND] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2(TrialStruct, ParamsTrial, ...
    ignoreLMANexpt, songasrend, singleRendOnly, timewindow_refsyl)
%% lt 11/17/19 - to restrict analyses to certain time windows within day.

if ~exist('timewindow_refsyl', 'var')
    timewindow_refsyl = [];
    % should be hour to hour, will only keep trials for which the reference
    % syllable time falls within this window. e.g. [7 9] is 7 to 9 am.
end

%% NOTE (7/17/19) - removed TrialStruct from output since here ffvals are flipped. this may mess up subsequent analses
%% 7/31/18 - lt diverged from lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse, here
% focusing on key extractions and analyses.


%% ##################################################################
%% ############################# PREPROCESSING

Numbirds = length(TrialStruct.birds);

%% ======= CONVERT TO USING SONG AS REND?

if songasrend==1
    TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_SongRend(TrialStruct);
end

%% ====== CONVERT ALL FF TO FF DEVIATION FROM BASELINE (IN LEARN DIR)
% AND DECIDE IF USE SONG AS REND

if ~isfield(TrialStruct, 'FFalreadyFlippedLearnDir')
    TrialStruct.FFalreadyFlippedLearnDir = 0;
end

if TrialStruct.FFalreadyFlippedLearnDir==0
    disp('NOTE TO SELF: ffvals will be automatically flipped to be dir of learning (so no need to flip again later)');
    disp('IMPORTANT - this will permanently affect TrialStruct!')
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
                
                assert(abs(targlearndir)==1);
                
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
    TrialStruct.FFalreadyFlippedLearnDir = 1;
    
end

if (0)
    %% removed below since did not seem tobe using this anywhere
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
% singleRendOnly=1; % 1: only takes one datapt (first nontarg after the referen)
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

% ================== include target syl?
includeTarg=0;

% ================ RUN
DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt, includeTarg);


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
end

%% ###########################################################
%% ###########################################################
%% [IMPORTANT]
%% COLLECT ALL DATA. ANALYSIS WITHOUT SEPARATING INTO HI AND LO DENS
% [ALTHOUGH BECAUSE GOT ALL DATA, CAN EASILY PARSE BY THAT]
%% ==================== [COLLECT - 1] collect deviations

% =============== BINS
% xedges = [-65 -45 -25 -5 5 25 45 65]; % minutes
xedges = [-80 -55 -30 -5 5 30 55 80]; % minutes
xedges = [-60 -31:5:-1 1:5:31 60]; % minutes
% xedges = [-50:15:-5 -0.5 0.5 5:15:50]; % minutes
xedges = [-120 -60 -30 -7 -0.5 0.5 7 30 60 120]; % minutes
xcenters = xedges(1:end-1)+diff(xedges)/2;
% xcenters = xedges(1:end-1)+binsize/2;

flipTargNontarg=0; % if 1, then for each targ/nontarg pair, flips the dataset
% so that is target dev relative to target pitch, as function of nontarget
% density [default =0];
flanktime_targ = 2.5; % minutes (will do +/- this number). This only used if

useMeanFFRef = 0; % 1: uses mean in window; 0; uses value for each rend

collectTarg =1;
tflankplot = 0.5;
twind_plot = [xedges(1) xedges(end)]./60;

TrialStruct = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S5(TrialStruct, ...
    ignoreLMANexpt, tflankplot, twind_plot, flanktime_targ, flipTargNontarg, ...
    useMeanFFRef, collectTarg, songasrend);


%% ====== TO COLLECT DATA TO PLOT (NEW VERSION, COLLECTS FIRST, THEN PLOTS)
% Instead of collecitng and plotting at same time

% ======== what data to take?
% dotrain = 1; % during training?
% singleRendOnly=1;
densitymethod = 'refperiod';
cutoffmethod = 'medianslice';
mintime_fromref = 5;
maxtime_fromref = 30; % minutes
includeTarg=1; % IMPORTANT: will need to parse by syl type in all foture analyses

% ================ RUN

DATBYREND = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse_S2(TrialStruct, ...
    singleRendOnly, densitymethod, cutoffmethod, mintime_fromref, ...
    maxtime_fromref, ignoreLMANexpt, includeTarg);


assert(all(isnan(DATBYREND.IsDurTrain) == isnan(DATBYREND.Density_isHigh)), 'then nan likely matches up');
assert(all(isnan(DATBYREND.Density_targ) == isnan(DATBYREND.Density_nontarg)), 'then nan likely matches up');


%% get within day time (0-24 hours)
DATBYREND.Tvals_withinday24hours = (DATBYREND.Tvals - floor(DATBYREND.Tvals))*24;

%% ============== [only keep if within time window]
if ~isempty(timewindow_refsyl)
    disp("restricting to within time window");
    % convert tvals (whcih are like 1.5 for noon of day 1) to 24 hr times
    DATBYREND.Tvals_withinday24hours = (DATBYREND.Tvals - floor(DATBYREND.Tvals))*24;
    indstmp = DATBYREND.Tvals_withinday24hours >= timewindow_refsyl(1) & DATBYREND.Tvals_withinday24hours<=timewindow_refsyl(2);
    DATBYREND = rmfield(DATBYREND, 'IsWN');
    DATBYREND = lt_structure_subsample_all_fields(DATBYREND, indstmp, 1);
end


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
title('num dat per rend')

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
onlyPlotFromGenStruct = 1; % i.e. ignore those from seqdeppitch

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

% --------- 4) only from genstruct



% --------- FINAL, convert to numerical indeices
indstoplot = find(indstoplot)';



%% =================== [PLOT] SHOW SAMPLE SIZE
if ~isempty(indstoplot)

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


end


%% ================= [SUMAMRY PLOTS] - GOOD
% -------------- INPUTPARAMS
getsame = 1;
gettrain = 1;
getsiglearn = 0;
birdstoget = 10;
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
    
    disp(tthis)
  
    plot(tthis, ffthis, 'x', 'Color',pcol);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % --------------- COLLECT FOR PLOT
    t_all = [t_all; tthis];
end

try
    linkaxes(hsplots, 'xy');
catch err
end


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


try
    linkaxes(hsplots, 'xy');
catch err
end


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

%% ========== 1) extract syls to plot

getsiglearn = 0;

% birdstoget = 13:17;
% birdstoget = [13 17];
birdstoget = 'notSDP';

syltype = 'targ';
% syltype = 'same';
% syltype = 'all';
% syltype = 'diff';

minbaserends = 100; % has at least this many baseline renditions.
minrends = 100;
% minbaserends = 5; % has at least this many baseline renditions.
% minrends = 5;

[Inds_sylcounter, ~, Inds_exptnum, ...
    Inds_IsSame, Inds_IsTarg] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
    syltype, getsiglearn, birdstoget, minbaserends, minrends);


%% =========== plot
close all;
gettrain = 2;
% 0: only base
% 1: only train;
% 2: both train and base, side by side
plotraw = 1;
edgelist = [2]; % list of edges to use
minrendsinbin = 3;
colorscheme = 'learnsig';

% choices: bird; learnsig

% --- smooth
plotsmooth = 1;

% ---- log time units?
logtime = 1;

% ---- hand enter edges (in units of min)
xedges_hand = []; % leave expty if want to automatically
% xedges_hand = [0 0.32 2 xedges(end)]; % leave expty if want to automatically
% xedges_hand = [0 0.32 xedges(end)]; % leave expty if want to automatically
% xedges_hand = [0 1 xedges(end)]; % leave expty if want to automatically
% segment into 3 bins (before value in edgelist above, then split rest by
% median)

[OUTSTRUCT] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_3(DATBYREND, TrialStruct, ...
    Inds_sylcounter, gettrain, plotraw, edgelist, minrendsinbin, colorscheme, ...
    plotsmooth, logtime, xedges_hand);


%% ================== FOR EACH SYL, COMPARE DIFF TYPE TO SAME TYPE.

% ================== 1) get list of experimnts
if (1)
%     % --- DON"T CARE ABOUT GOOD LERANING [DEFAULT]
%     getsiglearn = 0;
%     birdstoget = 'notSDP';
%     syltype = 'all';
%     minrends = 10;
%     minbaserends = 0; % has at least this many baseline renditions.
%     [Inds_sylcounter, Inds_birdnum, Inds_exptnum, ...
%         Inds_IsSame, Inds_IsTarg] = ...
%         lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
%         syltype, getsiglearn, birdstoget, minbaserends, minrends);
    
    % --- DON"T CARE ABOUT GOOD LERANING [USE FOR if not care about empty
    % gaps --  i.e. take all birds]
    getsiglearn = 0;
    birdstoget = 'all_preferNoSDP';
    syltype = 'all';
    minrends = 10;
    minbaserends = 0; % has at least this many baseline renditions.
    [Inds_sylcounter, Inds_birdnum, Inds_exptnum, ...
        Inds_IsSame, Inds_IsTarg] = ...
        lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
        syltype, getsiglearn, birdstoget, minbaserends, minrends);
else
    % ---- ONLY EXTRACT NOTARG THAT ARE GOOD LEARNING
    birdstoget = 'notSDP';
    minrends = 0;
    % targ
    getsiglearn = 0;
    syltype = 'targ';
    
    [Inds_sylcounter1, Inds_birdnum1, Inds_exptnum1, ...
        Inds_IsSame1, Inds_IsTarg1] = ...
        lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
        syltype, getsiglearn, birdstoget, minbaserends, minrends);
    
    % same
    getsiglearn = 1;
    syltype = 'same';
    
    [Inds_sylcounter2, Inds_birdnum2, Inds_exptnum2, ...
        Inds_IsSame2, Inds_IsTarg2] = ...
        lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
        syltype, getsiglearn, birdstoget, minbaserends, minrends);
    
    % diff
    getsiglearn = 0;
    syltype = 'diff';
    
    [Inds_sylcounter3, Inds_birdnum3, Inds_exptnum3, ...
        Inds_IsSame3, Inds_IsTarg3] = ...
        lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
        syltype, getsiglearn, birdstoget, minbaserends, minrends);
    
    Inds_sylcounter = [Inds_sylcounter1 Inds_sylcounter2 Inds_sylcounter3];
    Inds_birdnum = [Inds_birdnum1 Inds_birdnum2 Inds_birdnum3];
    Inds_exptnum = [Inds_exptnum1 Inds_exptnum2 Inds_exptnum3];
    Inds_IsSame = [Inds_IsSame1 Inds_IsSame2 Inds_IsSame3];
    Inds_IsTarg = [Inds_IsTarg1 Inds_IsTarg2 Inds_IsTarg3];
end


% ================= 3) EXTRACT DATA FOR ALL SYLS
% DEFAULT:
% gettrain = 1;
% plotraw = 0;
% edgelist = [2]; % list of edges to use. each one goes in edges
% %     = [0 edgelist(i) median_of_rest last];
% minrendsinbin = 4;
% colorscheme = 'learnsig';
% plotsmooth = 0;
% logtime = 1;
% xedges_hand = []; % leave expty if want to automatically
% % xedges_hand = [0 1 5 xedges(end)]; % leave expty if want to automatically
% normtobase = 0; % 2: then norms to catch
% use2bins = 0; % then will take only "early" and "late". if 0. then takes 3 bins .

% MODIFIED (to plot multiple bins, hand entered edges):
% gettrain = 1;
% plotraw = 0;
% edgelist = [2]; % list of edges to use
% minrendsinbin = 1;
% colorscheme = 'learnsig';
% plotsmooth = 0;
% logtime = 1;
% xedges_hand = [0 2 5 10 20]; % leave expty if want to automatically [should be [0 stuff]
% % xedges_hand = [0 5 10 20 60]; % leave expty if want to automatically [should be [0 stuff]
% % xedges_hand = [0 2 20]; % leave expty if want to automatically [should be [0 stuff]
% normtobase = 0; % 2: then norms to catch
% use2bins = 0; % then will take only "early" and "late". if 0. then takes 3 bins .
% if ~isempty(xedges_hand)
%     edgelist = [2]; % doesn't matter, but should be length 2, since will iterate, but will overwrite the actual value
% end

% MODIFIED - for getting all datapoints (don't care if there is gap)
gettrain = 1;
plotraw = 0;
edgelist = [2]; % list of edges to use
minrendsinbin = 3;
colorscheme = 'learnsig';
plotsmooth = 0;
logtime = 0;
xedges_hand = [0 5 10 20 40 120]; % leave expty if want to automatically [should be [0 stuff]
xedges_hand = [0 2 5 10 20 60 120]; % leave expty if want to automatically [should be [0 stuff]
% xedges_hand = [0 10 20 60]; % leave expty if want to automatically [should be [0 stuff]
% xedges_hand = [0 5 10 20 60]; % leave expty if want to automatically [should be [0 stuff]
% xedges_hand = [0 2 20]; % leave expty if want to automatically [should be [0 stuff]
normtobase = 0; % 1, minus base. 2: then norms to catch
use2bins = 0; % then will take only "early" and "late". if 0. then takes 3 bins .
if ~isempty(xedges_hand)
    edgelist = [2]; % doesn't matter, but should be length 2, since will iterate, but will overwrite the actual value
end


% [for looking at first 2 hours]
gettrain = 1;
plotraw = 0;
edgelist = [2]; % list of edges to use
minrendsinbin = 5;
colorscheme = 'learnsig';
plotsmooth = 0;
logtime = 0;
xedges_hand = [0 30 60 90 120]; % leave expty if want to automatically [should be [0 stuff]
% xedges_hand = [0 10 20 60]; % leave expty if want to automatically [should be [0 stuff]
% xedges_hand = [0 5 10 20 60]; % leave expty if want to automatically [should be [0 stuff]
% xedges_hand = [0 2 20]; % leave expty if want to automatically [should be [0 stuff]
normtobase = 0; % 1, minus base. 2: then norms to catch
use2bins = 0; % then will take only "early" and "late". if 0. then takes 3 bins .
if ~isempty(xedges_hand)
    edgelist = [2]; % doesn't matter, but should be length 2, since will iterate, but will overwrite the actual value
end

[OUTSTRUCT] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_3(DATBYREND, TrialStruct, ...
    Inds_sylcounter, gettrain, plotraw, edgelist, minrendsinbin, colorscheme, ...
    plotsmooth, logtime, xedges_hand, normtobase, use2bins);

assert(size(OUTSTRUCT.FFdevMean_binned,1) == length(Inds_birdnum))
%% ================= 2) iterate over all experiments
nbirds = max(Inds_birdnum);
nexpts = max(Inds_exptnum);

figcount=1;
subplotrows=3;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

FFdevMeanAll = struct;
TdevMeanAll = struct;
LearnDirAll = [];
count = 1;

for i = 1:nbirds
    for ii=1:nexpts
        
        % ============ SYLS FOR THIS EXPT
        if ~any(Inds_birdnum==i & Inds_exptnum==ii)
            continue
        end
        
        % =============== FOR THIS EXPERIMENT, PLOT TARG, SAME, DIFF
        bname = TrialStruct.birds(i).birdname;
        ename = TrialStruct.birds(i).exptnum(ii).exptname;
        learndir = TrialStruct.birds(i).exptnum(ii).targlearndir;
        
        % ----------- target
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['targ']);
        ylabel([bname '-' ename]);
        
        pcol = 'k';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==1;
        
        t = OUTSTRUCT.Tmean_binned(indsthis, :);
        ff = OUTSTRUCT.FFdevMean_binned(indsthis, :);
        ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
        
        lt_plot(t, ff, {'Errors', ffsem, 'Color', pcol});
        lt_plot_zeroline;
        lt_plot_text(t(2), ff(2)+15, ['learndir:' num2str(learndir)], 'm');
        
        % ----------- SAME
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['same']);
        ylabel([bname '-' ename]);
        
        pcol = 'b';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==1;
        
        t = OUTSTRUCT.Tmean_binned(indsthis, :);
        ff = OUTSTRUCT.FFdevMean_binned(indsthis, :);
        ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
        
        lt_plot(t, ff, {'Errors', ffsem, 'Color', pcol});
        lt_plot_zeroline;
        
        
        % --------------- DIFF
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['diff']);
        ylabel([bname '-' ename]);
        
        pcol = 'r';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==0;
        
        t = OUTSTRUCT.Tmean_binned(indsthis, :);
        ff = OUTSTRUCT.FFdevMean_binned(indsthis, :);
        ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
        
        lt_plot(t, ff, {'Errors', ffsem, 'Color', pcol});
        lt_plot_zeroline;
        
        
        % ----------- OVERLAID
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['ALL']);
        ylabel([bname '-' ename]);
        
        % targ
        pcol = 'k';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==1;
        t = OUTSTRUCT.Tmean_binned(indsthis, :);
        ff = OUTSTRUCT.FFdevMean_binned(indsthis, :);
        ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
        plot(t, ff, 'o', 'Color', pcol);
        plot(nanmean(t,1), nanmean(ff,1), 'LineWidth', 2, 'Color', pcol);
        FFdevMeanAll(count).targ = nanmean(ff,1)';
        TdevMeanAll(count).targ = nanmean(t, 1)';
        
        % ----------- SAME
        pcol = 'b';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==1;
        
        t = OUTSTRUCT.Tmean_binned(indsthis, :);
        ff = OUTSTRUCT.FFdevMean_binned(indsthis, :);
        ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
        plot(t, ff, 'o', 'Color', pcol);
        plot(nanmean(t,1), nanmean(ff,1), 'LineWidth', 2, 'Color', pcol);
        FFdevMeanAll(count).same = nanmean(ff,1)';
        TdevMeanAll(count).same = nanmean(t, 1)';
        
        % --------------- DIFF
        pcol = 'r';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==0;
        
        t = OUTSTRUCT.Tmean_binned(indsthis, :);
        ff = OUTSTRUCT.FFdevMean_binned(indsthis, :);
        ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
        plot(t, ff, 'o', 'Color', pcol);
        plot(nanmean(t,1), nanmean(ff,1), 'LineWidth', 2, 'Color', pcol);
        FFdevMeanAll(count).diff = nanmean(ff,1)';
        TdevMeanAll(count).diff = nanmean(t, 1)';
        
        % ---
        lt_plot_zeroline;
        
        % ============== COLLECT TARG LEARN DIR
        LearnDirAll = [LearnDirAll; learndir];
        
        
        % ====================== OVERLAY SMOOTHED FUNCTIONS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['smoothed']);
        
        % targ
        pcol = 'k';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==1;
        
        for j=find(indsthis)
            if isempty(OUTSTRUCT.Tsm_alltrials{j})
                continue
            end
            t = OUTSTRUCT.Tsm_alltrials{j}.Mean;
            ff = OUTSTRUCT.FFsm_alltrials{j}.Mean;
            plot(t, ff, 'Color', pcol);
        end
        % same
        pcol = 'b';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==1;
        
        for j=find(indsthis)
            if isempty(OUTSTRUCT.Tsm_alltrials{j})
                continue
            end
            t = OUTSTRUCT.Tsm_alltrials{j}.Mean;
            ff = OUTSTRUCT.FFsm_alltrials{j}.Mean;
            plot(t, ff, 'Color', pcol);
        end
        
        % diff
        pcol = 'r';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==0;
        
        for j=find(indsthis)
            if isempty(OUTSTRUCT.Tsm_alltrials{j})
                continue
            end
            t = OUTSTRUCT.Tsm_alltrials{j}.Mean;
            ff = OUTSTRUCT.FFsm_alltrials{j}.Mean;
            plot(t, ff, 'Color', pcol);
        end
        
        
        % ========================= PLOT DISTRIBUTIONS OF TIME DEVIATIONS
        % targ
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %         hsplots = [hsplots hsplot];
        title(['n adjacent rends']);
        
        pcol = 'k';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==1;
        
        tdevs = cell2mat(OUTSTRUCT.Tdev_alltrials(indsthis));
        if ~isempty(tdevs)
            [~, Xcenters] = lt_plot_histogram(tdevs, '' , 1, 1, '', 1, pcol);
        end
        
        % SAME
        pcol = 'b';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==1;
        
        tdevs = cell2mat(OUTSTRUCT.Tdev_alltrials(indsthis));
        if ~isempty(tdevs)
            lt_plot_histogram(tdevs, Xcenters , 1, 1, '', 1, pcol);
        end
        
        % DIFF
        pcol = 'r';
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==0;
        
        tdevs = cell2mat(OUTSTRUCT.Tdev_alltrials(indsthis));
        if ~isempty(tdevs)
            lt_plot_histogram(tdevs, Xcenters , 1, 1, '', 1, pcol);
        end
        
        %- --- overlay edges
        t_edges = median(OUTSTRUCT.Xedgesall(indsthis, :),1);
        for j=1:length(t_edges)
            line([t_edges(j) t_edges(j)], ylim);
        end
        
        
        % ================= COLLECT
        count = count+1;
        
    end
end


linkaxes(hsplots, 'xy');




% ============== SUMMARY PLOT, OVERLAY ALL TARG, SAME, DIFF.
lt_figure; hold on;

% ------- each expt
% lt_subplot(3,2,1); hold on;
lt_subplot(2,4,1); hold on;
title(['means of syls']);
xlabel('bin center (time)');
ylabel('mean FFdev');

% --- targ
pcol = 'k';
fname = 'targ';
for j=1:length(FFdevMeanAll)
    x = TdevMeanAll(j).(fname);
    y = FFdevMeanAll(j).(fname);
    plot(x,y,'-', 'Color', pcol);
end
ymean = nanmean([FFdevMeanAll.(fname)],2);
xmean = nanmean([TdevMeanAll.(fname)],2);
ysem = lt_sem([FFdevMeanAll.(fname)]');
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', pcol});
% -- test if deviate from bin1
tmp = [FFdevMeanAll.(fname)];
for j=2:size(tmp,1)
    if all(isnan(tmp(j,:)))
        continue
    end
   p = signrank(tmp(j,:), tmp(1,:));
   if p<0.2
   lt_plot_text(xmean(j), ymean(j), ['(vsbin1)p=' num2str(p)], pcol);
   end
end
p = signrank(nanmean(tmp(2:end,:)), tmp(1,:));
if p<0.2
    lt_plot_text(nanmean(xmean(2:end)), nanmean(ymean(2:end))+20, ['(vsbin1)p=' num2str(p)], pcol);
end
% -- test if deviate from 0
for j=2:size(tmp,1)
    if all(isnan(tmp(j,:)))
        continue
    end
   p = signrank(tmp(j,:), 0);
   if p<0.2
   lt_plot_text(xmean(j), ymean(j), ['(vs0)p=' num2str(p)], pcol);
   end
end
p = signrank(nanmean(tmp(2:end,:)), 0);
if p<0.2
    lt_plot_text(nanmean(xmean(2:end)),  nanmean(ymean(2:end))+20, ['(vs0)p=' num2str(p)], pcol);
end
% ---

% ---- same
pcol = 'b';
fname = 'same';
for j=1:length(FFdevMeanAll)
    x = TdevMeanAll(j).(fname);
    y = FFdevMeanAll(j).(fname);
    plot(x,y,'-', 'Color', pcol);
end
ymean = nanmean([FFdevMeanAll.(fname)],2);
xmean = nanmean([TdevMeanAll.(fname)],2);
ysem = lt_sem([FFdevMeanAll.(fname)]');
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', pcol});
% -- test if deviate from bin1
tmp = [FFdevMeanAll.(fname)];
for j=2:size(tmp,1)
    if all(isnan(tmp(j,:)))
        continue
    end
   p = signrank(tmp(j,:), tmp(1,:));
   if p<0.2
   lt_plot_text(xmean(j), ymean(j), ['(vsbin1)p=' num2str(p)], pcol);
   end
end
p = signrank(nanmean(tmp(2:end,:)), tmp(1,:));
if p<0.2
    lt_plot_text(nanmean(xmean(2:end)), nanmean(ymean(2:end))+20, ['(vsbin1)p=' num2str(p)], pcol);
end
% -- test if deviate from 0
for j=2:size(tmp,1)
    if all(isnan(tmp(j,:)))
        continue
    end
   p = signrank(tmp(j,:), 0);
   if p<0.2
   lt_plot_text(xmean(j), ymean(j), ['(vs0)p=' num2str(p)], pcol);
   end
end
p = signrank(nanmean(tmp(2:end,:)), 0);
if p<0.2
    lt_plot_text(nanmean(xmean(2:end)),  nanmean(ymean(2:end))+20, ['(vs0)p=' num2str(p)], pcol);
end
% ---

% ---- diff
pcol = 'r';
fname = 'diff';
for j=1:length(FFdevMeanAll)
    x = TdevMeanAll(j).(fname);
    y = FFdevMeanAll(j).(fname);
    plot(x,y,'-', 'Color', pcol);
end
ymean = nanmean([FFdevMeanAll.(fname)],2);
xmean = nanmean([TdevMeanAll.(fname)],2);
ysem = lt_sem([FFdevMeanAll.(fname)]');
lt_plot(xmean, ymean, {'Errors', ysem, 'Color', pcol});
% -- test if deviate from bin1
tmp = [FFdevMeanAll.(fname)];
for j=2:size(tmp,1)
    if all(isnan(tmp(j,:)))
        continue
    end
   p = signrank(tmp(j,:), tmp(1,:));
   if p<0.2
   lt_plot_text(xmean(j), ymean(j), ['(vsbin1)p=' num2str(p)], pcol);
   end
end
p = signrank(nanmean(tmp(2:end,:)), tmp(1,:));
if p<0.2
    lt_plot_text(nanmean(xmean(2:end)), nanmean(ymean(2:end))+20, ['(vsbin1)p=' num2str(p)], pcol);
end
% -- test if deviate from 0
for j=2:size(tmp,1)
    if all(isnan(tmp(j,:)))
        continue
    end
   p = signrank(tmp(j,:), 0);
   if p<0.2
   lt_plot_text(xmean(j), ymean(j), ['(vs0)p=' num2str(p)], pcol);
   end
end
p = signrank(nanmean(tmp(2:end,:)), 0);
if p<0.2
    lt_plot_text(nanmean(xmean(2:end)),  nanmean(ymean(2:end))+20, ['(vs0)p=' num2str(p)], pcol);
end
% ---

% ---
axis tight;
lt_plot_zeroline;


% ============ COMPARE TARG TO SAME
% lt_subplot(3,2,2); hold on;
% ylabel('same - targ (hz)');
% title('each expt, mean over syls');
% xlabel('time dev');
% 
% y = [FFdevMeanAll.same] - [FFdevMeanAll.targ];
% x = [TdevMeanAll.same];
% 
% plot(x, y, '-ob');
% lt_plot_zeroline;
% 
% % --- significance
% lt_subplot(3,2,2); hold on;
lt_subplot(2,4,2); hold on;
ylabel('targ - same (hz)');
title('each expt, mean over syls');
xlabel('time dev');

% y = [FFdevMeanAll.targ] - [FFdevMeanAll.same];
y = [FFdevMeanAll.same] - [FFdevMeanAll.targ];
x = [TdevMeanAll.same];

plot(x, y, '-ob');
lt_plot_zeroline;

ymean = nanmean(y,2);
ysem = lt_sem(y');
x = nanmean(x,2);
lt_plot_bar(x, ymean', {'Errors', ysem, 'FaceAlpha', 0.2})

% --- significance
% - for each time bin
for j=1:size(y,1)
    if all(isnan(y(j,:)))
        continue
    end
%    [~, p]= ttest(y(j,:));
   [p]= signrank(y(j,:));
       lt_plot_text(x(j), max(y(j,:)), ['p(vs0)=' num2str(p)], 'm', 8);
end
% -- compare each bin to 1
for j=2:size(y,1)
    if all(isnan(y(j,:)))
        continue
    end
%    [~, p]= ttest(y(j,:));
   [p]= signrank(y(1,:), y(j,:));
       lt_plot_text(x(j), 1.1*max(y(j,:)), ['p(vsbin1)=' num2str(p)], 'r', 8);
end
p = signrank(y(1,:), nanmean(y([2:3],:),1));
lt_plot_text(x(2), y(2), ['(Bin2/3vsBin1)p=' num2str(p)], 'c', 8);



% =========== COMBINE LAST 2 BINS INTO ONE
% lt_subplot(3,2,5); hold on;
lt_subplot(2,4,5); hold on;
ylabel('same - targ (hz)');
title('each expt, mean over syls');
xlabel('time dev');

y = [FFdevMeanAll.same] - [FFdevMeanAll.targ];
ytmp = [y(1,:); nanmean(y([2 3], :))];
x = [TdevMeanAll.same];
xtmp = [x(1,:); nanmean(x([2 3], :))];

plot(xtmp, ytmp, '-ob');
lt_plot_zeroline;
% --pval = 
for j=1:size(ytmp,1)
    if all(isnan(ytmp(j,:)))
        continue
    end
%    [~, p]= ttest(y(j,:));
   [p]= signrank(ytmp(j,:));
       lt_plot_text(x(j), max(ytmp(j,:)), ['p(vs0)=' num2str(p)], 'm', 8);
end
p = signrank(ytmp(1,:), ytmp(2,:));
lt_plot_text(x(2), ytmp(2), ['(Bin2vsBin1)p=' num2str(p)], 'r', 8);


% ========== FOR BOTH TARG AND SAME, SUIBNTRACT DIFF TYPE
% lt_subplot(3,2,4); hold on;
lt_subplot(2,4,4); hold on;
title('after minus diff type');

% --- targ
ytarg = [FFdevMeanAll.targ] - [FFdevMeanAll.diff];
xtarg = [TdevMeanAll.targ];

plot(xtarg, ytarg, '-k');
lt_plot(nanmean(xtarg,2), nanmean(ytarg, 2), {'Errors', lt_sem(ytarg'), 'Color', 'k'});

% -- same
ysame = [FFdevMeanAll.same] - [FFdevMeanAll.diff];
xsame = [TdevMeanAll.same];

plot(xsame, ysame, '-b');
lt_plot(nanmean(xsame,2), nanmean(ysame, 2), {'Errors', lt_sem(ysame'), 'Color', 'b'});
lt_plot_zeroline;



% =========== ACCOUNT FOR TOTAL LEARNING -
% lt_subplot(3,2,3); hold on;
lt_subplot(2,4,3); hold on;
title('k=targ; bu=same');
xlabel('ffdev (bin1)');
ylabel('ffdev (bin2)');

% ---- REGRESS CHANGE IN BIN 1 VS CHANGE IN BIN 2
bin1_inds = 1;
if use2bins==1
    bin2_inds=2;
else
    bin2_inds = 2:3; % average over if multiple
end

% --------------- TARG
pcol = 'k';
ffmat = [FFdevMeanAll.targ];

ff_bin1 = mean(ffmat(bin1_inds,:), 1);
ff_bin2 = mean(ffmat(bin2_inds, :), 1);
% plot(ff_bin1, ff_bin2, 'o', 'Color', pcol);
lt_plot(ff_bin1, ff_bin2, {'Color', pcol})

% --------------- SAME
pcol = 'c';
ffmat = [FFdevMeanAll.same];

ff_bin1 = mean(ffmat(bin1_inds,:), 1);
ff_bin2 = mean(ffmat(bin2_inds, :), 1);
lt_plot(ff_bin1, ff_bin2, {'Color', pcol});

% ---- line connecting
for j=1:length(FFdevMeanAll)
    
    x1 = mean(FFdevMeanAll(j).targ(bin1_inds));
    x2 = mean(FFdevMeanAll(j).same(bin1_inds));
    
    y1 = mean(FFdevMeanAll(j).targ(bin2_inds));
    y2 = mean(FFdevMeanAll(j).same(bin2_inds));
    
    line([x1 x2], [y1 y2], 'Color', 'm');
end

lt_plot_makesquare_plot45line(gca, 'r');



% ====================================
% lt_subplot(3,2,6); hold on;
lt_subplot(2,4,6); hold on;
title('deviations, late minus early');
xlabel('targ --- same');

Y = {};

% --------------- TARG
ffmat = [FFdevMeanAll.targ];
% ffmat = [FFdevMeanAll.targ] - [FFdevMeanAll.diff];

ff_bin1 = nanmean(ffmat(bin1_inds,:), 1);
ff_bin2 = nanmean(ffmat(bin2_inds, :), 1);

Y{1} = [ff_bin2 - ff_bin1]';


% --------------- SAME
ffmat = [FFdevMeanAll.same];
% ffmat = [FFdevMeanAll.same] - [FFdevMeanAll.diff];

ff_bin1 = nanmean(ffmat(bin1_inds,:), 1);
ff_bin2 = nanmean(ffmat(bin2_inds, :), 1);

Y{2} = [ff_bin2 - ff_bin1]';

% --- plot
x = [1 2];
Y = cell2mat(Y);
plot(x, Y, '-k')
lt_plot(x+0.2, mean(Y,1), {'Errors', lt_sem(Y)});
xlim([0 3]);
lt_plot_zeroline;

p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'vs', 1);

%% ====================== PLOT ALL SAME TYPE SYLS (BY SYL NOT BUT EXPT)

FFall = [];
FFall_targMinusSame = [];
Tall = [];

% tally to get N.
b_targ=[];
s_targ=[];
b_same =[];
s_same = [];

for i = 1:nbirds
    for ii=1:nexpts
        
        % ============ SYLS FOR THIS EXPT
        if ~any(Inds_birdnum==i & Inds_exptnum==ii)
            continue
        end
        
        % =============== FOR THIS EXPERIMENT, PLOT TARG, SAME, DIFF
        bname = TrialStruct.birds(i).birdname;
        ename = TrialStruct.birds(i).exptnum(ii).exptname;
        
        
        %% SAME MINUS DIFF
        % =============== FIRST, AVERAGE OVER ALL THE DIFF TYPES IN ORDER
        % TO NORMALIZE
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==0;
        
        ffnorm = nanmean(OUTSTRUCT.FFdevMean_binned(indsthis, :), 1);
        
        
        % ============= for each same type, collect after subtracting diff type mean
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==1;
        
        for j = find(indsthis)
            t = OUTSTRUCT.Tmean_binned(j, :);
            ff = OUTSTRUCT.FFdevMean_binned(j, :);
            %         ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
            
            % ------- subtract diff tyupe
            ff = ff - ffnorm;
            FFall = [FFall; ff];
            
        end
        
        %% TARG MINUS SAME
        % ============= SAME, FOR TARG MINUS SAME
                % =============== FIRST, AVERAGE OVER ALL THE DIFF TYPES IN ORDER
        % TO NORMALIZE
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==0 & Inds_IsSame==1;
        
        ffnorm = nanmean(OUTSTRUCT.FFdevMean_binned(indsthis, :), 1);
        
        
        % ============= for each same type, collect after subtracting diff type mean
        indsthis = Inds_birdnum==i & Inds_exptnum==ii & ...
            Inds_IsTarg==1 & Inds_IsSame==1;
        
        for j = find(indsthis)
%             t = OUTSTRUCT.Tmean_binned(j, :);
            ff = OUTSTRUCT.FFdevMean_binned(j, :);
            %         ffsem = OUTSTRUCT.FFdevSEM_binned(indsthis, :);
            
            % ------- subtract diff tyupe
            ff = ff - ffnorm;
            FFall_targMinusSame = [FFall_targMinusSame; ff];
            
        end

    end
end

lt_figure; hold on;
% ====== % all same tupe syls
lt_subplot(3,2,1); hold on;
title('targ');
ylabel('hz');

Y = OUTSTRUCT.FFdevMean_binned(Inds_IsSame==1 & Inds_IsTarg==1, :);

x = 1:size(Y,2);
plot(x, Y', '-k');

lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y)});
lt_plot_zeroline;
% - stats
for j=1:size(Y,2)
    if sum(~isnan(Y(:,j)))>1
       p = signrank(Y(:,j));
    if p<0.1
        lt_plot_text(x(j), max(Y(:,j)), ['(vs0)p=' num2str(p)], 'm', 8)
    end
    end
end
B = Inds_birdnum(Inds_IsSame==1 & Inds_IsTarg==1);
E = Inds_exptnum(Inds_IsSame==1 & Inds_IsTarg==1);
S = Inds_sylcounter(Inds_IsSame==1 & Inds_IsTarg==1);
Nbirds = length(unique(B));
Nsyls = length(unique(S));
Nexpt = length(unique(lt_tools_grp2idx({B', E'})));
lt_plot_annotation(1, ["Nbirds=" num2str(Nbirds) ", Nexpt=" num2str(Nexpt) ",Nsyls=" num2str(Nsyls)]);

% ====== % all same tupe syls
lt_subplot(3,2,2); hold on;
title('same');
ylabel('hz');

Y = OUTSTRUCT.FFdevMean_binned(Inds_IsSame==1 & Inds_IsTarg==0, :);
x = 1:size(Y,2);
plot(x, Y', '-b');

lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y)});

lt_plot_zeroline;
% - stats
for j=1:size(Y,2)
    if sum(~isnan(Y(:,j)))>1
       p = signrank(Y(:,j));
       [~, p] = ttest(Y(:,j));
    if p<0.1
        lt_plot_text(x(j), max(Y(:,j)), ['(vs0)p=' num2str(p)], 'm', 8)
    end
    end
end
B = Inds_birdnum(Inds_IsSame==1 & Inds_IsTarg==0);
E = Inds_exptnum(Inds_IsSame==1 & Inds_IsTarg==0);
S = Inds_sylcounter(Inds_IsSame==1 & Inds_IsTarg==0);
Nbirds = length(unique(B));
Nsyls = length(unique(S));
Nexpt = length(unique(lt_tools_grp2idx({B', E'})));
lt_plot_annotation(1, ["Nbirds=" num2str(Nbirds) ", Nexpt=" num2str(Nexpt) ",Nsyls=" num2str(Nsyls)]);


% ====== % all diff tupe syls
lt_subplot(3,2,3); hold on;
title('diff');
ylabel('hz');

Y = OUTSTRUCT.FFdevMean_binned(Inds_IsSame==0 & Inds_IsTarg==0, :);
x = 1:size(Y,2);
plot(x, Y', '-r');
lt_plot(x+0.2, nanmean(Y,1), {'Errors', lt_sem(Y)});
lt_plot_zeroline;
% - stats
for j=1:size(Y,2)
    if sum(~isnan(Y(:,j)))>1
       p = signrank(Y(:,j));
    if p<0.1
        lt_plot_text(x(j), max(Y(:,j)), ['(vs0)p=' num2str(p)], 'm', 8)
    end
    end
end
B = Inds_birdnum(Inds_IsSame==0 & Inds_IsTarg==0);
E = Inds_exptnum(Inds_IsSame==0 & Inds_IsTarg==0);
S = Inds_sylcounter(Inds_IsSame==0 & Inds_IsTarg==0);
Nbirds = length(unique(B));
Nsyls = length(unique(S));
Nexpt = length(unique(lt_tools_grp2idx({B', E'})));
lt_plot_annotation(1, ["Nbirds=" num2str(Nbirds) ", Nexpt=" num2str(Nexpt) ",Nsyls=" num2str(Nsyls)]);


% ====== 2) subtract diff type
lt_subplot(3,2,4); hold on;
title('same, subtract diff');
ylabel('hz');

x = 1:size(FFall,2);
plot(x, FFall', '-b');
ymean = nanmean(FFall, 1);
ysem = lt_sem(FFall);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'b'});
lt_plot_zeroline;

for ii = 1:size(FFall,2)
    if all(isnan(FFall(:, ii)))
        continue
    end
    pp = signrank(FFall(:,ii));
    lt_plot_text(x(ii), 2*ymean(ii), ['p=' num2str(pp)], 'm', 8);    
end



% ====== 2) subtract diff type
lt_subplot(3,2,5); hold on;
title('targ, subtract same');
ylabel('hz');

x = 1:size(FFall_targMinusSame,2);
plot(x, FFall_targMinusSame', '-k');
ymean = nanmean(FFall_targMinusSame, 1);
ysem = lt_sem(FFall_targMinusSame);
lt_plot(x+0.2, ymean, {'Errors', ysem, 'Color', 'k'});
lt_plot_zeroline;

for ii = 1:size(FFall_targMinusSame,2)
    if all(isnan(FFall_targMinusSame(:, ii)))
        continue
    end
    pp = signrank(FFall_targMinusSame(:,ii));
    lt_plot_text(x(ii), 2*ymean(ii), ['p=' num2str(pp)], 'm', 8);    
end


%% ======= SUBTRACT SAME TYPE FROM TARGET



%% ################## plots, comparing training to baseline
if (0) % IN PROGRESS!!!
    gettrain = 2;
    % 0: only base
    % 1: only train;
    % 2: both train and base, side by side
    plotraw = 1;
    edgelist = [3]; % list of edges to use
    minrendsinbin = 4;
    colorscheme = 'learnsig';
    % choices: bird; learnsig
    
    [OUTSTRUCT] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_5(DATBYREND, TrialStruct, ...
        Inds_sylcounter, gettrain, plotraw, edgelist, minrendsinbin, colorscheme);
end


%% ################ plot all raw data (including those that were filtered out for suymmary above)


plotraw = 1;
% edgelist = [2 3 4 5]; % list of edges to use
edgelist = [3]; % list of edges to use
minrendsinbin = 0;

gettrain = 1;
% 0: only base
% 1: only train;
% 2: both train and base, side by side

getsiglearn = 0;

% birdstoget = 13:17;
% birdstoget = 1;
birdstoget = [];

colorscheme = 'learnsig';
% choices: bird; learnsig

%     plotraw = 0;
%     edgelist = [2 3 4 5]; % list of edges to use

% syltype = 'same';
syltype = 'all';
% syltype = 'diff';

lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_1(DATBYREND, TrialStruct, syltype, gettrain, ...
    getsiglearn, birdstoget, plotraw, edgelist, minrendsinbin, colorscheme);

%% ############### 2) plot
% I.E. for each syl, get divide into first and second half, matching sample
% sizes
% 3 bins: 0-4min, the rest divided into 2

plotraw = 1;
% edgelist = [2 3 4 5]; % list of edges to use
edgelist = [3]; % list of edges to use
minrendsinbin = 4;

gettrain = 1;
% 0: only base
% 1: only train;
% 2: both train and base, side by side

getsiglearn = 0;

% birdstoget = 13:17;
% birdstoget = 1;
birdstoget = 'notSDP';

colorscheme = 'learnsig';
% choices: bird; learnsig

%     plotraw = 0;
%     edgelist = [2 3 4 5]; % list of edges to use

% syltype = 'same';
syltype = 'same';
% syltype = 'diff';

OUTSTRUCT = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_1(DATBYREND, TrialStruct, syltype, gettrain, ...
    getsiglearn, birdstoget, plotraw, edgelist, minrendsinbin, colorscheme);




%% ==================== [SUMMARY PLOT, INPUT NEW BINEDGES]
% I.E. for each syl, get divide into first and second half, matching sample
% sizes
% 3 bins: 0-4min, the rest divided into 2

% =====================
xedgethis = [0:4:32];
plotAllBirds = 1; % then plots everything in DATBYREND

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
    
    if plotAllBirds~=1
    if ~isempty(birdstoget)
        if ismember(bnum, birdstoget)==0
            continue
        end
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



%% ================================================
%% ======================== WN FEEDBACK/CATCH, EFFECTS

getsiglearn = 0;

% birdstoget = 13:17;
% birdstoget = [13 17];
birdstoget = 'notSDP';

syltype = 'targ';
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


%% ===============

figcount=1;
subplotrows=4;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

istrain = 1;
skipifnocatch = 0;
plotlogtime = 1;

binrun = 25; % smooth

pcolAll = {'r', 'b', 'm', 'c'};
pmarkAll = {'x', 'x', 'o', 'o'};

xbinedges = [3]; % in min. just put dividers - will fill in the ends

FFmeanAll = nan(length(Inds_sylcounter), 4);
for i=1:length(Inds_sylcounter)
    ss = Inds_sylcounter(i);
    
    indsthis = DATBYREND.Sylcounter==ss & DATBYREND.IsDurTrain==istrain ...
        & ~cellfun(@isempty, DATBYREND.FF_dev);
    
    if ~any(indsthis)
        continue
    end
    
    ffthis = DATBYREND.FF_dev(indsthis);
    tthis = DATBYREND.Time_dev(indsthis);
    isWN = DATBYREND.IsWN(indsthis);
    isCatch = DATBYREND.IsCatch(indsthis);
    
    if ~any(isWN)
        disp('NO WN!!!')
        pause
        continue
    end
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    %         title(['targ']);
    bname = TrialStruct.birds(unique(DATBYREND.Birdnum(indsthis))).birdname;
    ename = TrialStruct.birds(unique(DATBYREND.Birdnum(indsthis))).exptnum(unique(DATBYREND.Exptnum(indsthis))).exptname;
    title([bname '-' ename]);
    
    % ============== get bin edges for this syl
    xbinedge_this = [min(cell2mat(tthis))*60*24 xbinedges max(cell2mat(tthis))*60*24];
    if plotlogtime ==1
        xbinedge_this = log10(xbinedge_this);
        %         tthis=log10(tthis);
    end
    
    
    % ############################################ PLOT
    YsmAll = cell(1,4);
    TsmAll = cell(1,4);
    
    YbinAll = cell(1,4);
    YbinSem = cell(1,4);
    TbinAll = cell(1,4);
    
    YmeanAll = nan(1,4);
    % ========== 1) WN HIT (NOTCATCH)
    indtmp = isWN==1 & isCatch==0;
    pcol = 'r';
    pmark = 'x';
    x = 1;
    
    ttmp = cell2mat(tthis(indtmp))*60*24;
    fftmp = cell2mat(ffthis(indtmp));
    plot(ttmp, fftmp, pmark, 'Color', pcol);
    
    % --- get smoothed
    if singleRendOnly==1
        [~, indsort] = sort(ttmp);
        ttmp = ttmp(indsort);
        fftmp = fftmp(indsort);
        
        TsmAll{x} = lt_running_stats(ttmp, binrun);
        YsmAll{x} = lt_running_stats(fftmp, binrun);
    else
        TsmAll{x} = [];
        YsmAll{x} = [];
    end
    % --- get binned
    tbins = discretize(ttmp, xbinedge_this);
    [ffbinmean, ffbinsem] = grpstats(fftmp, tbins, {'mean', 'sem'});
    tbinmean = unique(tbins(~isnan(tbins)));
    YbinAll{x} = ffbinmean;
    YbinSem{x} = ffbinsem;
    TbinAll{x} = tbinmean;
    % --- get mean
    YmeanAll(x) = mean(fftmp);
    
    % ========== 2) ESCAPE (NOTCATCH)
    indtmp = isWN==0 & isCatch==0;
    pcol = 'b';
    pmark = 'x';
    x = 2;
    
    ttmp = cell2mat(tthis(indtmp))*60*24;
    fftmp = cell2mat(ffthis(indtmp));
    plot(ttmp, fftmp, pmark, 'Color', pcol);
    
    % --- get smoothed
    if singleRendOnly==1
        [~, indsort] = sort(ttmp);
        ttmp = ttmp(indsort);
        fftmp = fftmp(indsort);
        
        TsmAll{x} = lt_running_stats(ttmp, binrun);
        YsmAll{x} = lt_running_stats(fftmp, binrun);
    else
        TsmAll{x} = [];
        YsmAll{x} = [];
    end
    % --- get binned
    tbins = discretize(ttmp, xbinedge_this);
    [ffbinmean, ffbinsem] = grpstats(fftmp, tbins, {'mean', 'sem'});
    tbinmean = unique(tbins(~isnan(tbins)));
    YbinAll{x} = ffbinmean;
    YbinSem{x} = ffbinsem;
    TbinAll{x} = tbinmean;
    % --- get mean
    YmeanAll(x) = mean(fftmp);
    
    % ============================
    lt_plot_zeroline;
    
    % ##################### CATCH/WN?
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    if any(isCatch)
        
        % ========== 3) WN HIT (CATCH)
        indtmp = isWN==1 & isCatch==1;
        pcol = 'r';
        pmark = 'o';
        x = 3;
        
        ttmp = cell2mat(tthis(indtmp))*60*24;
        fftmp = cell2mat(ffthis(indtmp));
        plot(ttmp, fftmp, pmark, 'Color', pcol);
        
        % --- get smoothed
        if singleRendOnly==1
            [~, indsort] = sort(ttmp);
            ttmp = ttmp(indsort);
            fftmp = fftmp(indsort);
            
            TsmAll{x} = lt_running_stats(ttmp, binrun);
            YsmAll{x} = lt_running_stats(fftmp, binrun);
        else
            TsmAll{x} = [];
            YsmAll{x} = [];
        end
        % --- get binned
        tbins = discretize(ttmp, xbinedge_this);
        [ffbinmean, ffbinsem] = grpstats(fftmp, tbins, {'mean', 'sem'});
        tbinmean = unique(tbins(~isnan(tbins)));
        YbinAll{x} = ffbinmean;
        YbinSem{x} = ffbinsem;
        TbinAll{x} = tbinmean;
        % --- get mean
        YmeanAll(x) = mean(fftmp);
        
        % ========== 4) ESCAPE (CATCH)
        indtmp = isWN==0 & isCatch==1;
        pcol = 'b';
        pmark = 'o';
        
        x = 4;
        
        ttmp = cell2mat(tthis(indtmp))*60*24;
        fftmp = cell2mat(ffthis(indtmp));
        plot(ttmp, fftmp, pmark, 'Color', pcol);
        
        % --- get smoothed
        if singleRendOnly==1
            [~, indsort] = sort(ttmp);
            ttmp = ttmp(indsort);
            fftmp = fftmp(indsort);
            
            TsmAll{x} = lt_running_stats(ttmp, binrun);
            YsmAll{x} = lt_running_stats(fftmp, binrun);
        else
            TsmAll{x} = [];
            YsmAll{x} = [];
        end
        % --- get binned
        tbins = discretize(ttmp, xbinedge_this);
        [ffbinmean, ffbinsem] = grpstats(fftmp, tbins, {'mean', 'sem'});
        tbinmean = unique(tbins(~isnan(tbins)));
        YbinAll{x} = ffbinmean;
        YbinSem{x} = ffbinsem;
        TbinAll{x} = tbinmean;
        % --- get mean
        YmeanAll(x) = mean(fftmp);
        
        
        % ============================
        lt_plot_zeroline;
    end
    FFmeanAll(i, :) = YmeanAll;
    
    % ############################# OVERLAY SMOOTHED FOR ALL TYPES
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    for j=1:length(YsmAll)
        try
            t = TsmAll{j}.Mean;
            ff = YsmAll{j}.Mean;
            ffsem = YsmAll{j}.SEM;
            
            shadedErrorBar(t, ff, ffsem, {'Color', pcolAll{j}},1);
        catch err
        end
    end
    
    % ================= OVERLAY BINNED
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    for j=1:length(YbinAll)
        try
            t = TbinAll{j};
            ff = YbinAll{j};
            ffsem = YbinSem{j};
            
            lt_plot(t, ff, {'Errors', ffsem, 'Color', pcolAll{j}});
        catch
        end
    end
    lt_plot_zeroline;
    xlim([0 4]);
    
    
end

linkaxes(hsplots, 'xy');

% ==================== summary
lt_figure; hold on;


lt_subplot(3,2,1); hold on;
xlabel('HIT(NC) -- ESCAPE(NC) -- HIT(C) -- ESCAPE(NC)');
Y = FFmeanAll(~any(isnan(FFmeanAll)'), :);
x = 1:size(Y,2);
plot(x, Y, '-k');
xlim([0 5]);
lt_plot(x+0.2, mean(Y,1), {'Errors', lt_sem(Y)});

% ----------- subtract catch
lt_subplot(3,2,2); hold on;
xlabel('HIT(NC-C) -- ESCAPE(NC-C)');
Y = FFmeanAll(~any(isnan(FFmeanAll)'), :);
Y = [Y(:,1)-Y(:,3) Y(:,2)-Y(:,4)];
x = 1:size(Y,2);
plot(x, Y, '-k');
xlim([0 5]);
lt_plot(x+0.2, mean(Y,1), {'Errors', lt_sem(Y)});




%% ######################## PREDICTING LEARNING BY LOCAL LEARNING, HIT, AND CATCH
% DATTMP = struct;
% % NOTE: DO NOT CHANGE ORDER!! of following 3 extractions...
% use_targ_locallearn = 1;
% use_nminus2 = 0; % if 1, then use n-2 as local learning. if 0. thne use n-1.
%
% if ~strcmp(syltype, 'targ')
%     % then should use minus 1, so deviations will be agaiunst your own syl
%     use_nminus2 = 0;
% end
%
% % ------------------ TRAINING, NOTCATCH
istrain = 1;
iscatch = 0;
DATTMP.Train = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
    TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn);

% % ------------------ TRAINING, CATCH
% istrain = 1;
% iscatch = 1;
% DATTMP.Catch = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
%     TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn);
%
% % ------------------ BASELINE, C/NC
% istrain = 0;
% iscatch = [0 1];
% DATTMP.Base = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_6(DATBYREND, ...
%     TrialStruct, Inds_sylcounter, istrain, iscatch, use_targ_locallearn);
%
% % ------------------ use the local learning for the syllable of interest,
% % or for the target syllable in that experiment?
%
%
%
% %% ########################### PLOT (absolute balue of local earning - i.e.
% % more dwivation from mmedian predict more laernign?)
%
% figcount=1;
% subplotrows=4;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
%
% takeabs = 1; % takes absolute value of local learn (i.e. might predict strong increase in next
% % trial when current trial is more deviated from local.
% fnames = fieldnames(DATTMP);
% for j=1:length(DATTMP.Train)
%
%     for fname = fnames'
%         fname = fname{1};
%         % ================
%         [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%         hsplots = [hsplots hsplot];
% %         title([bname '-' ename ',' fname]);
%         xlabel('local learn (n minus n-1)');
%         ylabel('ff dev (n+1 minus n-1)');
%
%         try
%             if use_targ_locallearn==1
%             x = DATTMP.(fname)(j).learnlocal;
%             else
%             x = DATTMP.(fname)(j).learnlocal;
%             end
%             if use_nminus2==1
%             y = DATTMP.(fname)(j).ffdev_first + DATTMP.(fname)(j).learnlocal;
%             else
%             y = DATTMP.(fname)(j).ffdev_first;
%             end
%             if takeabs==1
%                 x = abs(x);
%             end
%             plot(x, y, 'x');
%             lt_plot_makesquare_plot45line(gca, 'b');
%             lt_regress(y, x, 0, 0, 1, 1, 'r', 1);
%         catch err
%         end
%     end
% end
%
% linkaxes(hsplots, 'xy');
%
%
%
% %% ==================================== EXTRACT SLOPES (TRAIN, CATCH, BASE)
% % =========== 1) COLLECT
% plotON = 1;
% plotWNhits = 1;
% [SlopesAll, Learn_And_FFdevAll] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_7(DATTMP, ...
%     use_nminus2, plotON, [], plotWNhits);
%
%
% % ########################### PLOT (separately plots positive and negative
% % deviations)
% if use_nminus2==0
%     lt_figure;
%     lt_plot_text(0, 0.5, 'NOTE: using (n-2)-(n-1) for ffdev');
% end
%
% % ========== SUMMARY PLOT
% lt_figure; hold on;
% excludeIfMissCatch = 1;
%
% % =========
% lt_subplot(3,2,1); hold on;
% xlabel('NEG(train) - NEG(catch) | POS(train) - POS(catch)');
% ylabel('slope (trial n+2 vs trial n+1)');
%
% Y = SlopesAll(:, [1 2 4 5]);
% if excludeIfMissCatch==1
%     Y = Y(~any(isnan(Y)'), :);
% end
% x = 1:size(Y,2);
% plot(x,Y, '-ok');
% xlim([0 5]);
% lt_plot_zeroline;
%
% % ----------
% lt_subplot(3,2, 2); hold on;
% xlabel('NEG -- POS -- AVERAGE(POS,NEG)');
% ylabel('slope (subtract catch, pos = more learning than catch)');
% Ythis = [Y(:,1)-Y(:,2) Y(:,3)-Y(:,4)];
% Ythis(:,1) = -Ythis(:,1); % flip so that positive is more learning than catch.
% Ythis = [Ythis mean(Ythis, 2)]; % take average for 3rd column.
% x = 1:size(Ythis,2);
% plot(x, Ythis, '-ok');
% xlim([0 3]);
% lt_plot_zeroline;
%
% % ============ REPLOT ALL DATA, USING MEASURE WHERE UP IS MORE LEARNING
% Yadaptive = [1-SlopesAll(:,1:3) SlopesAll(:,4:6)];
%
% % --------------- 1) all slopes
% lt_subplot(3,2,3); hold on;
% excludeIfMissCatch = 1;
% xlabel('NEG(train) - NEG(catch) - NEG(base) | POS(train) - POS(catch) - POS(catch)');
% ylabel('slope (trial n+2 vs trial n+1) [adaptive dir]');
%
% Y = Yadaptive;
% if excludeIfMissCatch==1
%     Y = Y(~any(isnan(Y)'), :);
% end
% x = 1:size(Y,2);
% for j=1:size(Y,1)
%     plot(x, Y(j,:), '-ok');
% end
% xlim([0 7]);
%
% % ---------------- 2) MEAN DEVIATIONS
% lt_subplot(3,2,4); hold on;
% excludeIfMissCatch = 1;
% xlabel('NEG(train) - NEG(catch) - NEG(base) | POS(train) - POS(catch) - POS(catch)');
% ylabel('mean ffdev (trial n+2 minus n)');
%
% Y = Learn_And_FFdevAll;
% if excludeIfMissCatch==1
%     Y = Y(~any(isnan(SlopesAll)'), :);
% end
%
% functmp = @(x)(mean(x(:,2)));
% FFdev_mean = cell2mat(cellfun(functmp, Y, 'UniformOutput', 0));
%
% x = 1:size(FFdev_mean,2);
% plot(x, FFdev_mean, '-ok');
% xlim([0 7]);
% lt_plot_zeroline;
%
%
%
% %% =============== COLLECT SLOPES, SEPARATING BY TIME BIN OF DEVIATION
% tbinedges = [2]; % minutes, will fill in the edges
%
% % ========== COLLECT MULTIPLE TIME BINS AND PLOT
% plotON = 0;
% twind = [2 60];
% [SlopesAll, Learn_And_FFdevAll] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_7(DATTMP, ...
%     use_nminus2, plotON, twind);
%
%
% %% ============= CROSS CORRELATION OF DEVIATIONS
%
% maxwind = 50;
% fnames = fieldnames(DATTMP);
%
% CCAll = cell(length(DATTMP.Train), 3); % [syls x [train, catch, base]];
% for j=1:length(DATTMP.Train)
%     rowcount = 1;
%     for fname = fnames'
%         fname = fname{1};
%
%         if length(DATTMP.(fname))>=j
%         y1 = DATTMP.(fname)(j).learnlocal;
%         y2 = DATTMP.(fname)(j).ffdev_first;
%
%         if any(isnan(y1)) | any(isnan(y2))
%             assert(sum(isnan(y1))/length(y1) < 0.2, 'too much nans...');
%             assert(sum(isnan(y2))/length(y2) < 0.2, 'too much nans...');
%
%             indstmp = any(isnan([y1 y2])');
%             y1(indstmp) = [];
%             y2(indstmp) = [];
%         end
%
%         [cc, lags] = xcov(y1, y2, 20, 'coeff');
%
%
%         else
%             cc = [];
%             lags = [];
%         end
%
%         % ============
%         CCAll{j, rowcount} = cc;
%         rowcount = rowcount+1;
%     end
% end
%
% % ============== IMAPORTANT
% lt_figure;
% lt_plot_text(0, 0.5, 'NOTE: do not plot catch, since they are not temporally sequenced...');
%
%
% % ================ plot
% figcount=1;
% subplotrows=4;
% subplotcols=2;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
%
% onlykeepfull=1;
% if onlykeepfull==1
%    CCAll(any(cellfun(@isempty, CCAll)'), :) = [];
% end
%
%
% % ====== 1) overlay all xcovs
% pcollist = {'r', 'm', 'k'};
% for j=1:size(CCAll,1)
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots =[hsplots hsplot];
%
% % ---- train
%     indtmp = 1;
%     pcol = 'r';
%
%     cc = CCAll{j, indtmp};
%     plot(lags, cc, '-', 'Color', pcol);
%
% %     % ---- catch
% %     indtmp = 2;
% %     pcol = 'm';
% %
% %     cc = CCAll{j, indtmp};
% %     plot(lags, cc, '-', 'Color', pcol);
%
%     % ---- base
%     indtmp = 3;
%     pcol = 'k';
%
%     cc = CCAll{j, indtmp};
%     plot(lags, cc, '-', 'Color', pcol);
%
%     % ======
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
% end
%
% linkaxes(hsplots, 'xy');
%
% lt_figure; hold on;
% title('r=train; k=base');
%     CCAlltmp = cellfun(@transpose, CCAll, 'UniformOutput', 0);
% for i=[1 3]
%     lt_subplot(2,2,i); hold on;
%     xlabel('local laerning (trial n+1 minus n)');
%     ylabel('ffdev (usually trial n+2 minus trial n+1)');
%
%     pcol = pcollist{i};
%
%     ccall = cell2mat(CCAlltmp(:, i));
% %
% %     ccall = reshape(cell2mat(CCAll(:,i)), size(CCAll,1), []);
%
%     plot(lags, ccall, '-k');
%     ymean = mean(ccall,1);
%     ysem = lt_sem(ccall);
%     x = 1:length(ymean);
%
%     lt_plot(lags, ymean, {'Errors', ysem, 'Color', pcol});
%
%     % ======
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
% end
%
%     % ======== subtract  base from train
%     ccall_train = cell2mat(CCAlltmp(:, 1));
%     ccall_base = cell2mat(CCAlltmp(:, 3));
%
%     ccnorm = ccall_train - ccall_base;
%
%     lt_subplot(2,2,4); hold on;
%     xlabel('local laerning (trial n+1 minus n)');
%     ylabel('ffdev (usually trial n+2 minus trial n+1)');
%     title('train minus base');
%     pcol = 'b';
%
%     plot(lags, ccnorm, '-k');
%     ymean = mean(ccnorm,1);
%     ysem = lt_sem(ccnorm);
%     x = 1:length(ymean);
%
%     lt_plot(lags, ymean, {'Errors', ysem, 'Color', pcol});
%
%
%
%
%
%
%
%
%






















