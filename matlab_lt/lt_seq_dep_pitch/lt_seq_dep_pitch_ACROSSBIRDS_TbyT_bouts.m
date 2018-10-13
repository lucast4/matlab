function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_bouts(DATBYREND, TrialStruct, addtoibitime, ...
    userankcorr)
%% lt 10/6/18 - song bout structure (i.e. short interbout intervals vs. longer

% addtoibitime = 3; % minutes, to add on to each specific ibi.

%% hand code isis for each experiment

InterboutList = {...
    {'pu53wh88', 'SeqDepPitchShift3', 1.15}, ...
    {'pu53wh88', 'SeqDepPitchShift', 1}, ...
    {'pu11wh87', 'SeqDepPitchLMAN', 1.1}, ...
    {'pu11wh87', 'SeqDepPitchShift2', 0.85}, ...
    {'bk34bk68', 'SeqDepPitchLMAN3', 1}, ...
    {'wh25pk77', 'SeqDepPitchLMAN', 1.1}, ...
    {'wh6pk36', 'LMANlearn2', 1.75}, ...
    {'bu77wh13', 'LMANlearn1', 1.05}, ...
    {'or74bk35', 'LMANneural2', 1.6}, ...
    {'pu69wh78', 'RALMANlearn1', 1.15}, ...
    {'wh44wh39', 'RALMANlearn2', 1.8}, ...
    {'wh44wh39', 'RALMANlearn3', 1.5}, ...
    }; % in minutes


%% EXTRACT EXPERIMENTS
getsiglearn = 0;
birdstoget = 'notSDP';
syltype = 'all';
minrends = 10;
minbaserends = 0; % has at least this many baseline renditions.

[Inds_sylcounter, Inds_birdnum, Inds_exptnum, ...
    IsSame, IsTarg] = ...
    lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_2(DATBYREND, TrialStruct, ...
    syltype, getsiglearn, birdstoget, minbaserends, minrends);



%% ========= FIGURE OUT ISI DIVIDER FOR EACH EXPERIMENT

plotlog = 1;

nbirds = max(Inds_birdnum);
nexpts = max(Inds_exptnum);

figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

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
        
        % ========== figure out threshold time.
        ibi_this = [];
        for tmp = InterboutList
            if strcmp(tmp{1}{1}, bname) & strcmp(tmp{1}{2}, ename)
                ibi_this = tmp{1}{3} + addtoibitime;
            end
        end
        
        indsthis = DATBYREND.Birdnum==i & DATBYREND.Exptnum==ii;
        
        tdevs = cell2mat(DATBYREND.Time_dev(indsthis)).*24*60;
        
        % ---- 1) not log10
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([bname '-' ename]);
        xlabel('time devs (all syls combined)');
        lt_plot_histogram(tdevs, '', 1, 0, '', '', 'k');
        axis tight;
        line([ibi_this ibi_this], ylim, 'Color', 'r');
        
        % --- 2) log10
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([bname '-' ename]);
        lt_plot_histogram(log10(tdevs), '', 1, 0, '', '', 'k');
        axis tight;
        line(log10([ibi_this ibi_this]), ylim, 'Color', 'r');
        
        % ---- 3) trials (dots)
        indtmp = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==1);
        tvals = TrialStruct.birds(i).exptnum(ii).sylnum(indtmp).Tvals.*24.*60;
        ffvals = TrialStruct.birds(i).exptnum(ii).sylnum(indtmp).FFvals;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([bname '-' ename]);
        ylabel('targ');
        plot(tvals, ffvals, 'ok');
        axis tight;
        
        
    end
end


%% EXTRACT FOR EACH SYL 2 ARRAYS (ffbals for t shorter and longer than the ibi)

FFvalsPairsAll = {};

LongInter_Tdev = {};
LongInter_FFdev = {};
LongInter_NumShortPre = {};
% LongInter_NumShortPre = {};
for i=1:length(Inds_sylcounter)
    sylthis = Inds_sylcounter(i);
    
    % =========== EXTRACT DATA FOR THIS SYL
    indtmp = DATBYREND.Sylcounter==sylthis & DATBYREND.IsDurTrain==1 & DATBYREND.IsCatch==0; % rends
    bnum = unique(DATBYREND.Birdnum(indtmp));
    enum = unique(DATBYREND.Exptnum(indtmp));
    
    bname = TrialStruct.birds(bnum).birdname;
    ename = TrialStruct.birds(bnum).exptnum(enum).exptname;
    %         learndir = TrialStruct.birds(bnum).exptnum(enum).targlearndir;
    
    % --------- what is ibi?
    ibi_this = [];
    for tmp = InterboutList
        if strcmp(tmp{1}{1}, bname) & strcmp(tmp{1}{2}, ename)
            ibi_this = tmp{1}{3} + addtoibitime;
        end
    end
    
    
    % ----------- EXTRACT
    t = cell2mat(DATBYREND.Time_dev(indtmp))*24*60;
    ff = cell2mat(DATBYREND.FF_dev(indtmp));
    % split into 2 bins
    FFvals = {ff(t<ibi_this), ff(t>ibi_this)};
    
    % ========== STORE
    FFvalsPairsAll = [FFvalsPairsAll; FFvals];
    
    
    % ########################### EXTRACT BETWEEN EDGES OF BOUTSEQUENCES
    % FOR EACH LONGER DURATION INTERVAL, COUNT THE NUMBER OF SHORTER
    % INTERVALS THAT OCCURED BEFORE IT.
    t = cell2mat(DATBYREND.Time_dev(indtmp))*24*60;
    ff = cell2mat(DATBYREND.FF_dev(indtmp));
    
    inds_endofbouts = find(t>ibi_this);
    tdevall = [];
    ffdevall = [];
    numpreall = [];
    for j=2:length(inds_endofbouts)
        indtmp = inds_endofbouts(j);
        % -- collect dat
        tdevall = [tdevall; t(indtmp)];
        ffdevall = [ffdevall; ff(indtmp)];
        numpre = inds_endofbouts(j) - inds_endofbouts(j-1) -1; % number of short intervals between this long and previous long
        numpreall = [numpreall; numpre];
    end
LongInter_Tdev = [LongInter_Tdev; tdevall];
LongInter_FFdev = [LongInter_FFdev; ffdevall];
LongInter_NumShortPre = [LongInter_NumShortPre; numpreall];
    
end



%% ============= ONE PLOT FOR EACH EXPT (TAKES ALL DATAPOINTS ACROSS SYLS FOR A GIVEN SYL TYPE)

FFtarg = [];
FFsame = [];
FFdiff = [];
onlyifhaveallsyltypes = 1; % NOT YET DONE!!!!
for i = 1:nbirds
    for ii=1:nexpts
        
        % ============ SYLS FOR THIS EXPT
        indthis = Inds_birdnum==i & Inds_exptnum==ii;
        if ~any(indthis)
            continue
        end
        
        % =================== TARG
        indthis = Inds_birdnum==i & Inds_exptnum==ii & IsSame==1 & IsTarg==1;
        ffmean = cellfun(@mean, FFvalsPairsAll(indthis, :));
        ffsem = cellfun(@lt_sem, FFvalsPairsAll(indthis, :));
        
        FFtarg = [FFtarg; ffmean];
        
        % =================== SAME
        indthis = Inds_birdnum==i & Inds_exptnum==ii & IsSame==1 & IsTarg==0;
        tmp = FFvalsPairsAll(indthis, :);
        tmp = {cell2mat(tmp(:,1)), cell2mat(tmp(:,2))}; % combine all renditions across syls
        ffmean = cellfun(@nanmean, tmp);
        ffsem = cellfun(@lt_sem, tmp);
        
        FFsame = [FFsame; ffmean];
        
        % =================== DIFF
        indthis = Inds_birdnum==i & Inds_exptnum==ii & IsSame==0 & IsTarg==0;
        tmp = FFvalsPairsAll(indthis, :);
        tmp = {cell2mat(tmp(:,1)), cell2mat(tmp(:,2))}; % combine all renditions across syls
        ffmean = cellfun(@nanmean, tmp);
        ffsem = cellfun(@lt_sem, tmp);
        
        FFdiff = [FFdiff; ffmean];
        
    end
end

% ======== PLOT
lt_figure; hold on;

lt_subplot(3,2,1); hold on;
pcol = 'k';
title('TARG');
xlabel('short - long (inter bout interval)');
ylabel('ff change');
Y = FFtarg;
x = [1 2];
plot(x, Y', '-o', 'Color', pcol);
xlim([0 3]);
lt_plot_zeroline;

lt_subplot(3,2,2); hold on;
pcol = 'b';
title('SAME');
xlabel('short - long (inter bout interval)');
ylabel('ff change');
Y = FFsame;
x = [1 2];
plot(x, Y', '-o', 'Color', pcol);
xlim([0 3]);
lt_plot_zeroline;

lt_subplot(3,2,3); hold on;
pcol = 'r';
title('DIFF');
xlabel('short - long (inter bout interval)');
ylabel('ff change');
Y = FFdiff;
x = [1 2];
plot(x, Y', '-o', 'Color', pcol);
xlim([0 3]);
lt_plot_zeroline;

lt_subplot(3,2,4); hold on;
title('same, subtract targ');
x = [1 2];
Y = FFsame - FFtarg;
plot(x, Y', '-o', 'Color', 'b');
xlim([0 3]);
lt_plot_zeroline;

lt_subplot(3,2,5); hold on;
title('same, subtract diff');
x = [1 2];
Y = FFsame - FFdiff;
plot(x, Y', '-o', 'Color', 'b');
xlim([0 3]);
lt_plot_zeroline;


lt_subplot(3,2,6); hold on;
title('targ, subtract diff');
x = [1 2];
Y = FFtarg - FFdiff;
plot(x, Y', '-o', 'Color', 'b');
xlim([0 3]);
lt_plot_zeroline;


%% ============= PREDICT CHANGE AT LONGER GAP BASED ON STUFF OCCURING PREVIOUS TO THIS GAP?
figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

slopes_targ = [];
slopes_same = [];
slopes_same_expt = [];
% ======================== PLOT ALL SYLS
for i = 1:nbirds
    for ii=1:nexpts
        
        % ============ SYLS FOR THIS EXPT
        indthis = Inds_birdnum==i & Inds_exptnum==ii;
        if ~any(indthis)
            continue
        end
        
        % =================== TARG
        indthis = Inds_birdnum==i & Inds_exptnum==ii & IsSame==1 & IsTarg==1;
        pcol = 'k';
        numpreall = LongInter_NumShortPre{indthis};
        ff = LongInter_FFdev{indthis};
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([bname '-' ename '[TARG]']);

        [b] = lt_regress(ff, numpreall, 1, 0, 1, 1, pcol);
%         lt_regress(ff, log10(numpreall+0.01), 1, 0, 1, 1, pcol);
    if userankcorr==1
                slopes_targ = [slopes_targ; corr(ff, numpreall, 'type', 'Spearman')];

    else
        slopes_targ = [slopes_targ; b(2)];
    end
    
        % ================== SAME
        indthisall = find(Inds_birdnum==i & Inds_exptnum==ii & IsSame==1 & IsTarg==0);
        slopestmp = [];
        for indthis = indthisall
            pcol = 'b';
            numpreall = LongInter_NumShortPre{indthis};
            ff = LongInter_FFdev{indthis};
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([bname '-' ename '[SAME]']);
            
            [b] = lt_regress(ff, numpreall, 1, 0, 1, 1, pcol);
            %         lt_regress(ff, log10(numpreall+0.01), 1, 0, 1, 1, pcol);
            
            if userankcorr==0
            slopes_same= [slopes_same; b(2)];
            slopestmp = [slopestmp; b(2)];
            elseif userankcorr==1                
            slopes_same= [slopes_same; corr(ff, numpreall, 'type', 'Spearman')];
            slopestmp = [slopestmp; corr(ff, numpreall, 'type', 'Spearman')];
            end
        end
        slopes_same_expt = [slopes_same_expt; nanmean(slopestmp)];
    end
end


lt_figure; hold on;

% +====== 
lt_subplot(3,2,1); hold on;
title('targ(k), same(b) [all syls]');

%  targ
pcol = 'k';
y = slopes_targ;
x = '';
lt_plot_histogram(y, x, 1, 1, '', 1, pcol);

%  same
pcol = 'b';
y = slopes_same;
x = '';
lt_plot_histogram(y, x, 1, 1, '', 1, pcol);

% ---- compare
lt_subplot(3,2,2); hold on
title('each expt one val, mean of slopes');
xlabel('targ -- same');
ylabel('slope (ffdev vs. numpreshort');
x = [1 2];
y = [slopes_targ slopes_same_expt];
plot(x, y', '-ok');
xlim([0 3]);

