%% ======= [PLOT] FOR EACH SYL, PLOT RAW AND DISTRIBUTIONS OF DENSITIES

% ============== DURINGIN TRAINNIG?
dotrain = 1;


% ======= only take data up to first nontarg after the referenc?
singleRendOnly=0;
% if 0, then takes all rends

% ====== how to compute density?
densitymethod = 'refperiod';
% refperiod: will use reference period
% entirebin: will use time from rendtion to maxtime_fromref (see below)
% beforefirstnontarg = will get number that is >= ref time and < time of first target post-reference.

if singleRendOnly==1
    % needs to be this to make sense ...
    densitymethod = 'beforefirstnontarg';
end

% ============ method for decideing hi and lo density trials
cutoffmethod = 'medianslice';
% medianslice: for each nontarg bin, finds median for targ
% medianoverall: overall median of targ

% ====== what time period to plot (locked to ref period) to summarize?
mintime_fromref = 5;
maxtime_fromref = 30; % minutes

% ----------- min num trials in this bin...
minrends_inbin = 20;


% ====
% singleRendOnly=0;
% densitymethod = 'refperiod';
% cutoffmethod = 'medianslice';
% mintime_fromref = 5;
% maxtime_fromref = 30; % minutes
% minrends_inbin = 20;

% ===========
nrendstoplot = 20;

if plotEachSylRaw==2
    % larger subplots
    figcount=1;
    subplotrows=2;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
else
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
end

Xall_LoDens = {};
Xall_HiDens = {};

Yall_LoDens = {};
Yall_HiDens = {};

Ybin_LoDens = [];
Ybin_HiDens = [];

nrendsinbinALL=[];

Nbybin_HiDens.targ = {};
Nbybin_HiDens.nontarg = {};
Nbybin_LoDens.targ = {};
Nbybin_LoDens.nontarg = {};
SigLearnNontarg = [];

for i=1:Numbirds
    Numexpt = length(TrialStruct.birds(i).exptnum);
    birdname = TrialStruct.birds(i).birdname;
    
    for ii=1:Numexpt
        exptname = TrialStruct.birds(i).exptnum(ii).exptname;
        
        % ---------- SKIP IF NO DATA
        if isempty(TrialStruct.birds(i).exptnum(ii).sylnum)
            disp(['[no DATA] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
            continue
        end
        
        % --------- ignore if lMAN?
        if ignoreLMANexpt==1
            if isfield(TrialStruct.birds(i).exptnum(ii), 'LMANinactivated')
                isLMAN = TrialStruct.birds(i).exptnum(ii).LMANinactivated;
                if isLMAN==1
                    disp(['[is LMAN] skipping ' TrialStruct.birds(i).exptnum(ii).exptname]);
                    continue
                end
            end
        end
        
        % ===================== INFOR ABOUT TARGET
        indtarg = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==1);
        tvals_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).Tvals;
        ffvals_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indtarg).FFvals;
        
        % ############################################ SAME
        issame=1;
        sylinds_same = find([TrialStruct.birds(i).exptnum(ii).sylnum.INFO_similar]==issame ...
            & [TrialStruct.birds(i).exptnum(ii).sylnum.INFO_istarget]==0);
        
        % ===================== GO THRU ALL NONTARG SYLS. FOR EACH REND COLLECT
        % DATA
        for j=1:length(sylinds_same)
            
            indthis = sylinds_same(j);
            sylthis = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).syl;
            
            
            % ================== NOTE DOWN AMOUNT OF LEARNING FOR THIS SYL
            
            WNon = TrialStruct.birds(i).exptnum(ii).WNontime;
            ttmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
            fftmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals;
            indstrain = ttmp>WNon;
            [b,bint] = lt_regress(fftmp(indstrain), ttmp(indstrain), 0);
            isSigLearn = sign(bint(2,1))==sign(bint(2,2)); % if CI of slope does not overlap 0;
            
            % ============= EXTRACT FOR THIS SYL
            if strcmp(densitymethod, 'refperiod')
                % ---------------------- 1) NTARGS VS. NNONTARGS IN REF WINDOW
                ntarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NTargRendsInRef;
                nnontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NNonTargRendsInRef;
                
            elseif strcmp(densitymethod, 'entirebin')
                %                 maxtime_day = maxtime_fromref/(60*24);
                mintime_day = 0/(60*24);
                maxtime_day = maxtime_fromref/(60*24);
                
                % ---------------- COUNT NUMBER OF targ and nontarg between
                % ref syl and end of bin [NOTE: end is potentialy not exact
                % end of last bin...]
                timedevs_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
                timedevs_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs;
                
                nrends = length(timedevs_nontarg);
                
                ntarg = nan(nrends,1);
                nnontarg = nan(nrends,1);
                for rr =1:nrends
                    if isempty(timedevs_nontarg{rr})
                        continue
                    end
                    %                     nnontarg(rr) = sum(timedevs_nontarg{rr}>=0 & timedevs_nontarg{rr}<maxtime_day);
                    %                     ntarg(rr) = sum(timedevs_targ{rr}>=0 & timedevs_targ{rr}<maxtime_day);
                    nnontarg(rr) = sum(timedevs_nontarg{rr}>=mintime_day & ...
                        timedevs_nontarg{rr}<maxtime_day);
                    ntarg(rr) = sum(timedevs_targ{rr}>=mintime_day & ...
                        timedevs_targ{rr}<maxtime_day);
                end
                
                % ===================== COMPARE IT TO REFPERIOD
                % ---------------------- 1) NTARGS VS. NNONTARGS IN REF WINDOW
                ntarg_REF = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NTargRendsInRef;
                nnontarg_REF = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_NNonTargRendsInRef;
                
            elseif strcmp(densitymethod, 'beforefirstnontarg')
                
                % ------------------ COUNT RENDITIONS FROM REF (INCLUSIVE)
                % UP TO FIRST NONTARG POST REF
                timedevs_nontarg = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
                timedevs_targ = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs;
                
                nrends = length(timedevs_nontarg);
                ntarg = nan(nrends,1);
                nnontarg = nan(nrends,1);
                for rr =1:nrends
                    if isempty(timedevs_nontarg{rr})
                        continue
                    end
                    
                    % ---- what is time of first nontarg post ref>?
                    timeoffirstNONTARG = timedevs_nontarg{rr}(find(timedevs_nontarg{rr}>0, 1, 'first'));
                    if isempty(timeoffirstNONTARG)
                        timeoffirstNONTARG=0;
                    end
                    
                    nnontarg(rr) = sum(timedevs_nontarg{rr}>=0 & timedevs_nontarg{rr}<timeoffirstNONTARG);
                    ntarg(rr) = sum(timedevs_targ{rr}>=0 & timedevs_targ{rr}<timeoffirstNONTARG);
                end
            end
            assert(all(isnan(ntarg)==isnan(nnontarg)));
            
            
            % ---------------------- 1B) GET CUTOFF LINE
            if strcmp(cutoffmethod, 'medianslice')
                %  for each value of nontarg, get median value of targ, and use
                %  that as threshold
                cutoffs_y_targ = grpstats(ntarg, nnontarg, {'median'});
                cutoffs_x_nontarg = unique(nnontarg);
                cutoffs_x_nontarg = unique(nnontarg(~isnan(nnontarg)));
                
                [~, indstmp] = ismember(nnontarg, cutoffs_x_nontarg);
                indstmp(indstmp==0) = 1; % temporarily make 1, since will convert to nan in a bit
                % ---- for each rend, determine whether it was high or low
                % density (in terms of targ syls)
                ishighdensity = ntarg>cutoffs_y_targ(indstmp);
                ishighdensity = single(ishighdensity);
                ishighdensity(isnan(ntarg)) = nan;
                
            elseif strcmp(cutoffmethod, 'medianoverall')
                
                cutoffthis = nanmedian(ntarg);
                
                ishighdensity = ntarg>cutoffthis;
                ishighdensity = single(ishighdensity);
                ishighdensity(isnan(ntarg)) = nan;
            end
            TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_IsHighDensTarg = ishighdensity;
            
            
            
            
            % ============================= PLOT (SHOWING SEPARTE HIGH AND
            % LOW DENSITY RENDS)
            if plotEachSylRaw==1 | plotEachSylRaw==2
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title('density of targ syl (r=high, b=lo)');
                ylabel('N targs in ref');
                xlabel('N nontags in ref');
                title([birdname '-' exptname '-' sylthis]);
                
                % 1) add some jitter
                nnontarg = nnontarg + 0.4*rand(size(nnontarg)) -0.2;
                ntarg = ntarg + 0.4*rand(size(nnontarg)) -0.2;
                
                % 2) separte colors for high and low density
                % - high
                indstmp = ishighdensity==1;
                pcol = 'r';
                plot(nnontarg(indstmp), ntarg(indstmp), 'x', 'Color', pcol);
                xmean = mean(nnontarg(indstmp)); % overlay mean and sd
                xstd = std(nnontarg(indstmp));
                ymean = mean(ntarg(indstmp));
                ystd = std(ntarg(indstmp));
                lt_plot(xmean, ymean, {'Errors', ystd, 'Color', pcol});
                line(xmean+[-xstd xstd], [ymean ymean], 'Color', pcol);
                
                % - lo
                indstmp = ishighdensity==0;
                pcol = 'b';
                plot(nnontarg(indstmp), ntarg(indstmp), 'x', 'Color', pcol);
                % ---- format
                lt_plot_makesquare_plot45line(gca, 'k', -1);
                xmean = mean(nnontarg(indstmp)); % overlay mean and sd
                xstd = std(nnontarg(indstmp));
                ymean = mean(ntarg(indstmp));
                ystd = std(ntarg(indstmp));
                lt_plot(xmean, ymean, {'Errors', ystd, 'Color', pcol});
                line(xmean+[-xstd xstd], [ymean ymean], 'Color', pcol);
                
                % ======================= COMPARE THIS METHOD TO REFERENCE
                % METHOD
                %                 strcmp(densitymethod, 'entirebin')
            end
            
            
            % ========================= PLOT RANDOM INDIVIDUAL TRIALS,
            % REPRESENTING NRENDS IN TARG AND NONTARG, AND OTHER TRIALS ...
            
            
            
            % ---------------------- 2) Overlay entire learning trajectory
            % with renditions from which datapoints were extracted
            % IN PROGRESS
            if plotEachSylRaw==2
                
                % ------------ targ context
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                
                % --------------- note down times of all rends that were
                % used
                hiDens_tmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_IsHighDensTarg;
                timedevNT_tmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev;
                tvals_tmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).Tvals;
                ffvals_tmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).FFvals;
                
                plot(tvals_targ, ffvals_targ, '.k');
                plot(ttmp, fftmp, '.b'); % ALL RENDS
                plot(ttmp(hiDens_tmp==1), fftmp(hiDens_tmp==1), '.r'); % HI DENS
                plot(ttmp(hiDens_tmp==0), fftmp(hiDens_tmp==0), '.c'); % HI DENS
                
                edge_time = xedges(find(xedges<maxtime_fromref, 1, 'last')); % bin centers might be lower than this (i.e. thsi is the right edge of the last bin to collect data from)
                
                for kk=1:length(tvals_tmp)
                    
                    if isempty(timedevNT_tmp{kk})
                        continue
                    end
                    
                    tthistmp = tvals_tmp(kk);
                    ffthistmp = ffvals_tmp(kk);
                    flankright = tthistmp+flanktime_targ/(60*24);
                    flankleft = tthistmp-flanktime_targ/(60*24);
                    
                    % ---- WINDOW for collecting all values precedinga dn
                    % post
                    datedges = tthistmp + twind_plot./(24);
                    line([datedges], [ffthistmp ffthistmp], ...
                        'Color', [0.6 0.6 0.6], 'LineStyle', '--');
                    
                    
                    % ----- window for counting number of rends [in plot
                    % that will make in next section]
                    edgeright = tthistmp + edge_time/(60*24);
                    line([tthistmp edgeright], [ffthistmp ffthistmp], ...
                        'Color', 'm', 'LineStyle', ':')
                    
                    
                    % ----- reference window
                    line([flankleft flankright], [ffthistmp ffthistmp], ...
                        'Color', 'c')
                    
                    
                    % ---- plot number of left and right rends collected
                    nrighttmp = sum(timedevNT_tmp{kk} > 0);
                    nmidtmp = sum(timedevNT_tmp{kk} == 0);
                    nlefttmp = sum(timedevNT_tmp{kk} < 0);
                    lt_plot_text(edgeright, ffthistmp, ...
                        ['n=' num2str(nlefttmp) '-' num2str(nmidtmp) '-' num2str(nrighttmp)], [0.6 0.6 0.6], 5);
                end
            end
            
            
            % ================================= 3) plot syl locked deviation
            % *** DUR TRAIN
            % -------------------- PARAMS:
            
            % --------------------- RUNS
            if plotEachSylRaw==1 | plotEachSylRaw==2
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                %             hsplots = [hsplots hsplot];
                xlabel('time dev from ref point (min)');
                title('by targ density (r=hi; b = lo)');
            end
            
            densitylevels = [0 1];
            for kk = densitylevels
                % --- high density
                hidensity = kk;
                
                % ---------------- COLLECT AND PLOT
                indstmp = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).IsDurTrain==dotrain ...
                    & TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_IsHighDensTarg==hidensity;
                if hidensity==1
                    pcol = 'r';
                else
                    pcol = 'b';
                end
                
                timedev = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDev(indstmp);
                ffdev = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_FFDev(indstmp);
                timedev_TARG = TrialStruct.birds(i).exptnum(ii).sylnum(indthis).RefLocked_TimeDevTargs(indstmp);
                
                
                % ################# OVERLAY BINNED MEANS
                % - 1) COLLECT ALL DATA INTO VECTOR
                nrends = length(timedev);
                timedev_all = [];
                ffdev_all = [];
                timedev_targ_all = [];
                for k=1:nrends
                    assert(~any(isnan(timedev{k})), 'then cant count rends using this ... some nan...');
                    
                    if singleRendOnly==1
                        % then only keeps timepoints up to first rendition
                        % of nontarget after reference point
                        timeoffirstNONTARG = timedev{k}(find(timedev{k}>0, 1, 'first'));
                        if isempty(timeoffirstNONTARG)
                            % then means there is no post-ref trial...
                            timeoffirstNONTARG = 1;
                        end
                        % ------ NONTARG (pare down) - all rends AT OR
                        % BEFORE
                        indstmptmp = timedev{k}<=timeoffirstNONTARG;
                        timedev{k} = timedev{k}(indstmptmp);
                        ffdev{k} = ffdev{k}(indstmptmp);
                        % -------- TARG (pare down) - all rends BEFORE
                        indstmptmp = timedev_TARG{k}<timeoffirstNONTARG;
                        timedev_TARG{k} = timedev_TARG{k}(indstmptmp);
                    end
                    timedev_all = [timedev_all; timedev{k}*(24*60)];
                    ffdev_all = [ffdev_all; ffdev{k}];
                    timedev_targ_all = [timedev_targ_all; timedev_TARG{k}*(24*60)];
                end
                
                
                if plotEachSylRaw==1 | plotEachSylRaw==2
                    % OLD METHOD, maintains within-trial correlations
                    %                     if length(timedev)>nrendstoplot
                    %                         rendstoplot = randperm(length(timedev), nrendstoplot);
                    %                     else
                    %                         rendstoplot = 1:length(timedev);
                    %                     end
                    %                     for k=rendstoplot
                    %                         t = timedev{k}*(24*60);
                    %                         ff = ffdev{k};
                    %                         plot(t, ff, 'x', 'Color', pcol);
                    %                     end
                    % NEW MOETHOD -required if using first syl only.
                    if length(timedev_all)>200
                        indstmptmp = randperm(length(timedev_all),200);
                    else
                        indstmptmp =1:length(timedev_all);
                    end
                    t = timedev_all(indstmptmp);
                    ff = ffdev_all(indstmptmp);
                    plot(t, ff, 'x', 'Color', pcol);
                    
                    
                end
                
                % ---------------- slide into bins
                xbins = discretize(timedev_all, xedges);
                
                if usemedian==1
                    [ymean, ysem] = grpstats(ffdev_all, xbins, {'median', 'sem'});
                else
                    [ymean, ysem] = grpstats(ffdev_all, xbins, {'mean', 'sem'});
                end
                X = xcenters(unique(xbins));
                
                
                % =================================== Get sample size per bin
                if ~isempty(ffdev_all)
                    % ------- NONTARGET
                    tmp = tabulate(xbins);
                    xtmp = tmp(tmp(:,2)~=0,1); assert(all(xtmp==unique(xbins)), 'not aligned...');
                    N_bybin_nontarg = tmp(tmp(:,2)~=0,2);
                    N_bybin_nontarg = N_bybin_nontarg./nrends; % to get sample size per rendition
                    
                    % -------- TARGET
                    xbins_TARG = discretize(timedev_targ_all, xedges);
                    X_TARG = xcenters(unique(xbins_TARG));
                    % -
                    tmp = tabulate(xbins_TARG);
                    xtmp = tmp(tmp(:,2)~=0,1); assert(all(xtmp==unique(xbins_TARG)), 'not aligned...');
                    N_bybin_TARG = tmp(tmp(:,2)~=0,2);
                    N_bybin_TARG = N_bybin_TARG./nrends;
                    
                else
                    N_bybin_nontarg = nan;
                    X = nan;
                    
                    N_bybin_TARG = nan;
                    X_TARG = nan;
                end
                
                
                % ------------------- PLOT BINNED DATA
                if plotEachSylRaw==1 | plotEachSylRaw==2
                    lt_plot(X, ymean, {'Errors', ysem, 'MarkerFaceColor', pcol, 'Color', 'k'});
                    xlim([-45 45]);
                    lt_plot_zeroline;
                    lt_plot_zeroline_vert;
                end
                
                
                % ========================== COLLECT
                
                % -------------- figure out number of trials in desired bin
                binsthis = find(xcenters>mintime_fromref & xcenters<maxtime_fromref);
                indsinbin = ismember(xbins, binsthis); % renditions that are in desired windwo
                
                % -------- only continue if enough renditions in bin
                nrendsinbin = sum(indsinbin);
                nrendsinbinALL = [nrendsinbinALL nrendsinbin];
                
                
                if nrendsinbin<minrends_inbin
                    lt_plot_annotation(1, ['insuff rends in bin,n=' num2str(nrendsinbin)], pcol);
                    ytosave = nan;
                    xtosave = nan;
                    ybinnedtosave = nan;
                else
                    ytosave = ymean;
                    xtosave = X;
                    ybinnedtosave = mean(ffdev_all(indsinbin)); % all renditions in desired bins
                end
                
                % ================== CONFIRM THAT DISTRIBUTIONS OF NONTARGS
                % ARE SIMILAR
                if hidensity==1
                    Xall_HiDens = [Xall_HiDens; xtosave];
                    Yall_HiDens = [Yall_HiDens; ytosave];
                    Ybin_HiDens = [Ybin_HiDens; ybinnedtosave];
                    SigLearnNontarg = [SigLearnNontarg; isSigLearn];
                    
                    Nbybin_HiDens.targ = [Nbybin_HiDens.targ; [X_TARG; N_bybin_TARG']];
                    Nbybin_HiDens.nontarg = [Nbybin_HiDens.nontarg; [X; N_bybin_nontarg']];
                else
                    Xall_LoDens = [Xall_LoDens; xtosave];
                    Yall_LoDens = [Yall_LoDens; ytosave];
                    Ybin_LoDens = [Ybin_LoDens; ybinnedtosave];
                    
                    Nbybin_LoDens.targ = [Nbybin_LoDens.targ; [X_TARG; N_bybin_TARG']];
                    Nbybin_LoDens.nontarg = [Nbybin_LoDens.nontarg; [X; N_bybin_nontarg']];
                end
                
                
            end
            
            
            % ============================== PLOT SAMPEL SIZES PER BIN
            if plotEachSylRaw==1 | plotEachSylRaw==2
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                %             hsplots = [hsplots hsplot];
                xlabel('time dev from ref point (min)');
                title('by targ density (r=hi; b = lo)');
                ylabel('N (solid=targ; dash=ntarg');
                
                % -- nontarg, hi dens
                pcol = 'r';
                dattmp = Nbybin_HiDens.nontarg{end};
                plot(dattmp(1,:), dattmp(2,:), '--o', 'Color', pcol);
                
                % -- nontarg, lo dens
                pcol = 'b';
                dattmp = Nbybin_LoDens.nontarg{end};
                plot(dattmp(1,:), dattmp(2,:), '--o', 'Color', pcol);
                
                % -- targ, hi dens
                pcol = 'r';
                dattmp = Nbybin_HiDens.targ{end};
                plot(dattmp(1,:), dattmp(2,:), '-o', 'Color', pcol);
                
                % -- targ, lo dens
                pcol = 'b';
                dattmp = Nbybin_LoDens.targ{end};
                plot(dattmp(1,:), dattmp(2,:), '-o', 'Color', pcol);
                
                % ---
                xlim([-60 60]);
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
            end
            
            
            if plotEachSylRaw==2
                pause
                close all;
            end
            
        end
    end
end


% =========== histogram of rends in analysis bin
if (0)
    lt_figure; hold on;
    % lt_plot_histogram(nrendsinbinALL);
    lt_plot_cdf(nrendsinbinALL);
end




%% ============ [PLOT] SUMMARY - ff dev for cases with high and low density targ syl
removeIfNotHigherTargDensity = 1; % if 0 still removes for the final plot.
numcases = length(Xall_LoDens);

% ================ TO COUNT NUMBER RENDS FROM REF POINT TO END OF BIN THAT
% IS DIRECTLY PREDECING THE LAST BIN
edge_time = xedges(find(xedges<maxtime_fromref, 1, 'last')); % bin centers might be lower than this (i.e. thsi is the right edge of the last bin to collect data from)


% ======================== COMPARE SAMPLE SIZES
lt_figure; hold on;
NtrialAll = [];
SigLearnAll = [];
for i=1:numcases
    
    if all(isnan(Xall_LoDens{i})) | all(isnan(Xall_HiDens{i}))
        disp(i);
        continue
    end
    
    NTHIS = [];
    % ---- targ, HI
    fname = 'targ';
    x = Nbybin_HiDens.(fname){i}(1,:);
    y = Nbybin_HiDens.(fname){i}(2,:);
    
    indsthis = x>=0 & x<edge_time;
    Ntrial = sum(y(indsthis));
    
    NTHIS = [NTHIS Ntrial];
    
    % ---- targ, LO
    fname = 'targ';
    x = Nbybin_LoDens.(fname){i}(1,:);
    y = Nbybin_LoDens.(fname){i}(2,:);
    
    indsthis = x>=0 & x<edge_time;
    Ntrial = sum(y(indsthis));
    
    NTHIS = [NTHIS Ntrial];
    
    % ---- nontarg, HI
    fname = 'nontarg';
    x = Nbybin_HiDens.(fname){i}(1,:);
    y = Nbybin_HiDens.(fname){i}(2,:);
    
    indsthis = x>=0 & x<edge_time;
    Ntrial = sum(y(indsthis));
    
    NTHIS = [NTHIS Ntrial];
    
    % ---- nontarg, LO
    fname = 'nontarg';
    x = Nbybin_LoDens.(fname){i}(1,:);
    y = Nbybin_LoDens.(fname){i}(2,:);
    
    indsthis = x>=0 & x<edge_time;
    Ntrial = sum(y(indsthis));
    
    NTHIS = [NTHIS Ntrial];
    
    % ================== COLLECT
    NtrialAll = [NtrialAll; NTHIS];
    SigLearnAll = [SigLearnAll; SigLearnNontarg(i)];
end

% ============ what are bad inds?
indsbad = find(NtrialAll(:,1) < 1.05*NtrialAll(:,2));

if (0)
    tmp1 = NtrialAll(:,1)./NtrialAll(:,3)
    tmp2 = NtrialAll(:,2)./NtrialAll(:,4)
    indsbad = find(tmp1<tmp2);
end

if removeIfNotHigherTargDensity==1
    NtrialAll(indsbad,:) = [];
    SigLearnAll(indsbad) = [];
end



% ================ PLOT
NtrialAll = log2(NtrialAll);
if (1)
    lt_subplot(2,2,1); hold on;
    Ncell = {NtrialAll(:,1), NtrialAll(:,2), NtrialAll(:,3), NtrialAll(:,4)};
    lt_plot_MultDist(Ncell, [1 2 3 4]);
    
    xlabel('targ/hi -- targ/lo -- nontarg/hi -- nontarg/lo');
    ylabel('sample sizes');
end

lt_subplot(2,2,2); hold on;
boxplot(NtrialAll);

lt_subplot(2,2,3); hold on;
title('targ');
xlabel('lo density');
ylabel('hi density');
plot(NtrialAll(:,2), NtrialAll(:,1), 'ok')
lt_plot_makesquare_plot45line(gca, 'b', -1);

lt_subplot(2,2,4); hold on;
title('NONtarg');
xlabel('lo density');
ylabel('hi density');
plot(NtrialAll(:,4), NtrialAll(:,3), 'ok')
lt_plot_makesquare_plot45line(gca, 'b', -1);

p1 = signrank(NtrialAll(:,1), NtrialAll(:,2));
lt_plot_text(1, max(NtrialAll(:)), ['p=' num2str(p1)], 'r');
p2 = signrank(NtrialAll(:,3), NtrialAll(:,4));
lt_plot_text(3, max(NtrialAll(:)), ['p=' num2str(p2)], 'r');



% ####################################################### PLOT SUMMARY
lt_figure;

xcollect.lo = [];
xcollect.hi = [];

ycollect.lo = [];
ycollect.hi = [];

counter.lo = [];
counter.hi = [];

Y_lohi = [];
count = 1;

for i=1:numcases
    
    if all(isnan(Xall_LoDens{i})) | all(isnan(Xall_HiDens{i}))
        disp(i);
        continue
    end
    
    % ########################## 1) PLOT BINNED RESPOSNES, LOCKED TO REFERENCE
    lt_subplot(3,1,1); hold on;
    title('all syllables');
    xlabel('time rel nontarg ref (min)');
    ylabel('nontarg ff dev from ref (hz)');
    
    % ===================================== low density
    pcol = 'b';
    x = Xall_LoDens{i};
    y = Yall_LoDens{i};
    plot(x,y, '-x', 'Color', pcol);
    
    % ----- COLLECT, FOR GLOBAL MEAN
    xcollect.lo = [xcollect.lo; x'];
    ycollect.lo = [ycollect.lo; y];
    counter.lo = [counter.lo; count*ones(size(x'))];
    
    
    % ===================================== hi density
    pcol = 'r';
    x = Xall_HiDens{i};
    y = Yall_HiDens{i};
    plot(x,y, '-x', 'Color', pcol);
    
    
    % ----- COLLECT, FOR GLOBAL MEAN
    xcollect.hi = [xcollect.hi; x'];
    ycollect.hi = [ycollect.hi; y];
    counter.hi = [counter.hi; count*ones(size(x'))];
    
    % -----
    
    % ------ formats
    xlim([-25 25]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    count = count+1;
    % ================ COLLECT VALUES FOR BINNED MEAN
    Y_lohi = [Y_lohi; ...
        Ybin_LoDens(i) Ybin_HiDens(i)];
end

if removeIfNotHigherTargDensity==1
    Y_lohi(indsbad,:) = [];
    
    indstoremove = ismember(counter.lo, indsbad);
    xcollect.lo(indstoremove) = [];
    ycollect.lo(indstoremove) = [];
    
    indstoremove = ismember(counter.hi, indsbad);
    xcollect.hi(indstoremove) = [];
    ycollect.hi(indstoremove) = [];
    
end


% ============================ OVERLAY MEAN ACROSS SYLS
lt_subplot(3,1,2); hold on;
title('all syllables');
xlabel('time rel nontarg ref (min)');
ylabel('nontarg ff dev from ref (hz)');

% --------------- 1) low
pcol ='b';
fname = 'lo';

x = xcollect.(fname);
y = ycollect.(fname);

[ymean, ysem]=grpstats(y, x, {'mean', 'sem'});
xx = unique(x);
lt_plot(xx, ymean, {'Errors', ysem, 'Color', pcol, 'MarkerFaceColor', 'k'});

% --------------- 1) hi
pcol ='r';
fname = 'hi';

x = xcollect.(fname);
y = ycollect.(fname);

[ymean, ysem]=grpstats(y, x, {'mean', 'sem'});
xx = unique(x);
lt_plot(xx, ymean, {'Errors', ysem, 'Color', pcol, 'MarkerFaceColor', 'k'});

% ------ formats
xlim([-25 25]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ######################## 2) FOR EACH CASE, COLLECT PAIRED DATA
% (hi density, lo density ...)
lt_subplot(3,2,5); hold on;
xlabel('Low targ dens -- Hi Targ dens');
ylabel('ff dev, nontarg, in bin');

plot([1 2], Y_lohi', '-', 'Color', [0.7 0.7 0.7]);
lt_plot([1.1 2.1], mean(Y_lohi), {'Errors', lt_sem(Y_lohi)});

% --- annotate cases that do not ahve difference in targ sampel size for
% high density vs low density
if removeIfNotHigherTargDensity==0
    for j=indsbad'
        plot(2.2, Y_lohi(j,2), 'rx');
        plot([1 2], Y_lohi(j,:), '-r');
    end
end

[~, p] = ttest(Y_lohi(:,1), Y_lohi(:,2));
lt_plot_pvalue(p, 'ttest', 1);

lt_plot_zeroline;


% ################################ ANNOTATE BY SIGNIFICANT LEARNING?
% ============ SIGNIFICANT
ytoplot= Y_lohi(logical(SigLearnAll),:);
plot([3 4], ytoplot', '-b');
lt_plot([3 4]+0.2, mean(ytoplot), {'Errors', lt_sem(ytoplot), 'Color', 'b'});
[~, p] = ttest(ytoplot(:,1), ytoplot(:,2));
lt_plot_text(3, max(ytoplot(:))*1.1, ['SIGlearn, p=' num2str(p)], 'b');

% ============ NOT SIGNIFICANT
ytoplot= Y_lohi(~logical(SigLearnAll),:);
plot([6 7], ytoplot', '-m');
lt_plot([6 7]+0.2, mean(ytoplot), {'Errors', lt_sem(ytoplot), 'Color', 'm'});
[~, p] = ttest(ytoplot(:,1), ytoplot(:,2));
lt_plot_text(6, max(ytoplot(:))*1.1, ['NoSIGlearn, p=' num2str(p)], 'm');


if removeIfNotHigherTargDensity==0
    % ==================== same, but removing syls without diff in density.
    lt_subplot(3,2,6); hold on;
    title('removing syls w/o diff in density');
    xlabel('Low targ dens -- Hi Targ dens');
    ylabel('ff dev, nontarg, in bin');
    
    Y_lohi_CUT = Y_lohi;
    Y_lohi_CUT(indsbad,:) = [];
    
    plot([1 2], Y_lohi_CUT', '-', 'Color', [0.7 0.7 0.7]);
    lt_plot([1.1 2.1], mean(Y_lohi_CUT), {'Errors', lt_sem(Y_lohi_CUT)});
    
    [~, p] = ttest(Y_lohi_CUT(:,1), Y_lohi_CUT(:,2));
    lt_plot_pvalue(p, 'ttest', 1);
    
    lt_plot_zeroline;
end


