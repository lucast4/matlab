function [OUTSTRUCT] = lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_3(DATBYREND, TrialStruct, ...
    Inds_sylcounter, gettrain, plotraw, edgelist, minrendsinbin, colorscheme, ...
    plotsmooth, logtime, xedges_hand, normtobase, use2bins)

% normtobase =1; % then subtracts baseline (matches) the time bins used by data
% if 2: then norms to catch songs.

% leave empty to not do.
% if no base trials, then wil throw out data entirely.
minbase_rendsperbin = minrendsinbin;

if ~exist('normtobase', 'var')
    normtobase = [];
end

if isempty(normtobase)
    normtobase = 0;
end

if ~exist('use2bins', 'var')
    use2bins = [];
end
if isempty(use2bins)
    use2bins = 0;
end

%% smoothing params


nbin_sm = 20;

%% different ways of entering birsd to get

% birdstoget = 13:17;
% birdstoget = 1;
% birdstoget = 'notSDP'; % e.g. generalization struct
% birdstoget = 'SDP'; % only from SDP experiments


%% params
%     plotraw = 1;
%     edgelist = [2 3 4 5]; % list of edges to use
%     edgelist = [3]; % list of edges to use
%     minrendsinbin = 4;
%     syltype; one of 'same', 'diff', 'targ'
%     gettrain = 2;
%         % 0: only base
%         % 1: only train;
%         % 2: both train and base, side by side
%     getsiglearn = 0;
%     birdstoget = 13:17;
%     birdstoget = 1;
%     colorscheme = 'learnsig';
%     % choices: bird; learnsig

%% lt 7/31/18 -

% divides data into 3 bins. first bin right edge is hand coded. the second
% 2 bins divide the rest of the data up evenly.

%     plotraw = 0; % if 1 then plots for each syl. if 0, then just one raw
%     plot for each value of edgelist\

%     edgelist = [2 3 4 5]; % list of edges to use. each one goes in edges
%     = [0 edgelist(i) median_of_rest last];

%%

if length(edgelist)>1
    % then don't want to plot too many raw plots
    plotraw = 0;
end

% ------- whether to plot only base (0), train (1), or both (1,2)
% ---- if both, then summary will combine both. need to modify this.
if gettrain==0
    GetTrainList = [0];
elseif gettrain==1
    GetTrainList = [1];
elseif gettrain ==2
    GetTrainList = [0 1];
end

%%

if plotsmooth==1
    plotraw=1;
end

%%
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

hsplots_timehist = [];

% --- collect
%%
for edge1 = edgelist
    
    %     Xall = [];
    %     Yall = [];
    %     Yall_sem = [];
    %     Nall = [];
    %     Xedgesall = [];
    
    nrows = length(Inds_sylcounter)*length(GetTrainList);
    
    if isempty(xedges_hand)
        if use2bins==1
            ncols = 3;
        elseif use2bins==0
            ncols = 4; % default bins, see below.
        end
    else
        ncols = length(xedges_hand);
    end
    
    Xall = nan(nrows, ncols);
    Yall = nan(nrows, ncols);
    Yall_sem = nan(nrows, ncols);
    Nall = nan(nrows, ncols);
    Xedgesall = nan(nrows, ncols);
    Tdev_alltrials = cell(nrows, 1);
    FFsm_alltrials = cell(nrows,1);
    Tsm_alltrials = cell(nrows,1);
    
    %     Xall = cell(nrows, 1);
    %     Yall = cell(nrows, 1);
    %     Yall_sem = cell(nrows, 1);
    %     Nall = cell(nrows, 1);
    %     Xedgesall = cell(nrows, 1);
    rowcount = 1;
    if logtime==1
        edge1 = log10(edge1);
    end
    
    for ss = Inds_sylcounter
        
        numGoodCases = 0;
        for gg = GetTrainList
            indsthis = DATBYREND.IsDurTrain==gg & DATBYREND.Sylcounter==ss;
            
            
            if ~any(indsthis)
                rowcount = rowcount+1;
                continue
            end
            
            % ---- want this syl?
            bnum = unique(DATBYREND.Birdnum(indsthis));
            enum = unique(DATBYREND.Exptnum(indsthis));
            snum = unique(DATBYREND.Sylnum(indsthis));
            siglearn = unique(DATBYREND.SigLearn(indsthis));
            
            
            % ============================ PLOT RAW DAT, AND SMOOTHED RUNNING
            tthis = cell2mat(DATBYREND.Time_dev(indsthis));
            tthis = tthis*(24*60);
            ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
            
            
            % ========== NOTE DOWN LEARNING RATE (SHIFT PER REND)
            % ---- [sum of ff deviations]/[sum of time deviations]
            t_tot = sum(tthis);
            ff_tot = sum(ffthis);
            
            
            % ====================== COLLECT BINNED DATA
            if logtime==1
                tthis = log10(tthis);
            end
            
            
            % --------- get 3 bins
            if isempty(xedges_hand)
                if use2bins==1
                    xedgethis = [min(tthis) edge1 max(tthis)];
                else
                    xedgethis = [min(tthis) edge1 prctile(tthis(tthis>edge1), [50]) max(tthis)];
                end
            else
                xedgethis = xedges_hand;
                if logtime==1
                    xedgethis = [min(tthis) log10(xedgethis(2:end))];
                    xedgethis(2) = max(xedgethis(1:2));
                end
            end
            

            % -------------------------- BIN DATA
            tbins = discretize(tthis, xedgethis);
            
            % ------ collect binned balues
            BinID = unique(tbins);
            tbinned = grpstats(tthis, tbins, {'mean'});
            ffbinned = grpstats(ffthis, tbins, {'mean'});
            ffbinned_sem = grpstats(ffthis, tbins, {'sem'});
            nbinned = tabulate(tbins);
            nbinned = nbinned(:,2);
            
            
            % ------------- COLLECT
            tmptmp = tabulate(tbins);
            if min(tmptmp(:,2))< minrendsinbin
                plotflag=1;
            else
                plotflag = 0;
            end
            
            
            %% ========================= NORMALIZE TO BASE?
            if normtobase~=0
                indsbase = [];
                if normtobase==1
                    indsbase = DATBYREND.IsDurTrain==0 & DATBYREND.Sylcounter==ss;
                elseif normtobase ==2
                    indsbase = DATBYREND.IsDurTrain==0 & DATBYREND.Sylcounter==ss & DATBYREND.IsCatch==1;
                end
                if isempty(indsbase) % skip this syl if lacks base data (and desire base)
                    rowcount = rowcount+1;
                    continue
                end
                
                % ---- EXTRACT BASE DATA
                tthisBASE = cell2mat(DATBYREND.Time_dev(indsbase));
                tthisBASE = tthisBASE*(24*60);
                if logtime==1
                    tthisBASE = log10(tthisBASE);
                end
                
                % -------- BIN DATA
                tbinsBASE = discretize(tthisBASE, xedgethis);
                
                ffthisBASE = cell2mat(DATBYREND.FF_dev(indsbase));
                ffbinnedBASE = grpstats(ffthisBASE, tbinsBASE, {'mean'});
                
                if length(unique(tbinsBASE(~isnan(tbinsBASE)))) ~= length(unique(tbins))
                    rowcount = rowcount+1;
                    continue
                end
                
                
                % make sure all bins are matched between dat and base.
                if ~all(unique(tbinsBASE(~isnan(tbinsBASE))) == unique(tbins))
                    rowcount = rowcount+1;
                    continue
                end
                
                % make sure pass minumum samp size
                if min(grpstats(ffthisBASE, tbinsBASE, {'numel'}))<minrendsinbin
                    rowcount = rowcount+1;
                    continue
                end
                
                
                % =================== FINALLY, perform normalization
                ffbinned = ffbinned - ffbinnedBASE;
                
            end
            
            
            %% =================== PLOT
            birdname = TrialStruct.birds(bnum).birdname;
            exptname = TrialStruct.birds(bnum).exptnum(enum).exptname;
            sylthis = TrialStruct.birds(bnum).exptnum(enum).sylnum(snum).syl;
            
            if strcmp(colorscheme, 'bird')
                disp('COLOR BK, bird scheme not yet done ...');
                pcol = 'k';
            elseif strcmp(colorscheme, 'learnsig');
                if siglearn==1
                    pcol = 'b';
                elseif siglearn==0
                    pcol = 'r';
                end
            end
            
            
            
            % ============ PLOT RAW
            if plotraw==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' exptname '-' sylthis]);
                ylabel(['durWN?: ' num2str(gg)]);
                
                lt_plot(tbinned, ffbinned, {'Errors', ffbinned_sem, 'Color', 'k'});
                
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
                
                
                % ---------- put lines for bin edges
                for j=1:length(xedgethis)
                    line([xedgethis(j) xedgethis(j)], ylim, 'Color', 'r');
                end
                
                % --------- put individual datapoints
                plot(tthis, ffthis, '.', 'Color', pcol);
                
                % --- flag showing low N
                if plotflag==1
                    lt_plot_annotation(1, 'low N! - excluding', 'c');
                end
                
                % ----- annotate total elarnig
                lt_plot_annotation(2, ['[' num2str(ff_tot) 'hz]/[' num2str(t_tot) 'sec]'], 'r');
            end
            
            
            % ================ running smoother
            % - first sort
            [~, indsort] = sort(tthis);
            tthis = tthis(indsort);
            ffthis = ffthis(indsort);
            if length(tthis)<1.5*nbin_sm
                rowcount = rowcount+1;
                continue
            end
            
            % --- get runing
            tthis_sm = lt_running_stats(tthis, nbin_sm);
            ffthis_sm = lt_running_stats(ffthis, nbin_sm);
            
            if plotsmooth ==1
                % =================== PLOT
                shadedErrorBar(tthis_sm.Mean, ffthis_sm.Mean, ffthis_sm.SEM, {'Color', pcol}, 1);
                %                 lt_plot(tthis.Median, ffthis.Mean, {'Errors', ffthis.SEM, 'Color', pcol});
                
                lt_plot_zeroline;
                lt_plot_zeroline_vert;
            end
            
            
            
            % ########################################################
            % ######## HISTOGRAMS OF TIME INTERVALS
            if plotraw==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots_timehist = [hsplots_timehist hsplot];
                title([birdname '-' exptname '-' sylthis]);
                ylabel(['durWN?: ' num2str(gg)]);
                xlabel('time intervals');
                
                lt_plot_histogram(tthis);
                for j=1:length(xedgethis)
                    line([xedgethis(j) xedgethis(j)], ylim, 'Color', 'r');
                end
            end
            
            %% ==== collect?
            if (0)
                
                if plotflag==1 % then don't collect
                    continue
                end
                %             tmp = size(Yall,1);
                Xall = [Xall; tbinned'];
                Yall = [Yall; ffbinned'];
                Yall_sem = [Yall_sem; ffbinned_sem'];
                Nall = [Nall; nbinned'];
                Xedgesall = [Xedgesall; xedgethis];
                %             if size(Yall,1) == tmp
                %                 keyboard
                %             end
            else
                
                if plotflag==1 % then don't collect
                    rowcount = rowcount+1;
                    
                    continue
                end
                
                
                %                 Xall{rowcount} = tbinned';
                %                 Yall{rowcount}  = ffbinned';
                %                 Yall_sem{rowcount}  = ffbinned_sem';
                %                 Nall{rowcount}  = nbinned';
                %                 Xedgesall{rowcount}  = xedgethis;
                Xall(rowcount, BinID) = tbinned';
                Yall(rowcount, BinID)  = ffbinned';
                Yall_sem(rowcount, BinID)  = ffbinned_sem';
                Nall(rowcount, BinID)  = nbinned(nbinned>0)';
                Xedgesall(rowcount, :)  = xedgethis;
                Tdev_alltrials{rowcount} = tthis;
                
                Tsm_alltrials{rowcount} = tthis_sm;
                FFsm_alltrials{rowcount} = ffthis_sm;
                
                
                rowcount = rowcount+1;
            end
            numGoodCases = numGoodCases+1;
        end
        
        
        % ###################################### COMPARE BASELINE TO TRAIN
        if numGoodCases==2
            if length(GetTrainList)==2
                
                % ----
                assert(GetTrainList(1) ==0);
                assert(GetTrainList(2) == 1); % for below, assumes is firsrt base then train
                
                % =============== 1) COMPARE FF DEVIATIONS (MEAN) IN DIFFERENT BINS
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' exptname '-' sylthis]);
                ylabel('[mean ffdev]k=base; r=train');
                xlabel('mean timedev');
                
                if (0)
                    xbase = Xall(end-1, :); %
                    xtrain = Xall(end, :);
                    
                    ybase = Yall(end-1, :);
                    ytrain = Yall(end, :);
                    
                    ysembase = Yall_sem(end-1, :);
                    ysemtrain = Yall_sem(end, :);
                else
                    xbase = Xall(rowcount-2, :); %
                    xtrain = Xall(rowcount-1, :);
                    
                    ybase = Yall(rowcount-2, :);
                    ytrain = Yall(rowcount-1, :);
                    
                    ysembase = Yall_sem(rowcount-2, :);
                    ysemtrain = Yall_sem(rowcount-1, :);
                    
                end
                lt_plot(xbase, ybase, {'Errors', ysembase, 'Color', 'k'});
                lt_plot(xtrain, ytrain, {'Errors', ysemtrain, 'Color', 'r'});
                lt_plot_zeroline;
                
                
                % =============== 1) COMPARE TOTAL FFDEV IN EACH BIN
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                %             hsplots = [hsplots hsplot];
                title([birdname '-' exptname '-' sylthis]);
                ylabel('[total ffdev] k=base; r=train');
                xlabel('mean timedev');
                if (0)
                    xbase = Xall(end-1, :); %
                    xtrain = Xall(end, :);
                    
                    ybase = Yall(end-1, :).*Nall(end-1, :);
                    ytrain = Yall(end, :).*Nall(end,:);
                else
                    xbase = Xall(rowcount-2, :); %
                    xtrain = Xall(rowcount-1, :);
                    ybase = Yall(rowcount-2, :).*Nall(rowcount-2, :);
                    ytrain = Yall(rowcount-1, :).*Nall(rowcount-1,:);
                end
                
                lt_plot_bar(xbase, ybase, {'Color', 'k'});
                lt_plot_bar(xtrain, ytrain, {'Color', 'r'});
                
                
                
                % ============== COMPARE MEAN FF DEV OVER ALL TIME DEV
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-' exptname '-' sylthis]);
                ylabel('[mean ffdev, all]');
                xlabel('BASE -- TRAIN');
                
                Ytmp = {};
                % ----- BASE
                indsthis = DATBYREND.IsDurTrain==0 & DATBYREND.Sylcounter==ss;
                ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
                
                Ytmp{1} = ffthis;
                
                % ---- TRAIN
                indsthis = DATBYREND.IsDurTrain==1 & DATBYREND.Sylcounter==ss;
                ffthis = cell2mat(DATBYREND.FF_dev(indsthis));
                
                Ytmp{2} = ffthis;
                
                lt_plot_MultDist(Ytmp, [1 2]);
                xlim([0 3]);
                lt_plot_zeroline;
                
                
                % ============== 2) OVERLAY RUNNING AVERAGES
                
            end
        end
    end
    
    if plotraw ==1
        try
            linkaxes(hsplots, 'xy');
            linkaxes(hsplots_timehist, 'xy');
        catch err
        end
    end
    
    
    %% plot summary for this edge value.
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['summary, each syl, edge1 = ' num2str(edge1) 'min']);
    
    for j=1:size(Xall,1)
        plot(Xall(j, :), Yall(j,:), '-o');
        %         plot(Xall{j}, Yall{j}, '-o');
    end
    
    Ymean = nanmean(Yall,1);
    Ysem = lt_sem(Yall);
    Xmean = nanmean(Xall, 1);
    lt_plot_bar(Xmean, Ymean, {'Errors', Ysem, 'Color', 'k'});
    
    % --- test each bin vs. 0
    for j=1:size(Yall,2)
        if any(~isnan(Yall(:,j)))
            %        [~, p]= ttest(Yall(:,j));
            p= signrank(Yall(:,j));
            lt_plot_text(Xmean(j), Ymean(j), [num2str(p)], 'm');
        end
    end
    
    % ---- test bin2 vs 3
    %     [~, p]= ttest(Yall(:,2), Yall(:,3));
    if size(Yall,2)>3
        p = signrank(Yall(:,2), Yall(:,3));
        lt_plot_text(mean(Xmean(2:3)), max(Yall(:,3)), ['2vs3, p=' num2str(p)], 'r');
    end
end



%% ========== save for output
OUTSTRUCT.Tmean_binned = Xall;
OUTSTRUCT.FFdevMean_binned = Yall;

OUTSTRUCT.FFdevSEM_binned  = Yall_sem;
OUTSTRUCT.N_binned  = Nall;
OUTSTRUCT.Xedgesall  = Xedgesall;
OUTSTRUCT.Tdev_alltrials = Tdev_alltrials;
OUTSTRUCT.Tsm_alltrials = Tsm_alltrials;
OUTSTRUCT.FFsm_alltrials = FFsm_alltrials;











