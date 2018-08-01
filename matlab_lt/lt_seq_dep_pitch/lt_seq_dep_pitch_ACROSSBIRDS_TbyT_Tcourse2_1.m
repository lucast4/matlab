function lt_seq_dep_pitch_ACROSSBIRDS_TbyT_Tcourse2_1(DATBYREND, TrialStruct, getsame, gettrain, ...
    getsiglearn, birdstoget, plotraw, edgelist, minrendsinbin, colorscheme)

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

%%
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

sylmax = max(DATBYREND.Sylcounter);

% --- collect
%%
for edge1 = edgelist
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
        xedgethis = [0 3 prctile(tthis(tthis>edge1), [50]) max(tthis)];
        tbins = discretize(tthis, xedgethis);
        
        % ------ collect binned balues
        tbinned = grpstats(tthis, tbins, {'mean'});
        ffbinned = grpstats(ffthis, tbins, {'mean'});
        ffbinned_sem = grpstats(ffthis, tbins, {'sem'});
        
        
        % ------------- COLLECT
        tmptmp = tabulate(tbins);
        if min(tmptmp(:,2))< minrendsinbin
            plotflag=1;
        else
            plotflag = 0;
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
        
        if plotraw==1
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
            
            % --- flag showing low N
            if plotflag==1
                lt_plot_annotation(1, 'low N! - excluding', 'c');
            end
            
        end
        
        %% ==== collect?
        if plotflag==1 % then don't collect
            continue
        end
        
        Xall = [Xall; tbinned'];
        Yall = [Yall; ffbinned'];
    end
    
    if plotraw ==1
        linkaxes(hsplots, 'x');
    end
    
    %% plot summary for this edge value.
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['summary, each syl, edge1 = ' num2str(edge1) 'min']);
    
    for j=1:size(Xall,1)
        plot(Xall(j, :), Yall(j,:), '-o');
    end
    Ymean = mean(Yall,1);
    Ysem = lt_sem(Yall);
    Xmean = mean(Xall, 1);
    lt_plot_bar(Xmean, Ymean, {'Errors', Ysem, 'Color', 'k'});
    
    % --- test each bin vs. 0
    for j=1:size(Yall,2)
        %        [~, p]= ttest(Yall(:,j));
        p= signrank(Yall(:,j));
        lt_plot_text(Xmean(j), Ymean(j), [num2str(p)], 'm');
    end
    
    % ---- test bin2 vs 3
    %     [~, p]= ttest(Yall(:,2), Yall(:,3));
    p = signrank(Yall(:,2), Yall(:,3));
    lt_plot_text(mean(Xmean(2:3)), max(Yall(:,3)), ['2vs3, p=' num2str(p)], 'r');
end
