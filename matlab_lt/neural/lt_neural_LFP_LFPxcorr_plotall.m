freqslist = LFPXCORR_freqsall{1};
% disp('note: need to replace freqslist with code that extracts from LFPXCORR_freqsall');
% pause;
ind_freq = find(freqslist==freqtoplot);
assert(length(ind_freq)==1, 'only does one freq at a time..');

%% === PLOT SUMMARY FOR EACH SWITCH, ALL SYLS
[indsgrp, inds_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
enametmp = '';
for j=1:length(inds_unique)
    
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    % === nontarg
    indsthis = indsgrp==inds_unique(j);
    sylnames = OUTSTRUCT.motifname(indsthis);
    istarg = OUTSTRUCT.istarg(indsthis);
    
    bname = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsthis))).birdname;
    ename = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsthis))).exptnum(unique(OUTSTRUCT.enum(indsthis))).exptname;
    swnum = unique(OUTSTRUCT.switch(indsthis));
    %     OUTSTRUCT.chanpair(indsthis,:)
    
    % =========== IF THIS IS A NEW EXPERIMENTS, THEN PAUSE AND LET ME LOOK
    % AT PREVIOUS EXPERIMENT...
    if (0) % pause by switch isntad, so ignroe this.
        if isempty(enametmp)
            enametmp = ename;
        elseif strcmp(enametmp, ename)
            % -- ignroe
        elseif ~strcmp(enametmp, ename)
            enametmp = ename;
            % -- pause and wait
            pause
            close all;
        end
    else
        if isempty(enametmp)
            enametmp = swnum;
        elseif enametmp == swnum
            % -- ignroe
        elseif ~(enametmp == swnum)
            enametmp = swnum;
            % -- pause and wait
            pause
            close all;
        end
    end
    
    
    % --
    ccmax_base = cellfun(@(x)x(:,ind_freq,1), LFPXCORR_Base(indsthis), 'UniformOutput', 0);
    cclag_base = cellfun(@(x)x(:,ind_freq,2), LFPXCORR_Base(indsthis), 'UniformOutput', 0);
    
    ccmax_wn = cellfun(@(x)x(:,ind_freq,1), LFPXCORR_WN(indsthis), 'UniformOutput', 0);
    cclag_wn = cellfun(@(x)x(:,ind_freq,2), LFPXCORR_WN(indsthis), 'UniformOutput', 0);
    
    % - 1 plot for each syl
    nsyls = length(ccmax_base);
    for s=1:nsyls
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        if istarg(s)==1
            title([sylnames{s} '[TARG]'], 'Color', 'r');
        else
            title(sylnames{s});
        end
        if s==1
            ylabel([bname '-' ename '-sw' num2str(swnum)]);
        else
            ylabel('max peak (cc, unbiased)');
        end
        xlabel('lag (sec)');
        
        
        % -- base
        pcol= 'k';
        ccthis = ccmax_base{s};
        lagthis = cclag_base{s};
        plot(lagthis, ccthis, '.', 'Color', pcol);
        xmean = mean(lagthis);
        xsem = lt_sem(lagthis);
        ymean = mean(ccthis);
        ysem = lt_sem(ccthis);
        lt_plot(xmean, ymean, {'Errors', ysem, 'Xerrors', xsem, 'Color', pcol});
        
        
        % -- wn
        pcol= 'r';
        ccthis = ccmax_wn{s};
        lagthis = cclag_wn{s};
        plot(lagthis, ccthis, '.', 'Color', pcol);
        xmean = mean(lagthis);
        xsem = lt_sem(lagthis);
        ymean = mean(ccthis);
        ysem = lt_sem(ccthis);
        lt_plot(xmean, ymean, {'Errors', ysem, 'Xerrors', xsem, 'Color', pcol});
    end
    
    %% == overlay all syls on same plot (plot mean base and dur WN).
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('magenta = targ');
    xlabel('lag');
    ylabel('ccmax');
    for s=1:nsyls
        
        X = [];
        Y = [];
        % -- base
        pcol= 'k';
        ccthis = ccmax_base{s};
        lagthis = cclag_base{s};
        
        xmean = mean(lagthis);
        xsem = lt_sem(lagthis);
        ymean = mean(ccthis);
        ysem = lt_sem(ccthis);
        %         lt_plot(xmean, ymean, {'Errors', ysem, 'Xerrors', xsem, 'Color', pcol});
        lt_plot(xmean, ymean, {'Color', pcol});
        X = [X, xmean];
        Y = [Y, ymean];
        
        % -- wn
        pcol= 'r';
        ccthis = ccmax_wn{s};
        lagthis = cclag_wn{s};
        
        xmean = mean(lagthis);
        xsem = lt_sem(lagthis);
        ymean = mean(ccthis);
        ysem = lt_sem(ccthis);
        %         lt_plot(xmean, ymean, {'Errors', ysem, 'Xerrors', xsem, 'Color', pcol});
        lt_plot(xmean, ymean, {'Color', pcol});
        X = [X, xmean];
        Y = [Y, ymean];
        
        % -- lines connectiong
        if istarg(s)==1
            line(X, Y, 'Color', 'm', 'LineWidth', 2);
        else
            line(X, Y, 'Color', 'k');
        end
    end
    YLIM = ylim;
    ylim([0 YLIM(2)]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    
    %% == overlay all syls on same plot (plot mean base and dur WN).
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('using means');
    xlabel('targ = red');
    ylabel('frac change in xcorr');
    
    ccmax_all = [cellfun(@mean, ccmax_base) cellfun(@mean, ccmax_wn)];
    
    ccmax_frac = (ccmax_all(:,2) - ccmax_all(:,1))./ccmax_all(:,1);
    x = 1:length(ccmax_frac);
    lt_plot_bar(x, ccmax_frac);
    lt_plot_bar(x(istarg==1), ccmax_frac(istarg==1), {'Color', 'r'});
    lt_plot_zeroline;
    
    set(gca, 'XTick', x, 'XTickLabel', sylnames);
    rotateXLabels(gca, 90);

    %% == overlay all syls on same plot (plot mean base and dur WN).
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('using means');
    xlabel('targ = red');
    ylabel('change in lag (sec)');
    
%     ccmax_all = [cellfun(@mean, ccmax_base) cellfun(@mean, ccmax_wn)];
    cclag_all = [cellfun(@mean, cclag_base) cellfun(@mean, cclag_wn)];
    
    ccmax_frac = (cclag_all(:,2) - cclag_all(:,1));
    x = 1:length(ccmax_frac);
    lt_plot_bar(x, ccmax_frac);
    lt_plot_bar(x(istarg==1), ccmax_frac(istarg==1), {'Color', 'r'});
    lt_plot_zeroline;
    
    set(gca, 'XTick', x, 'XTickLabel', sylnames);
    rotateXLabels(gca, 90);

    %% == overlay all syls on same plot (plot mean base and dur WN).
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('using medians');
    xlabel('targ = red');
    ylabel('frac change in xcorr');
    
    ccmax_all = [cellfun(@median, ccmax_base) cellfun(@median, ccmax_wn)];
%     cclag_all = [cellfun(@median, cclag_base) cellfun(@median, cclag_wn)];
    
    ccmax_frac = (ccmax_all(:,2) - ccmax_all(:,1))./ccmax_all(:,1);
    x = 1:length(ccmax_frac);
    lt_plot_bar(x, ccmax_frac);
    lt_plot_bar(x(istarg==1), ccmax_frac(istarg==1), {'Color', 'r'});
    lt_plot_zeroline;
    
    set(gca, 'XTick', x, 'XTickLabel', sylnames);
    rotateXLabels(gca, 90);

        %% == overlay all syls on same plot (plot mean base and dur WN).
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('using medians');
    xlabel('targ = red');
    ylabel('change in lag(sec)');
    
%     ccmax_all = [cellfun(@median, ccmax_base) cellfun(@median, ccmax_wn)];
    cclag_all = [cellfun(@median, cclag_base) cellfun(@median, cclag_wn)];
    
    ccmax_frac = (cclag_all(:,2) - cclag_all(:,1));
    x = 1:length(ccmax_frac);
    lt_plot_bar(x, ccmax_frac);
    lt_plot_bar(x(istarg==1), ccmax_frac(istarg==1), {'Color', 'r'});
    lt_plot_zeroline;
    
    set(gca, 'XTick', x, 'XTickLabel', sylnames);
    rotateXLabels(gca, 90);

end
