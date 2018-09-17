function lt_neural_v2_ANALY_FRsmooth_Plot2(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, corrwindow, dontclosefig, ignoreDiff)
%% LT 9/12/18 - plots all things for each syl (e.g. FR, deviations, etc)
if ~exist('dontclosefig', 'var')
    dontclosefig = 0;
end
    
%%
numbirds = max(OUTDAT.All_birdnum);
numexpts = max(OUTDAT.All_exptnum);
numswitch = max(OUTDAT.All_swnum);
numneur = max(OUTDAT.All_neurnum);

figcount=1;
subplotrows=6;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:numbirds
    for ii=1:numexpts
        for ss =1:numswitch
            
            
            if ~any(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss)
                continue
            end
            
            for nn=1:numneur
                
                bname = SwitchStruct.bird(i).birdname;
                ename = SwitchStruct.bird(i).exptnum(ii).exptname;
                
                % =========== targ syl
                istarg_this = 1;
                issame_this = 1;
                [fignums_alreadyused, hfigs, figcount] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, corrwindow);
                
                % ========== same types
                istarg_this = 0;
                issame_this = 1;
                [fignums_alreadyused, hfigs, figcount] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, corrwindow);
                
                % ========== diff types
                if ignoreDiff==0
                istarg_this = 0;
                issame_this = 0;
                [fignums_alreadyused, hfigs, figcount] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, corrwindow);
                end
            end
            
            pause;
            if dontclosefig==0
                close all;
            end
        end
    end
end


end



function [fignums_alreadyused, hfigs, figcount] = ...
    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, corrwindow)


indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==istarg_this & OUTDAT.All_issame==issame_this);
if isempty(indsthis)
    return
end



for j=indsthis'
    mm = OUTDAT.All_motifnum(j);
    motifthis = ...
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif_regexpr_str{mm};
    
    
    hsplots = [];
    hsplots2 = [];
    % ============ BASELINE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('baseline (r = hi)');
    hsplots = [hsplots hsplot];
    ylabel({[bname], [ename] ['sw' num2str(ss) ', neur' num2str(nn)]});
    indtmp = 1;
    
    frmat = OUTDAT.All_FRsmooth{j, indtmp};
    tbin = OUTDAT.All_FRsmooth_t{j, indtmp};
    ff = OUTDAT.All_FF{j, indtmp};
    thresh = median(ff);
    
    if ~all(isnan(ff))
        % -- LOW FF
        indbyff = ff<thresh;
        pcol = 'b';
        ymean = mean(frmat(:, indbyff), 2);
        ysem = lt_sem(frmat(:, indbyff)');
        shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
%         plot(tbin, OUTDAT.AllBase_FRsmooth_Lo{j}, '-m');
        ylo = ymean;
        
        % -- HI FF
        indbyff = ff>thresh;
        pcol = 'r';
        ymean = mean(frmat(:, indbyff), 2);
        ysem = lt_sem(frmat(:, indbyff)');
        shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
%         plot(tbin, OUTDAT.AllBase_FRsmooth_Hi{j}, '-m');
        yhi = ymean;
    else
        % - then don't separate by FF
        pcol = [0.7 0.7 0.7];
        ymean = mean(frmat, 2);
        ysem = lt_sem(frmat');
        shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
    end
    % --- format
    axis tight;
    YLIM = ylim;
    ylim([0 YLIM(2)]);
    line([0 0], ylim, 'Color', 'k');
    
    % ---- save baseline for later
    ymean_base = mean(frmat, 2);
    ybase_himinuslo = yhi - ylo;
    
    
    %     % =========== BASELINE, HIGH PITCH SUBTRACT LOW PITCH
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title('baseline (hi - lo)');
    %     hsplots2 = [hsplots2 hsplot];
    %
    %     plot(tbin, ybase_himinuslo, '-k', 'LineWidth', 2);
    %     % --- format
    %     axis tight;
    %     lt_plot_zeroline;
    %     lt_plot_zeroline_vert;
    
    % ========== TRAINING
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title({'training', 'base(k), early(m), late(r)'});
    hsplots = [hsplots hsplot];
    indtmp = 2;
    frmat = OUTDAT.All_FRsmooth{j, indtmp};
    tbin = OUTDAT.All_FRsmooth_t{j, indtmp};
    
    ntrial = size(frmat,2);
    midtrial = round(ntrial/2);
    
    % ---- OVERLAY MEAN OF BASELINE
    plot(tbin, ymean_base, '-k');
    
    % ---- EARLY TRAINING
    pcol = 'm';
    ymean = mean(frmat(:, 1:midtrial-1), 2);
    ysem = lt_sem(frmat(:, 1:midtrial-1)');
    shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
    ylearnearly = ymean - ymean_base;
    ylearnearly_sem = ysem;
    
    % ---- LATE TRAINING
    pcol = 'r';
    ymean = mean(frmat(:, midtrial:end), 2);
    ysem = lt_sem(frmat(:, midtrial:end)');
    shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
    % ---- SAVE FOR LATER
    ylearn = ymean - ymean_base;
    ylearn_sem = ysem;
    
    % --- format
    axis tight;
    YLIM = ylim;
    ylim([0 YLIM(2)]);
    line([0 0], ylim, 'Color', 'k');
    linkaxes(hsplots, 'xy');
    lt_plot_zeroline_vert;
    
    
    % ========================== DEV FROM BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots2 = [hsplots2 hsplot];
    title('train - base');
    % === early
    shadedErrorBar(tbin, ylearnearly, ylearnearly_sem, {'Color', 'm'},1);
    
    % =-=== late
    shadedErrorBar(tbin, ylearn, ylearn_sem, {'Color', 'r'}, 1);
    
    % === overlay baseline (hi minus low)
    plot(tbin, ybase_himinuslo, '--k', 'LineWidth', 2);
    
    axis tight;
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([corrwindow(1) corrwindow(1)], ylim, 'Color', 'b');
    line([corrwindow(2) corrwindow(2)], ylim, 'Color', 'b');
    
    % ========================= DERIVED MEASURE, TAKING INTO ACCOUNT DRIFT
    analytype = 'AllDevDiff_NotAbs';
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots2 = [hsplots2 hsplot];
    title(['NormForDrift: ' analytype]);
    
    y = OUTDAT.(analytype){j};
    plot(tbin, y, 'Color', 'r', 'LineWidth', 2);
    
    % --- overlay baseline
    plot(tbin, ybase_himinuslo, '--k', 'LineWidth', 2);
    
    % --
    axis tight;
    lt_plot_zeroline;
    linkaxes(hsplots2, 'xy');
    lt_plot_zeroline_vert;
    % --- put line for corr window
    line([corrwindow(1) corrwindow(1)], ylim, 'Color', 'b');
    line([corrwindow(2) corrwindow(2)], ylim, 'Color', 'b');
    
    % ====================== summarize corrlelation between learning
    % related change and baseline hi minus lo ff
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['correlation']);
    indstocorr = tbin>=corrwindow(1) & tbin<corrwindow(2);
    ylabel('corr');
    
    y1 = corr(ybase_himinuslo(indstocorr), OUTDAT.(analytype){j}(indstocorr));
    y2 = corr(ybase_himinuslo(indstocorr), ylearnearly(indstocorr));
    y3 = corr(ybase_himinuslo(indstocorr), ylearn(indstocorr));
    lt_plot_bar([1 2 3], [y1 y2 y3]);
    xlim([0 4]);
    ylim([-0.8 0.8]);
    lt_plot_zeroline;
    
    xlabthis = {analytype, 'diff(early)', 'diff(late)'};
    set(gca, 'XTickLabel', xlabthis);
    
    
    % ====================== PLOT FF
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    if istarg_this ==1
        pcol = 'k';
    elseif issame_this==1
        pcol = 'b';
    else
        pcol = 'r'; % diff
    end
    title(motifthis, 'Color', pcol);
    
    ff = cell2mat(OUTDAT.All_FF(j, :));
    t = cell2mat(OUTDAT.All_FF_t(j, :));
    indfirsttrain = length(OUTDAT.All_FF{j,1})+1;
    
    plot(t, ff, 'o', 'Color', pcol);
    
    % -- note down z score learning
    pitchz = OUTDAT.AllMinusBase_PitchZ(j, end);
    lt_plot_text(t(end), ff(end), ['z=' num2str(pitchz)], 'r');
    
    % ---- if is targ, determine direction of WN
    if istarg_this==1
        lconting = SwitchStruct.bird(i).exptnum(ii).switchlist(ss).learningContingencies;
        wndat = lconting{find(strcmp(lconting, motifthis))+1};
        lt_plot_text(t(indfirsttrain), min(ff), ['wn:' num2str(wndat)], 'k');
        %                     lt_plot_annotation(1, ['wn:' num2str(wndat)], 'r');
    end
    
    % ---
    axis tight
    tmp = mean(t(indfirsttrain-1:indfirsttrain)); % base-train divider
    line([tmp tmp], ylim);
    
    
    
end
end



