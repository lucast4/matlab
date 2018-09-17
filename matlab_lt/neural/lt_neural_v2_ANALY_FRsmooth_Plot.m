function lt_neural_v2_ANALY_FRsmooth_Plot(OUTDAT, MOTIFSTATS_Compiled, SwitchStruct)
%%



%%
numbirds = max(OUTDAT.All_birdnum);
numexpts = max(OUTDAT.All_exptnum);
numswitch = max(OUTDAT.All_swnum);
numneur = max(OUTDAT.All_neurnum);

figcount=1;
subplotrows=8;
subplotcols=3;
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
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct);
                
                % ========== same types
                istarg_this = 0;
                issame_this = 1;
                [fignums_alreadyused, hfigs, figcount] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct);
                
                % ========== diff types
                istarg_this = 0;
                issame_this = 0;
                [fignums_alreadyused, hfigs, figcount] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct);
                
            end
            
            pause; close all;
        end
    end
end


end

function [fignums_alreadyused, hfigs, figcount] = ...
    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
    subplotcols, MOTIFSTATS_Compiled, SwitchStruct)


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
    % ============ BASELINE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('baseline');
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
        
        % -- HI FF
        indbyff = ff>thresh;
        pcol = 'r';
        ymean = mean(frmat(:, indbyff), 2);
        ysem = lt_sem(frmat(:, indbyff)');
        shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
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
    
    
    
    % ========== TRAINING
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title({'training', 'base(k), early(m), late(r)'});
    hsplots = [hsplots hsplot];
    indtmp = 2;
    frmat = OUTDAT.All_FRsmooth{j, indtmp};
    tbin = OUTDAT.All_FRsmooth_t{j, indtmp};
    ff = OUTDAT.All_FF{j, indtmp};
    
    ntrial = size(frmat,2);
    midtrial = round(ntrial/2);
    
    % ---- OVERLAY MEAN OF BASELINE
    plot(tbin, ymean_base, '-k');
    
    % ---- EARLY TRAINING
    pcol = 'm';
    ymean = mean(frmat(:, 1:midtrial-1), 2);
    ysem = lt_sem(frmat(:, 1:midtrial-1)');
    shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
    
    % ---- LATE TRAINING
    pcol = 'r';
    ymean = mean(frmat(:, midtrial:end), 2);
    ysem = lt_sem(frmat(:, midtrial:end)');
    shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
    
    
    % --- format
    axis tight;
    YLIM = ylim;
    ylim([0 YLIM(2)]);
    line([0 0], ylim, 'Color', 'k');
    linkaxes(hsplots, 'xy');
    
    
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



