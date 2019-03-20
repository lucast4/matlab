function lt_neural_v2_ANALY_FRsmooth_Plot(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, epochtoplot, dopause)

%%
numbirds = max(OUTDAT.All_birdnum);
numexpts = max(OUTDAT.All_exptnum);
numswitch = max(OUTDAT.All_swnum);
numneur = max(OUTDAT.All_neurnum);




for i=1:numbirds
    for ii=1:numexpts
        for ss =1:numswitch
            
            if ~any(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss)
                continue
            end
            figcount=1;
            subplotrows=8;
            subplotcols=3;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            for nn=1:numneur
                
                if ~any(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                        & OUTDAT.All_neurnum==nn)
                    continue
                end
                
                bname = SwitchStruct.bird(i).birdname;
                ename = SwitchStruct.bird(i).exptnum(ii).exptname;
                hsplots = [];
                
                % =========== targ syl
                istarg_this = 1;
                issame_this = 1;
                [fignums_alreadyused, hfigs, figcount, hsplots] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, epochtoplot, hsplots);
                
                % ========== same types
                istarg_this = 0;
                issame_this = 1;
                [fignums_alreadyused, hfigs, figcount, hsplots] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, epochtoplot, hsplots);
                
                % ========== diff types
                istarg_this = 0;
                issame_this = 0;
                [fignums_alreadyused, hfigs, figcount, hsplots] = ...
                    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
                    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
                    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, epochtoplot, hsplots);
                
                if ~isempty(hsplots)
                    linkaxes(hsplots, 'xy');
                end
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            end
            
            if dopause==1
            pause; close all;
            end
        end
    end
end


end

function [fignums_alreadyused, hfigs, figcount, hsplots] = ...
    fn_plotsmooth(OUTDAT, bname, ename, i, ii, ss, nn, istarg_this, ...
    issame_this, fignums_alreadyused, hfigs, figcount, subplotrows, ...
    subplotcols, MOTIFSTATS_Compiled, SwitchStruct, epochtoplot, hsplots)


indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
    & OUTDAT.All_neurnum==nn & OUTDAT.All_istarg==istarg_this & OUTDAT.All_issame==issame_this);
if isempty(indsthis)
    return
end



for j=indsthis'
    mm = OUTDAT.All_motifnum(j);
    motifthis = ...
        MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.neurons(nn).motif_regexpr_str{mm};
    
    
    %     hsplots = [];
    % ========== TRAINING
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    if istarg_this ==1
        pcol = 'k';
    elseif issame_this==1
        pcol = 'b';
    else
        pcol = 'r'; % diff
    end
    title(motifthis, 'Color', pcol);
    ylabel([bname '-' ename '-sw' num2str(ss) '-n' num2str(nn)]);
    
    % ================= BASELINE
    indtmp = 1;
    
    frmat = OUTDAT.All_FRsmooth{j, indtmp};
    tbin = OUTDAT.All_FRsmooth_t{j, indtmp};
    
    % ---- save baseline for later
    indstmp = OUTDAT.AllBase_indsepoch{j};
    ymean_base = mean(frmat(:, indstmp), 2);
    ysem_base = lt_sem(frmat(:, indstmp)');
    % ---- OVERLAY MEAN OF BASELINE
    shadedErrorBar(tbin, ymean_base, ysem_base, {'Color', 'k'}, 1);
    %     plot(tbin, ymean_base, '-k');
    
    
    
    % ============= WN
    indtmp = 2;
    frmat = OUTDAT.All_FRsmooth{j, indtmp};
    pcol = 'r';
    
    % ---- TRAINING
    indstmp = OUTDAT.AllWN_indsepoch{j}(epochtoplot):OUTDAT.AllWN_indsepoch{j}(epochtoplot+1)-1;
    ymean = mean(frmat(:, indstmp), 2);
    ysem = lt_sem(frmat(:, indstmp)');
    shadedErrorBar(tbin, ymean, ysem, {'Color', pcol}, 1);
    
    
    % --- format
    axis tight;
    YLIM = ylim;
    ylim([0 YLIM(2)]);
    line([0 0], ylim, 'Color', 'k');
    %     linkaxes(hsplots, 'xy');
    
    
    
end
end



