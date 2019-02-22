function [lfpcollect, lfpx]= lt_neural_POPLEARN_SylLocked_Over_sub1(chanlist_lfp, ...
    bregion_lfp, DATLFP, indstoplot, xtoplot, fpass, plotraw)

%% 
lfpcollect = [];
for cc=1:length(chanlist_lfp)
    chanthis = chanlist_lfp(cc);
    bregionthis = bregion_lfp{cc};
    
    if strcmp(bregionthis, 'LMAN')
        pcol = [0.3 0.7 0.3];
    elseif strcmp(bregionthis, 'RA')
        pcol = [0.8 0.2 0.2];
    end
    
    % ---- extract lfp data
    lfpall = DATLFP.lfpall(indstoplot, DATLFP.lfpall_chans==chanthis);
    lfpx = DATLFP.t_lfp;
    indx = lfpx>=xtoplot(1) & lfpx<=xtoplot(2);
    
    % -- keep the correct times.
    lfpall = cellfun(@(x)x(indx), lfpall, 'UniformOutput', 0);
    lfpx = lfpx(indx);
    
    % =========== get average over trials of interest.
    lfpall = squeeze(lt_neural_Coher_Cell2Mat(lfpall));
    
    
    if (0)
        %% ==== testing filtering of lfp
        lfpall2 = lt_neural_filter(lfpall, 1500, 0, fpass(1), fpass(2));
        
        figcount=1;
        subplotrows=4;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots = [];
        
        for xx=1:8
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            indrand = randi(size(lfpall,2));
            
            plot(lfpx, lfpall(:,indrand), 'k');
            plot(lfpx, lfpall2(:,indrand), 'r');
            
        end
    end
    % ========== FILTER TO REMOVE LOW F FLULCTURATION
    lfpall = lt_neural_filter(lfpall, 1500, 0, fpass(1), fpass(2));
    
    ymean = mean(lfpall, 2);
    ysem = lt_sem(lfpall');
    
    % -- zscore
    ysem = ysem./std(ymean);
    ymean = (ymean-mean(ymean))./std(ymean);
    
    % ===== flip upside down
    ymeanplot = -ymean;
    
    % --------- SHIFT UP TO NOT OVERLAP NEURAL DATA
    ymeanplot = single(ymeanplot)+4;
    
    
    % ============= PLOT
    if plotraw==1
        shadedErrorBar(lfpx, ymeanplot, ysem, {'Color', pcol},1);
    end
    
    % ======== ANNOTATE
    % --- neuron ID and chan
    lt_plot_text(lfpx(end), ymeanplot(end), ['-ch' num2str(chanthis)], pcol, 8);
    
    % ============= COLLECT OVER CHANNELS
    lfpcollect = [lfpcollect; ymean'];
end