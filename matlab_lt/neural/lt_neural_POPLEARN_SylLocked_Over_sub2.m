function [frsmall, x] = lt_neural_POPLEARN_SylLocked_Over_sub2(neurlist, DAT, bregionlist, ...
    segglobal, xtoplot, plotraw, PARAMS, indstoplot, chanlist)


%%
frsmall = [];
    for nn=1:length(neurlist)
        assert(neurlist(nn)==DAT(nn).neurID_orig);
        seg = DAT(nn).SegmentsExtract;
        
        if strcmp(bregionlist{nn}, 'LMAN')
            pcol = [0.3 0.7 0.3];
        elseif strcmp(bregionlist{nn}, 'RA')
            pcol = [0.8 0.2 0.2];
        end
        
        % =========== extract smoothed fr
        seg = lt_neural_SmoothFR(seg, [], [], [], 0, segglobal);
        
        % =========== PLOT (SPIKES)
        %         lt_neural_PLOT_rasterline(seg.spk_Times, nn, pcol);
        
        % =========== PLOT (smoothed FR)
        x = seg(1).FRsmooth_xbin_CommonTrialDur - PARAMS.motif_predur;
        ymat = [seg(indstoplot).FRsmooth_rate_CommonTrialDur];
        indx = x>=xtoplot(1) & x<=xtoplot(2);
        x = x(indx);
        ymat = ymat(indx,:);
        ymean = mean(ymat,2);
        ysem = lt_sem(ymat');
        
        % -- zscore
        ysem = ysem./std(ymean);
        ymean = (ymean-mean(ymean))./std(ymean);
        if plotraw==1
            shadedErrorBar(x, ymean, ysem, {'Color', pcol},1);
        end
        
        % ======== ANNOTATE
        % --- neuron ID and chan
        chanthis = chanlist(nn);
        if plotraw==1
            lt_plot_text(x(end), ymean(end), ['n' num2str(neurlist(nn)) '-ch' num2str(chanlist(nn))], pcol, 8);
        end
        
        
        % ================ COLLECT ALL SMOOTHED FR
        frsmall = [frsmall; ymean'];
    end