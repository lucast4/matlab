function lt_neural_Coher_Summ_Plot2(MotifLists, BregionPair, MOTIFSTATS_Compiled, ...
    All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, All_chanpair, ...
    All_bregionpair, All_bregionpair_alphaorder, tbins, ffbins, ffbinsedges, ...
    normtomeancoh)
%% lt 10/12/18 - break out by motifs in order.

% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% tbins = PARAMS.tbins;
% ffbins = PARAMS.ffbins;
pcols = lt_make_plot_colors(length(ffbinsedges)-1, 1, [1 0 0]);


%%

func = @(X)length(X{2});
motiftot = sum(cellfun(func, MotifLists));
CohMeanAll = nan(length(tbins), length(ffbins), motiftot); % one for each motif
BirdnumAll = nan(motiftot, 1);
countmotif = 1;
for k=1:length(MotifLists)
    
    bthis = MotifLists{k}{1};
    bnumthis = find(strcmp({MOTIFSTATS_Compiled.birds.birdname}, bthis));
    motiflistthis = MotifLists{k}{2};
    
   
    
    figcount=1;
    subplotrows=2;
    subplotcols=length(motiflistthis);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ################### PLOT COHEROGRAMS
    % ================ plot each motif one by one
    for m=1:length(motiflistthis)
        motifthis = motiflistthis{m};
        indthis = All_birdnum==bnumthis & strcmp(All_bregionpair_alphaorder, BregionPair) & ...
            strcmp(All_motifname, motifthis);
        if ~any(indthis)
            countmotif = countmotif+1;
            continue
        end
        
        % ============ extract coherence
        cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indthis));
        
        if normtomeancoh ==1
            cohmat = lt_neural_Coher_QUICK_NormCohGram(cohmat);
        end
        
        % ################## PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(motifthis);
        ylabel(bthis);
        
        cohmean = nanmean(cohmat,3);
        imagesc(tbins, ffbins, cohmean');
        axis tight;
        line([0 0], ylim, 'Color', 'k');
        
        CohMeanAll(:,:,countmotif) = cohmean;
        BirdnumAll(countmotif) = bnumthis;
        countmotif = countmotif+1;
    end
    
    
    % ############# PLOT SMOOTHED POWER
    % ================ plot each motif one by one
    for m=1:length(motiflistthis)
        motifthis = motiflistthis{m};
        indthis = All_birdnum==bnumthis & strcmp(All_bregionpair_alphaorder, BregionPair) & ...
            strcmp(All_motifname, motifthis);
        if ~any(indthis)
            continue
        end
        
        % ============ extract coherence
        cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indthis));
        
        if normtomeancoh ==1
            cohmat = lt_neural_Coher_QUICK_NormCohGram(cohmat);
        end
        
        % ################## PLOT
        % ========= PLOT MEAN TRACE IN MULTIPLE FF BINS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(motifthis);
        ylabel(bthis);
        
        if normtomeancoh==1
            lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', [-0.1 0.1], 0, 0, ...
            ffbinsedges);
        else
            lt_neural_Coher_Plot(cohmat, tbins, ffbins, 2, '-', [0.2 0.8], 0, 0, ...
            ffbinsedges);
        end
        lt_plot_zeroline;
        line([0 0], ylim, 'Color', 'k');
    end
end


%% ======== for each bird, one summary across all syls
% show mean and each trace

lt_figure; hold on;

% ################# USNING MEAN COH
% ===== 1) mean coh
lt_subplot(3,3,1); hold on;
title('mean coh (each motif =1)');
cohmean = nanmean(CohMeanAll,3);
imagesc(tbins, ffbins, cohmean');
% colorbar;
axis tight;

% ==== 2) mean coh, ff bands
lt_subplot(3,3,2); hold on;
lt_neural_Coher_Plot(CohMeanAll, tbins, ffbins, 2, ':', [0.2 0.8], 0, 1, ...
    ffbinsedges);

lt_subplot(3,3,3); hold on;
lt_neural_Coher_Plot(CohMeanAll, tbins, ffbins, 2, '-', [0.2 0.8], 0, 0, ...
    ffbinsedges);


% ################ MODULATION OF COH (SUBTRACT MEAN COH)
CohMeanAll_norm = lt_neural_Coher_QUICK_NormCohGram(CohMeanAll);
% ===== 1) mean coh
lt_subplot(3,3,4); hold on;
title('mean coh (each motif =1)');
ylabel('heat = mean(acrosstime) subtracted');
cohmean = nanmean(CohMeanAll_norm,3);
imagesc(tbins, ffbins, cohmean');
% colorbar;
axis tight;

% ==== 2) mean coh, ff bands
lt_subplot(3,3,5); hold on;
lt_neural_Coher_Plot(CohMeanAll_norm, tbins, ffbins, 2, ':', [-0.2 0.2], 0, 1, ...
    ffbinsedges);
lt_subplot(3,3,6); hold on;
lt_neural_Coher_Plot(CohMeanAll_norm, tbins, ffbins, 2, '-', [-0.1 0.1], 0, 0, ...
    ffbinsedges);










