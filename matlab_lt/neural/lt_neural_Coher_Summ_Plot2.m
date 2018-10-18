function lt_neural_Coher_Summ_Plot2(MotifLists, BregionPair, MOTIFSTATS_Compiled, ...
    All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, All_chanpair, ...
    All_bregionpair, All_bregionpair_alphaorder, tbins, ffbins, ffbinsedges)
%% lt 10/12/18 - break out by motifs in order.

% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);


pcols = lt_make_plot_colors(length(ffbinsedges)-1, 1, [1 0 0]);

%%
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
            continue
        end
        % ============ extract coherence
        cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indthis));
        
        % ################## PLOT
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(motifthis);
        ylabel(bthis);
        
        cohmean = nanmean(cohmat,3);
        imagesc(tbins, ffbins, cohmean');
        axis tight;
        line([0 0], ylim, 'Color', 'k');
        
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
        
        % ################## PLOT
        % ========= PLOT MEAN TRACE IN MULTIPLE FF BINS
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(motifthis);
        ylabel(bthis);
        for kk=1:length(ffbinsedges)-1
            
            indsff = ffbins>ffbinsedges(kk) & ffbins<=ffbinsedges(kk+1);
            ffmean = nanmean(ffbins(indsff));
            % ff bins
            cohthis = squeeze(nanmean(cohmat(:, indsff, :), 2)); % first take mean over the ff bins
            %                cohthis = squeeze(cohmat(:, indsff, :));
            cohmean = nanmean(cohthis,2); % then take mean across trials
            cohsem = lt_sem(cohthis');
            if length(cohsem)==1
                plot(tbins, cohmean, 'Color', pcols{kk});
            else
                shadedErrorBar(tbins, cohmean, cohsem, {'Color', pcols{kk}}, 1);
            end
            lt_plot_text(tbins(end), cohmean(end), [num2str(ffmean)], pcols{kk}, 10);
        end
        axis tight;
        ylim([0.2 0.8]);
        line([0 0], ylim, 'Color', 'k');
    end
end

