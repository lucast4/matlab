function lt_neural_Coher_Summ_Plot1(pairthis, All_CohgramMean, All_birdnum, All_enum, All_setnum, All_motifname, ...
    All_chanpair, All_bregionpair, All_bregionpair_alphaorder,ffbinsedges)
% pairthis = 'LMAN-RA';
% ffbinsedges = [15 30 80 150]; % edges, to plot timecourse in frequency bands

%% lt 10/12/18 - plots summary over all data
% ============ 1) Lock to syllable onset, one plot for each brain region
% pair
lt_figure; hold on;
indsthis = strcmp(All_bregionpair_alphaorder, pairthis);
cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indsthis));

cohmean = nanmean(cohmat,3);
imagesc(tbins, ffbins, cohmean');
figure; plot(cohmean(:,3), 'ok-')



% ============ 2) ONE PLOT FOR EACH MOTIF AND PAIR TYPE
% --------------- INPUTS
pcols = lt_make_plot_colors(length(ffbinsedges)-1, 1, [1 0 0]);
figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ------------- figure out ff bins

maxbird = max(All_birdnum);
bregionlist = unique(All_bregionpair_alphaorder);
for i=1:maxbird
    
    motiflist = unique(All_motifname(All_birdnum==i));
    birdname = SummaryStruct.birds(i).birdname;
    
    for b=1:length(bregionlist)
        bregionthis = bregionlist{b};
        
        % ==== go thru all motifs.
        for j=1:length(motiflist)
            motifthis = motiflist{j};
            
            indsthis = All_birdnum==i & strcmp(All_motifname, motifthis) ...
                & strcmp(All_bregionpair_alphaorder, bregionthis);
            
            if ~any(indsthis)
                continue
            end
            
            % --------- get mean cohgram
            cohmat = lt_neural_Coher_Cell2Mat(All_CohgramMean(indsthis));
            
            % ========== PLOT MEAN COHGRAM
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([motifthis ',' bregionthis]);
            ylabel(birdname);
            
            cohmean = nanmean(cohmat,3);
            imagesc(tbins, ffbins, cohmean');
            axis tight;
            line([0 0], ylim, 'Color', 'k');
            
            % ========= PLOT MEAN TRACE IN MULTIPLE FF BINS
            %             for k=1:size(cohmat,2)
            %                 % ff bins
            %                cohthis = squeeze(cohmat(:, k, :));
            %                cohmean = nanmean(cohthis,2);
            %                cohsem = lt_sem(cohthis');
            %
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([motifthis ',' bregionthis]);
            
            for k=1:length(ffbinsedges)-1
                
                indsff = ffbins>ffbinsedges(k) & ffbins<=ffbinsedges(k+1);
                ffmean = nanmean(ffbins(indsff));
                % ff bins
                cohthis = squeeze(nanmean(cohmat(:, indsff, :), 2)); % first take mean over the ff bins
                %                cohthis = squeeze(cohmat(:, indsff, :));
                cohmean = nanmean(cohthis,2); % then take mean across trials
                cohsem = lt_sem(cohthis');
                if length(cohsem)==1
                    plot(tbins, cohmean, 'Color', pcols{k});
                else
                    shadedErrorBar(tbins, cohmean, cohsem, {'Color', pcols{k}}, 1);
                end
                lt_plot_text(tbins(end), cohmean(end), [num2str(ffmean)], pcols{k}, 10);
            end
            axis tight;
            ylim([0.2 0.8]);
            line([0 0], ylim, 'Color', 'k');
            
        end
    end
end