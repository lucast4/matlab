function lt_neural_Coher_Plot(cohmat, tbins, ffbins, plottype, linestyle, clim)
%% lt 10/9/18 - various plots for coherence
if ~exist('linestyle', 'var')
    linestyle = '-';
end
% plottype:
% 1 - imagesc
% 2 - diff frequency bands (line plot)

% plottype 1:
if ~exist('clim', 'var')
    clim = [0.1 0.8]; % limits for coherence
end

% plottype 2:
ffbinsedges = [15 30 80 150]; % edges, to plot timecourse in frequency bands

%%
if plottype==1
    cohmean = nanmean(cohmat,3);
    if isempty(clim)
    imagesc(tbins, ffbins, cohmean');    
    else
    imagesc(tbins, ffbins, cohmean', clim);
    end
    axis tight;
    line([0 0], ylim, 'Color', 'k');
    
elseif plottype==2
    pcols = lt_make_plot_colors(length(ffbinsedges)-1, 1, [1 0 0]);
    
    for k=1:length(ffbinsedges)-1
        
        indsff = ffbins>ffbinsedges(k) & ffbins<=ffbinsedges(k+1);
        ffmean = mean(ffbins(indsff));
        % ff bins
        cohthis = squeeze(nanmean(cohmat(:, indsff, :), 2)); % first take mean over the ff bins
        %                cohthis = squeeze(cohmat(:, indsff, :));
        cohmean = nanmean(cohthis,2); % then take mean across trials
        cohsem = lt_sem(cohthis');
        if length(cohsem)==1
            plot(tbins, cohmean, 'Color', pcols{k}, 'LineStyle', linestyle);
        else
            shadedErrorBar(tbins, cohmean, cohsem, {'Color', pcols{k}, 'LineStyle', linestyle}, 1);
        end
        lt_plot_text(tbins(end), cohmean(end), [num2str(ffmean)], pcols{k}, 10);
    end
    axis tight;
    ylim(clim);
    line([0 0], ylim, 'Color', 'k');
    
end
