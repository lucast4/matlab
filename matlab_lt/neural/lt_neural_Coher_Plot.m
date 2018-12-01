function lt_neural_Coher_Plot(cohmat, tbins, ffbins, plottype, linestyle, clim, ...
    plotsig, plotindivtraces, ffbinsedges)
%% lt 10/9/18 - various plots for coherence

% cohmat can be multiple trials (dim 3 is trials).

%%

if ~exist('plotindivtraces', 'var')
    plotindivtraces = [];
    % if 1, then plots all lines for plot versino 2. doesn't do anything
    % for v1
end
if ~exist('plotsig', 'var')
    % overlays * for p<0.05 (vs. 0);
    plotsig=0;
end
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
if ~exist('ffbinsedges', 'var')
    ffbinsedges = [];
end
if isempty(ffbinsedges)
ffbinsedges = [15 30 80 120]; % edges, to plot timecourse in frequency bands
% ffbinsedges = [15 30 80 150]; % edges, to plot timecourse in frequency bands
end
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
        ffmin = ffbinsedges(k);
        ffmax = ffbinsedges(k+1);
        
        % ff bins
        
        cohthis = squeeze(nanmean(cohmat(:, indsff, :), 2)); % first take mean over the ff bins
        if plotindivtraces==1
            % - plot each datapoint as line
            plot(tbins, cohthis, '-', 'Color', 0.6*pcols{k});
        else
            %                cohthis = squeeze(cohmat(:, indsff, :));
            cohmean = nanmean(cohthis,2); % then take mean across trials
            cohsem = lt_sem(cohthis');
            if length(cohsem)==1
                plot(tbins, cohmean, 'Color', pcols{k}, 'LineStyle', linestyle);
            else
                shadedErrorBar(tbins, cohmean, cohsem, {'Color', pcols{k}, 'LineStyle', linestyle}, 1);
            end
        end
        % -------------- if is multiple trials, then test deviation from 0
        if plotsig==1
            [h, p] = ttest(cohthis');
            try
            plot(tbins(h==1), cohmean(h==1)+0.02, '*', 'Color', pcols{k});    
            catch err
            plot(tbins(h==1), cohthis(h==1, 1)+0.02, '*', 'Color', pcols{k});
            end
        end
        % --
        try
        lt_plot_text(tbins(end), cohmean(end), [num2str(ffmin) '-' num2str(ffmax)], pcols{k}, 10);
        catch err
            lt_plot_text(tbins(end), cohthis(end, 1), [num2str(ffmin) '-' num2str(ffmax)], pcols{k}, 10);
        end
    end
    axis tight;
    ylim(clim);
    line([0 0], ylim, 'Color', 'k');
    
end
