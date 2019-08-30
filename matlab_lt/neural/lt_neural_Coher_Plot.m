function lt_neural_Coher_Plot(cohmat, tbins, ffbins, plottype, linestyle, clim, ...
    plotsig, plotindivtraces, ffbinsedges, colorbarloc, ffbinsedges_indstoplot, ...
    pcol, plotPval)
%% lt 10/9/18 - various plots for coherence
if ~exist('pcol', 'var')
    pcol = [];
end

if ~exist('plotPval', 'var')
    plotPval=1;
end
% cohmat can be multiple trials (dim 3 is trials).
%%

if ~exist('colorbarloc', 'var')
    colorbarloc = 'EastOutside';
end

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
% 3 - heat map of p values

% plottype 1:
if ~exist('clim', 'var')
    clim = [0.1 0.8]; % limits for coherence
    % elseif imsempty(clim)
    %     clim = [0.1 0.8];
end

% plottype 2:
if ~exist('ffbinsedges', 'var')
    ffbinsedges = [];
end
if isempty(ffbinsedges)
    ffbinsedges = [15 30 80 120]; % edges, to plot timecourse in frequency bands
    ffbinsedges = [20 30 80 120]; % edges, to plot timecourse in frequency bands
    % ffbinsedges = [15 30 80 150]; % edges, to plot timecourse in frequency bands
end

if ~exist('ffbinsedges_indstoplot', 'var')
    ffbinsedges_indstoplot = 1:(length(ffbinsedges)-1); % which one of the intervals to plot?
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
    %     colormap('spring');
    lt_plot_colormap('centered');
    colorbar(colorbarloc);
elseif plottype==2
    
    pcols = lt_make_plot_colors(length(ffbinsedges)-1, 1, [1 0 0]);
    if ~isempty(pcol)
        for j=1:length(pcols)
            pcols{j} = pcol;
        end
    end
    for k=ffbinsedges_indstoplot
        
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
%             [p] = signrank(cohthis');
% try
%             plot(tbins(p<0.05), cohmean(p<0.05)+0.02, '*', 'Color', pcols{k});
%             plot(tbins(p<0.005), cohmean(p<0.005)+0.02, '*', 'Color', pcols{k});
%             plot(tbins(p<0.0005), cohmean(p<0.005)+0.02, '*', 'Color', pcols{k});
% catch err
%     try
%             plot(tbins(p<0.05), cohthis(p<0.05, 1)+0.02, '*', 'Color', pcols{k});
%             plot(tbins(p<0.005), cohthis(p<0.005, 1)+0.02, '*', 'Color', pcols{k});
%             plot(tbins(p<0.0005), cohthis(p<0.005, 1)+0.02, '*', 'Color', pcols{k});
%     catch err
%     end
% end
            try
                plot(tbins(h==1), cohmean(h==1)+0.02, '*', 'Color', pcols{k});
            catch err
                plot(tbins(h==1), cohthis(h==1, 1)+0.02, '*', 'Color', pcols{k});
            end
            if plotPval ==1
                for j=1:length(p)
                    if p(j)<0.1
                        lt_plot_text(tbins(j), cohthis(j, 1), ['p=' num2str(p(j))], 'm', 8);
                    end
                end
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
    if ~isempty(clim)
        ylim(clim);
    end
    line([0 0], ylim, 'Color', 'k');
    
elseif plottype==3
    % FOR each t, ff bin, compute p-value (i.e. deviation from 0, distr
    % over dim 3)
    
    nrows = size(cohmat,1);
    ncols = size(cohmat,2);
    pmat = nan(nrows, ncols);
    for rr=1:nrows
        for cc=1:ncols
            tmp = squeeze(cohmat(rr, cc, :));
            if all(isnan(tmp))
                pmat(rr, cc) = nan;
            else
                pmat(rr, cc) = signrank(squeeze(cohmat(rr, cc, :)));
            end
        end
    end
    
    pmat = log10(pmat);
    
    imagesc(tbins, ffbins, pmat', [-4 0]);
    axis tight;
    line([0 0], ylim, 'Color', 'k');
    lt_plot_colormap('pval');
    %     colorbar('East');
    colorbar(colorbarloc);
    
    % --- mark places with pval <0.05
    pmatsig = pmat<log10(0.05);
    for rr=1:length(tbins)
        for cc=1:length(ffbins)
            if pmatsig(rr, cc)==1
                plot(tbins(rr), ffbins(cc), 'or');
            end
            
        end
    end
    
end
