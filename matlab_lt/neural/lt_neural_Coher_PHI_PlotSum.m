function lt_neural_Coher_PHI_PlotSum(OUTSTRUCT, PARAMS, PhiMeanAll, birdtoplot)
%% lt 1/8/19 - plots summary phi-gram and slices across all cases

% PhiMeanAll = OUTSTRUCT.PhiMeanBaseAll;

if ~isempty(birdtoplot)
    indsgood = OUTSTRUCT.bnum==birdtoplot;
    
    PhiMeanAll = PhiMeanAll(:,:, indsgood);
    
end
clear OUTSTRUCT % If want to use, then must subsample first...


%% ====== collapse data?

% if strcmp(datlevel, 

%%

t=PARAMS.tbins;
f = PARAMS.ffbins;
clim = [-pi pi];

%%
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% =========== GLOBAL MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('global mean');

phimean = nan(size(PhiMeanAll,1), size(PhiMeanAll,2));
for j=1:size(PhiMeanAll,1)
    for jj=1:size(PhiMeanAll,2)
        phimean(j, jj) = circ_mean(squeeze(PhiMeanAll(j, jj, :)));
    end
end

imagesc(t, f, phimean', clim);
colorbar;
lt_plot_colormap('centered');
axis tight;


% =========== GLOBAL MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('global mean [smaller clim range]');

phimean = nan(size(PhiMeanAll,1), size(PhiMeanAll,2));
for j=1:size(PhiMeanAll,1)
    for jj=1:size(PhiMeanAll,2)
        phimean(j, jj) = circ_mean(squeeze(PhiMeanAll(j, jj, :)));
    end
end

imagesc(t, f, phimean', [-1.5 1.5]);
lt_plot_colormap('centered');
colorbar;
axis tight;


% =========== GLOBAL MEAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('global mean [smaller clim range]');

phimean = nan(size(PhiMeanAll,1), size(PhiMeanAll,2));
for j=1:size(PhiMeanAll,1)
    for jj=1:size(PhiMeanAll,2)
        phimean(j, jj) = circ_mean(squeeze(PhiMeanAll(j, jj, :)));
    end
end

imagesc(t, f, phimean', [-0.5 0.5]);
lt_plot_colormap('centered');
colorbar;
axis tight;


% ============ PLOT ANGLE SPECTRUM AT DIFFERENT TIMEPOINTS
timebins = PARAMS.tbins(1):0.02:PARAMS.tbins(end);
tbininds = discretize(t, timebins);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('angle spectra at diff time bins');
plotcols = lt_make_plot_colors(length(timebins)-1, 1, [1 0 0]);

for i=1:length(timebins)-1
    indsthis = tbininds==i;
    
    phivec = circ_mean(phimean(indsthis, :));
    
    plot(f, phivec, '-', 'Color', plotcols{i});
    
    lt_plot_text(f(end), phivec(end), ['t=' num2str(timebins(i)) ' to ' num2str(timebins(i+1))], plotcols{i});
    
end
lt_plot_zeroline;
ylim([-pi pi]);

if (0) % OBSOLETE - above plots all on one plot. here plots them on separate plots.
    % ============ PLOT ANGLE SPECTRUM AT DIFFERENT TIMEPOINTS
    timebins = PARAMS.tbins(1):0.02:PARAMS.tbins(end);
    tbininds = discretize(t, timebins);
    for i=1:length(timebins)-1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        indsthis = tbininds==i;
        title(['t=' num2str(timebins(i)) ' to ' num2str(timebins(i+1))]);
        
        phivec = circ_mean(phimean(indsthis, :));
        
        plot(f, phivec, '-k');
        lt_plot_zeroline;
        ylim([-pi pi]);
        
    end
end
