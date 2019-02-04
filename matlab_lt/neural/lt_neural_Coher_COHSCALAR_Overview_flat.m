function lt_neural_Coher_COHSCALAR_Overview_flat(birdthis, OUTSTRUCT, FFcorrCoh, ...
    FFcorrCoh_pval, FFcorrCoh_pctileVsShuff, FFcorrCoh_zscoreVsShuff, FFcorrCoh_shuffCI, ...
    PARAMS, tlims, flims, combineChans)

%% lt 1/8/19 - flatten, for each time,ff bin get stats

YLIM = [18 120];

%% get just this bird

if ~isempty(birdthis)
    indsgood = OUTSTRUCT.bnum==birdthis;
    OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indsgood, 1);
    % FFcorrCoh = FFcorrCoh(:,:, indsgood);
    FFcorrCoh = FFcorrCoh(:,:, indsgood);
    FFcorrCoh_pval = FFcorrCoh_pval(:,:, indsgood);
    % FFcorrCoh_pctileVsShuff = FFcorrCoh_pctileVsShuff(:,:, indsgood);
    FFcorrCoh_zscoreVsShuff = FFcorrCoh_zscoreVsShuff(:,:, indsgood);
    % FFcorrCoh_shuffCI = FFcorrCoh_shuffCI(:,:,:, indsgood);
end

%% =========== PLOTS [NOT COLLAPSING ACROSS ANYTHING (E.G. CHANNELS/EXPERIEMNTS)
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ======= NUMBER SIGNIFCAIT POISTIVE AND NEGATIVE CORRS
nrows = size(FFcorrCoh,1);
ncols = size(FFcorrCoh,2);
nPosCorr = nan(nrows, ncols);
nNegCorr = nan(nrows, ncols);
for rr=1:nrows
    for cc=1:ncols
            
        tmp = FFcorrCoh(rr, cc, squeeze(FFcorrCoh_pval(rr, cc, :)<0.05));
        
        nPosCorr(rr, cc) = sum(tmp>0);
        nNegCorr(rr, cc) = sum(tmp<0);
        
    end
end



% ========= 1) NUMBER SIGNIFICANT POSITIVE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('n significant positive (all cases)');
Y = nPosCorr;

clim = [min(Y(:)) max(Y(:))];
lt_neural_Coher_Plot(Y, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_annotation(1, ['n=' num2str(size(FFcorrCoh,3))], 'm');
lt_plot_colormap('pval');
colorbar('EastOutside');
ylim(YLIM)

% ========= 1) NUMBER SIGNIFICANT POSITIVE [FRACTION]
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('n significant positive (all cases) [frac]');
Y = nPosCorr./size(FFcorrCoh,3);

clim = [min(Y(:)) max(Y(:))];
lt_neural_Coher_Plot(Y, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_annotation(1, ['n=' num2str(size(FFcorrCoh,3))], 'm');
lt_plot_colormap('pval');
colorbar('EastOutside');
ylim(YLIM)


% ========= 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('n significant negative (all cases)');
Y = nNegCorr;

clim = [min(Y(:)) max(Y(:))];
lt_neural_Coher_Plot(Y, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_annotation(1, ['n=' num2str(size(FFcorrCoh,3))], 'm');
lt_plot_colormap('pval');
colorbar('EastOutside');
ylim(YLIM)

% ========= 1) 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('n significant positive (all cases) [frac]');
Y = nNegCorr./size(FFcorrCoh,3);

clim = [min(Y(:)) max(Y(:))];
lt_neural_Coher_Plot(Y, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_annotation(1, ['n=' num2str(size(FFcorrCoh,3))], 'm');
lt_plot_colormap('pval');
colorbar('EastOutside');
ylim(YLIM)


%% ======= [MEAN ACROSS CHANNEL PAIRS FIRST?]

if combineChans>0
    if combineChans==1 % then one val for each unique motif + switch.
        [indsgrp, indsgrp_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum OUTSTRUCT.enum OUTSTRUCT.switch OUTSTRUCT.motifID_unique});
    elseif combineChans==2 % then one val for each unique motif
        [indsgrp, indsgrp_unique] = lt_tools_grp2idx({OUTSTRUCT.bnum OUTSTRUCT.motifID_unique});
    end
    
    nrows = size(FFcorrCoh,1);
    ncols = size(FFcorrCoh,2);
    outmat = nan(nrows, ncols, length(indsgrp_unique));
    outmat_z = nan(nrows, ncols, length(indsgrp_unique));
    for rr=1:nrows
        for cc=1:ncols
            
            tmp = grpstats(squeeze(FFcorrCoh(rr, cc,:)), indsgrp);
            outmat(rr, cc, :) = tmp;
            
            tmp = grpstats(squeeze(FFcorrCoh_zscoreVsShuff(rr, cc,:)), indsgrp);
            outmat_z(rr, cc, :) = tmp;
        end
    end
    
    FFcorrCoh = outmat;
    FFcorrCoh_zscoreVsShuff = outmat_z;
end


%% ====== [PLOT A FEW THINGS] [USING JUST CORRELATIONS VALUES]

figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ========= 1) MEAN CORRELATION
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, (not abs)');
cormean = nanmean(FFcorrCoh, 3);
ymax = max(abs(cormean(:)));
clim = [-ymax ymax];
lt_neural_Coher_Plot(FFcorrCoh, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_annotation(1, ['n=' num2str(size(FFcorrCoh,3))], 'm');
colorbar('EastOutside');
ylim(YLIM)

% ===== 2) P-VALUE OF CORRELATIONS (i.e. distribution deviation from 0?)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('[pval] mean cor coeff, (not abs)');
% cormean = nanmean(FFcorrCoh, 3);
% ymax = max(abs(cormean(:)));
% clim = [-ymax ymax];
lt_neural_Coher_Plot(FFcorrCoh, PARAMS.tbins, PARAMS.ffbins, 3, '', clim);
lt_plot_annotation(1, ['n=' num2str(size(FFcorrCoh,3))], 'm');
colorbar('EastOutside');
ylim(YLIM)


% ========= 1) MEAN ABS CORRELATION
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, abs');
cormean = nanmean(abs(FFcorrCoh), 3);
ymax = max(abs(cormean(:)));
ymin = min(abs(cormean(:)));
clim = [ymin ymax];
lt_neural_Coher_Plot(abs(FFcorrCoh), PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_colormap('pval')
colorbar('EastOutside');
ylim(YLIM)

% ======== 2) MEAN CORRELATION (ZSCORES)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, (zscore(vs shuff), not abs)');
cormean = nanmean(FFcorrCoh_zscoreVsShuff, 3);
ymax = max(abs(cormean(:)));
clim = [-ymax ymax];
lt_neural_Coher_Plot(FFcorrCoh_zscoreVsShuff, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
colorbar('EastOutside');
ylim(YLIM)


% ======== 2) MEAN CORRELATION, ABS(ZSCORES)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('mean cor coeff, abs [zscore]');
cormean = nanmean(abs(FFcorrCoh_zscoreVsShuff), 3);
ymax = max(abs(cormean(:)));
ymin = min(abs(cormean(:)));
clim = [ymin ymax];
lt_neural_Coher_Plot(abs(FFcorrCoh_zscoreVsShuff), PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
lt_plot_colormap('pval')
colorbar('EastOutside');
ylim(YLIM)



% % ===== 2) P-VALUE OF CORRELATIONS (i.e. distribution deviation from 0?)
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% title('[pval] mean cor coeff, (zscore, not abs)');
% lt_neural_Coher_Plot(FFcorrCoh_zscoreVsShuff, PARAMS.tbins, PARAMS.ffbins, 3, '', clim);
% colorbar('EastOutside');
% ylim(YLIM)




