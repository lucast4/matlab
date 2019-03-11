%% lt 11/6/18 - plots ...

%% ------------- spectrograms
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_base, tbins, ffbins, 1, '', clim);
title('chan1, base');
ylabel(ylabthis);
axis tight;
xlim(XLIM);
ylim(YLIM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_wn, tbins, ffbins, 1, '', clim);
title('chan1, WN');
axis tight;
xlim(XLIM);
ylim(YLIM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_base, tbins, ffbins, 1, '', clim);
title('chan2, base');
axis tight;
xlim(XLIM);
ylim(YLIM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_wn, tbins, ffbins, 1, '', clim);
title('chan2, WN');
axis tight;
xlim(XLIM);
ylim(YLIM);


%% ==================== WN MINUS BASE

dat1_wnminusbase = dat1_wn-dat1_base;
dat2_wnminusbase = dat2_wn-dat2_base;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_wnminusbase, tbins, ffbins, 1, '', clim);
title('chan1, WN-base');
ylabel(ylabthis);
axis tight;
xlim(XLIM);
ylim(YLIM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_wnminusbase, tbins, ffbins, 3, '', clim);
title('chan1, WN-base');
ylabel(ylabthis);
axis tight;
colorbar('EastOutside')
xlim(XLIM);
ylim(YLIM);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_wnminusbase, tbins, ffbins, 1, '', clim);
title('chan2, WN-base');
ylabel(ylabthis);
xlim(XLIM);
ylim(YLIM);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_wnminusbase, tbins, ffbins, 3, '', clim);
colorbar('EastOutside')
title('chan2, WN-base');
ylabel(ylabthis);
xlim(XLIM);
ylim(YLIM);


%% ------------ PLOT SPECTRA (PICK A TIMEPOINT FOR TBIN)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('chan1, (k=base; r=WN)');
xlabel('f');
ylabel('mean power frac');
tinds = tbins>=timewindowtoplot(1) & tbins<=timewindowtoplot(2);

specall = squeeze(nanmean(dat1_base(tinds, :, :),1));
spec_mean = nanmean(specall, 2);
spec_sem = lt_sem(specall');
shadedErrorBar(ffbins, spec_mean, spec_sem, {'Color', 'k'}, 1);

specall = squeeze(nanmean(dat1_wn(tinds, :, :),1));
spec_mean = mean(specall, 2);
spec_sem = lt_sem(specall');
shadedErrorBar(ffbins, spec_mean, spec_sem, {'Color', 'r'}, 1);
axis tight;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('chan2, (k=base; r=WN)');
xlabel('f');
ylabel('mean power frac');
tinds = tbins>=timewindowtoplot(1) & tbins<=timewindowtoplot(2);

specall = squeeze(nanmean(dat2_base(tinds, :, :),1));
spec_mean = nanmean(specall, 2);
spec_sem = lt_sem(specall');
shadedErrorBar(ffbins, spec_mean, spec_sem, {'Color', 'k'}, 1);

specall = squeeze(nanmean(dat2_wn(tinds, :, :),1));
spec_mean = mean(specall, 2);
spec_sem = lt_sem(specall');
shadedErrorBar(ffbins, spec_mean, spec_sem, {'Color', 'r'}, 1);
axis tight;


%% ------------ PLOT TIMECOURSE IN A FEW FREQUENCY BINS
if plotEachTrial==1

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_wn, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
title('chan1, WN');
ylabel(ylabthis);
axis tight;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_wn, tbins, ffbins, 2, '-', clim+[-addval addval], 1, 1);
title('chan2, WN');
axis tight;
end

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_base, tbins, ffbins, 2, '-', clim, 1, 0);
title('chan1, base');
axis tight;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_wn, tbins, ffbins, 2, '-', clim, 1, 0);
title('chan1, WN');
axis tight;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_base, tbins, ffbins, 2, '-', clim, 1, 0);
title('chan2, base');
axis tight;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_wn, tbins, ffbins, 2, '-', clim, 1, 0);
title('chan2, WN');
axis tight;


% ------ WN MINUS BASE
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat1_wnminusbase, tbins, ffbins, 2, '-', clim, 1, 0);
title('chan1, WN-base');
axis tight;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(dat2_wnminusbase, tbins, ffbins, 2, '-', clim, 1, 0);
title('chan2, WN-base');
axis tight;