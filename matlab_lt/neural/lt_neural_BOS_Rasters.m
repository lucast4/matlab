
% =============
nunits = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset,2);
nsongs = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset, 1);

figcount=1;
subplotrows=4;
subplotcols = 1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for nn=1:nunits
    
    
    spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset(:,nn);
    onsets = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets;
    offsets = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets;
    segextract = SummaryBOS.expt(i).DAT_bysongrend.SegextractFormat.unitnum(nn).segextract;
    
    % ====== 2) rasters
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('rasters [syls shaded]');
    ylabel(['unit' num2str(nn)]);
    for ss=1:nsongs
        
        % -- overlay onsets and offsets
        ons = onsets{ss};
        offs = offsets{ss};
        lt_neural_QUICK_PlotSylPatches(ons, offs, ss);
        lt_neural_PLOT_rasterline(spks{ss}, ss, 'r', 0);
        %
        %                        X=[segextract(j).WNonset_sec  segextract(j).WNoffset_sec ...
        %                     segextract(j).WNoffset_sec  segextract(j).WNonset_sec];
        %                 Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];
        
    end
    
    % ====== 3) mean fr
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('rasters [syls shaded]');
    %     segextract = segextract';
    segextract = lt_neural_SmoothFR(segextract, [], [], [], [], [], PARAMS.flanktotake(1));
    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    plot(x, frmat, 'Color', [0.7 0.7 0.7]);
    frmean = mean(frmat,2);
    frsem = lt_sem(frmat');
    shadedErrorBar(x, frmean, frsem, {'Color', 'k'},1);
    
    % sanity check
    if (0)
        figure; hold on;
        % ----
        fr = segextract(1).FRsmooth_rate_CommonTrialDur;
        t = segextract(1).FRsmooth_xbin_CommonTrialDur;
        %        t = t+SummaryBOS.expt(1).DAT_bysongrend.TEdges(1);
        plot(t, fr, '-k');
        spks = SummaryBOS.expt(1).DAT_bysongrend.SpkTime_RelSongOnset{1,nn};
        lt_neural_PLOT_rasterline(spks, 100, 'r', 0);
        
        % ----
        fr = segextract(1).FRsmooth_rate_CommonTrialDur;
        lt_neural_PLOT_rasterline(spks, 100, 'r', 0);
    end
    
end