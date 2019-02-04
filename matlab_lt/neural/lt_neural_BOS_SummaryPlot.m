
% =============
nunits = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset,2);
% nsongs = size(SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset, 1);

figcount=1;
subplotrows=4;
subplotcols = 1;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for nn=1:nunits
    
    chan = SummaryBOS.expt(i).channels(nn);
    bregion = SummaryBOS.expt(i).bregions{nn};
    spks = SummaryBOS.expt(i).DAT_bysongrend.SpkTime_RelSongOnset(:,nn);
    onsets = SummaryBOS.expt(i).DAT_bysongrend.SylOnsets;
    offsets = SummaryBOS.expt(i).DAT_bysongrend.SylOffsets;
    segextract = SummaryBOS.expt(i).DAT_bysongrend.SegextractFormat.unitnum(nn).segextract;
    bostype = SummaryBOS.expt(i).DAT_bysongrend.BOStype;
    BOSnames = PARAMS.BOSnames{strcmp(PARAMS.BOSbirdname, birdthis)};
    BOSlabels = PARAMS.BOSlabels{strcmp(PARAMS.BOSbirdname, birdthis)};
    
    % ======= collect smoothed fr for each trial
    segextract = lt_neural_SmoothFR(segextract, [], [], [], [], [], PARAMS.flanktotake(1));
    frmat = [segextract.FRsmooth_rate_CommonTrialDur];
    x = segextract(1).FRsmooth_xbin_CommonTrialDur;
    
    %     % ====== 2) rasters
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title('rasters [syls shaded]');
    %     ylabel(['unit' num2str(nn)]);
    %     for ss=1:nsongs
    %
    %        % -- overlay onsets and offsets
    %         ons = onsets{ss};
    %         offs = offsets{ss};
    %         lt_neural_QUICK_PlotSylPatches(ons, offs, ss);
    %        lt_neural_PLOT_rasterline(spks{ss}, ss, 'r', 0);
    % %
    % %                        X=[segextract(j).WNonset_sec  segextract(j).WNoffset_sec ...
    % %                     segextract(j).WNoffset_sec  segextract(j).WNonset_sec];
    % %                 Y=[-j-0.4 -j-0.4 -j+0.4 -j+0.4];
    %
    %     end
    
    % ================= ONE PLOT FOR EACH BOS TYPE
    nbostypes = max(bostype);
    pcols_bos = lt_make_plot_colors(nbostypes, 0,0);
    
    for bb=1:nbostypes
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['unit' num2str(nn) ' [ch' num2str(chan) '-' bregion '] -' BOSnames{bb}]);
        hsplots = [hsplots hsplot];
        
        indsthis = bostype==bb;
        
        frthis = frmat(:, indsthis);
        
        frmean = mean(frthis,2);
        frsem = lt_sem(frthis');
        if size(frthis,2)==1
            plot(x, frthis, 'Color', pcols_bos{bb});
        else
            shadedErrorBar(x, frmean, frsem, {'Color', pcols_bos{bb}},1);
        end
        lt_plot_zeroline;
        
        % ======== overlay song
        ons_mat = cell2mat(onsets(indsthis));
        if ~all(max(ons_mat)-min(ons_mat)<0.001)
            tmp = max(max(ons_mat)-min(ons_mat));
            lt_plot_annotation(1, ['some onsets up to: ' num2str(tmp) ' sec diff (across tirals)...'], 'm');
        end
        ons = median(ons_mat,1);
        
        offs_mat = cell2mat(offsets(indsthis));
        if ~all(max(offs_mat)-min(offs_mat)<0.001)
            tmp = max(max(offs_mat)-min(offs_mat));
            lt_plot_annotation(1, ['some oiffests up to: ' num2str(tmp) ' sec diff (across tirals)...'], 'm');
        end
        offs = median(offs_mat,1);
        
        lt_neural_QUICK_PlotSylPatches(ons, offs, median(frmean), 1);
        
        % ----------- PLOT LABELS
        labels = BOSlabels{bb};
        assert(length(labels)==length(ons), 'wierd');
        for k=1:length(labels)
            lt_plot_text(ons(k)+0.01, median(frmean)+5, labels(k), 'b');
        end
    end
    
    
    
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

linkaxes(hsplots, 'xy');
