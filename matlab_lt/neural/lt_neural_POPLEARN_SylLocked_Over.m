function lt_neural_POPLEARN_SylLocked_Over(OUTSTRUCT, OUTSTRUCT_XCOV, SwitchStruct, ...
    SwitchCohStruct, MOTIFSTATS_Compiled, MOTIFSTATS_pop, SummaryStruct, ...
    PARAMS, onlygoodexpt, xtoplot, plotraw)
%% lt 2/8/19 - plots all syl locked LFP and smoothed MU

epochtoplot = 'base_epoch';
% fpass = [18 35]; % for bandpass filtering LFP.
fpass = [12 150]; % for bandpass filtering LFP.

clim = [-1 1]; % for plotting heat maps
XLIM = [-0.06 0.06]; % for lags.

normtype = 'unbiased';
% normtype = 'coeff';

%%
if onlygoodexpt==1
    % ===== filter outstruct
    expttype = 'xcov_spikes';
    [OUTSTRUCT] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct, expttype);
    [OUTSTRUCT_XCOV] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, SwitchStruct, expttype);
end


%% ============ go thru each switch

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.switch, ...
    OUTSTRUCT_XCOV.motifnum});

% ==== go thru each motif from each switch...
figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



% =================== COLLECT THINGS ACROSS TRIALS
DATSTRUCT.bnum = [];
DATSTRUCT.enum = [];
DATSTRUCT.switch = [];
DATSTRUCT.motifID = [];
DATSTRUCT.motifnum = [];
DATSTRUCT.FRsmooth_dat = {};
DATSTRUCT.FRsmooth_neurID = {};
DATSTRUCT.FRsmooth_bregions = {};
DATSTRUCT.FRsmooth_t = {};
DATSTRUCT.FRsmooth_chanlist = {};

DATSTRUCT.LFP_dat = {};
DATSTRUCT.LFP_chans = {};
DATSTRUCT.LFP_bregions = {};
DATSTRUCT.LFP_t = {};


    DATSTRUCT.sylonsets = {};
    DATSTRUCT.syloffsets = {};

    
% ================ COLLECT (AND PLOT)
for i=1:length(indsgrpU)
    disp(i);
    
    indlist= find(indsgrp == indsgrpU(i));
    
    bnum = unique(OUTSTRUCT_XCOV.bnum(indlist));
    enum = unique(OUTSTRUCT_XCOV.enum(indlist));
    ss = unique(OUTSTRUCT_XCOV.switch(indlist));
    mm = unique(OUTSTRUCT_XCOV.motifnum(indlist));
    neurlist = unique(OUTSTRUCT_XCOV.neurpair(indlist, :));
    inds_base = OUTSTRUCT_XCOV.inds_base_epoch{indlist(1)};
    inds_WN = OUTSTRUCT_XCOV.inds_WN_epoch{indlist(1)};
    
    bregionlist = {SummaryStruct.birds(bnum).neurons(neurlist).NOTE_Location};
    chanlist = [SummaryStruct.birds(bnum).neurons(neurlist).channel];
    
    bname = SummaryStruct.birds(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    motifname = unique(OUTSTRUCT_XCOV.motifname(indlist));
    
    motifID = lt_neural_QUICK_MotifID(bname, motifname);
    % ==== which trials to plot?
    if strcmp(epochtoplot, 'base_epoch')
        indstoplot = inds_base;
    elseif strcmp(epochtoplot, 'base_all')
        indstoplot = OUTSTRUCT_XCOV.inds_base_allgood{indlist(1)};
    else
        disp('WHICH TRIALS TO PLOT?')
        pause;
    end
    
    %% ############ GO THRU ALL NEURONS AND PLOT
    if plotraw==1
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title({[bname '-' ename '-sw' num2str(ss)], [motifname{1}]});
    end
    
    nset = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm).neursetused;
    DAT = MOTIFSTATS_pop.birds(bnum).exptnum(enum).DAT.setnum(nset).motif(mm).SegExtr_neurfakeID;
    segglobal = DAT(1).SegmentsExtract;
    
    frsmall = [];
    for nn=1:length(neurlist)
        assert(neurlist(nn)==DAT(nn).neurID_orig);
        seg = DAT(nn).SegmentsExtract;
        
        if strcmp(bregionlist{nn}, 'LMAN')
            pcol = [0.3 0.7 0.3];
        elseif strcmp(bregionlist{nn}, 'RA')
            pcol = [0.8 0.2 0.2];
        end
        
        % =========== extract smoothed fr
        seg = lt_neural_SmoothFR(seg, [], [], [], 0, segglobal);
        
        % =========== PLOT (SPIKES)
        %         lt_neural_PLOT_rasterline(seg.spk_Times, nn, pcol);
        
        % =========== PLOT (smoothed FR)
        x = seg(1).FRsmooth_xbin_CommonTrialDur - PARAMS.motif_predur;
        ymat = [seg(indstoplot).FRsmooth_rate_CommonTrialDur];
        indx = x>=xtoplot(1) & x<=xtoplot(2);
        x = x(indx);
        ymat = ymat(indx,:);
        ymean = mean(ymat,2);
        ysem = lt_sem(ymat');
        
        % -- zscore
        ysem = ysem./std(ymean);
        ymean = (ymean-mean(ymean))./std(ymean);
        if plotraw==1
            shadedErrorBar(x, ymean, ysem, {'Color', pcol},1);
        end
        
        % ======== ANNOTATE
        % --- neuron ID and chan
        chanthis = chanlist(nn);
        if plotraw==1
            lt_plot_text(x(end), ymean(end), ['n' num2str(neurlist(nn)) '-ch' num2str(chanlist(nn))], pcol, 8);
        end
        
        
        % ================ COLLECT ALL SMOOTHED FR
        frsmall = [frsmall; ymean'];
    end
    
    FRsmooth_dat = frsmall;
    FRsmooth_neurID = neurlist;
    FRsmooth_chanlist = chanlist;
    FRsmooth_bregions = bregionlist;
    FRsmooth_t = x;
    
    %         % ================ FORMAT
    %     axis tight;
    %     lt_plot_zeroline_vert;
    %     grid on
    
    
    
    %% ########### GO THRU ALL CHANS, AND PLOT LFP
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    
    
    DATLFP = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(ss).motifnum(mm);
    assert(size(DATLFP.lfpall,2)==length(DATLFP.lfpall_chans), 'just to make sure I rememebered format correctly');
    assert(size(DATLFP.lfpall,1)==length(DATLFP.tvals), 'just to make sure I rememebered format correctly');
    assert(length(segglobal) == size(DATLFP.lfpall,1), 'making sure that lfp and seg match up');
    assert(all([segglobal.song_datenum] == DATLFP.tvals), 'making sure that lfp and seg match up');
    
    % --------- GET UNIQUE CHANEL
    [~, indtmp] = unique(chanlist);
    chanlist_lfp = chanlist(indtmp);
    bregion_lfp = bregionlist(indtmp);
    
    lfpcollect = [];
    for cc=1:length(chanlist_lfp)
        chanthis = chanlist_lfp(cc);
        bregionthis = bregion_lfp{cc};
        
        if strcmp(bregionthis, 'LMAN')
            pcol = [0.3 0.7 0.3];
        elseif strcmp(bregionthis, 'RA')
            pcol = [0.8 0.2 0.2];
        end
        
        % ---- extract lfp data
        lfpall = DATLFP.lfpall(indstoplot, DATLFP.lfpall_chans==chanthis);
        lfpx = DATLFP.t_lfp;
        indx = lfpx>=xtoplot(1) & lfpx<=xtoplot(2);
        
        % -- keep the correct times.
        lfpall = cellfun(@(x)x(indx), lfpall, 'UniformOutput', 0);
        lfpx = lfpx(indx);
        
        % =========== get average over trials of interest.
        lfpall = squeeze(lt_neural_Coher_Cell2Mat(lfpall));
        
        
        if (0)
            %% ==== testing filtering of lfp
            lfpall2 = lt_neural_filter(lfpall, 1500, 0, fpass(1), fpass(2));
            
            figcount=1;
            subplotrows=4;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            for xx=1:8
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                indrand = randi(size(lfpall,2));
                
                plot(lfpx, lfpall(:,indrand), 'k');
                plot(lfpx, lfpall2(:,indrand), 'r');
                
            end
        end
        % ========== FILTER TO REMOVE LOW F FLULCTURATION
        lfpall = lt_neural_filter(lfpall, 1500, 0, fpass(1), fpass(2));
        
        ymean = mean(lfpall, 2);
        ysem = lt_sem(lfpall');
        
        % -- zscore
        ysem = ysem./std(ymean);
        ymean = (ymean-mean(ymean))./std(ymean);
        
        % ===== flip upside down
        ymeanplot = -ymean;
        
        % --------- SHIFT UP TO NOT OVERLAP NEURAL DATA
        ymeanplot = single(ymeanplot)+4;
        
        
        % ============= PLOT
        if plotraw==1
            shadedErrorBar(lfpx, ymeanplot, ysem, {'Color', pcol},1);
        end
        
        % ======== ANNOTATE
        % --- neuron ID and chan
        lt_plot_text(lfpx(end), ymeanplot(end), ['-ch' num2str(chanthis)], pcol, 8);
        
        % ============= COLLECT OVER CHANNELS
        lfpcollect = [lfpcollect; ymean'];
    end
    
    LFP_dat = lfpcollect;
    LFP_chans = chanlist_lfp;
    LFP_bregions = bregion_lfp;
    LFP_t = lfpx;
    
    
    
    
    
    %% ###################### OVERLAY SYL PATCHES
    ons = reshape([segglobal.motifsylOnsets], [], length(segglobal))';
    offs = reshape([segglobal.motifsylOffsets], [], length(segglobal))';
    
    ons = median(ons,1) - PARAMS.motif_predur;
    offs = median(offs,1) - PARAMS.motif_predur;
    
    %     ons(ons<xtoplot(1) | ons>xtoplot(2)) = [];
    %     offs(offs<xtoplot(1) | offs>xtoplot(2)) = [];
    if plotraw==1
        YLIM = ylim;
        YLIM(1) = YLIM(2)-1;
        lt_neural_QUICK_PlotSylPatches(ons, offs, YLIM);
        
        
        % ================ FORMAT
        ylabel('MU(lower), LFP(flipped, upper) [zscore]');
        axis tight;
        xlim(xtoplot);
        lt_plot_zeroline_vert;
        %     grid on
    end
    
    
    
    %% ============ COLLECT ACROSS CASES;
    
    DATSTRUCT.bnum = [DATSTRUCT.bnum; bnum];
    DATSTRUCT.enum = [DATSTRUCT.enum; enum];
    DATSTRUCT.switch = [DATSTRUCT.switch ; ss];
    DATSTRUCT.motifID = [DATSTRUCT.motifID; motifID];
    
    DATSTRUCT.sylonsets = [DATSTRUCT.sylonsets; ons];
    DATSTRUCT.syloffsets = [DATSTRUCT.syloffsets; offs];

    DATSTRUCT.motifnum = [DATSTRUCT.motifnum; mm];
    DATSTRUCT.FRsmooth_dat = [DATSTRUCT.FRsmooth_dat; FRsmooth_dat'];
    DATSTRUCT.FRsmooth_neurID = [DATSTRUCT.FRsmooth_neurID; FRsmooth_neurID];
    DATSTRUCT.FRsmooth_bregions = [DATSTRUCT.FRsmooth_bregions; {FRsmooth_bregions}];
    DATSTRUCT.FRsmooth_t = [DATSTRUCT.FRsmooth_t; FRsmooth_t];
    DATSTRUCT.FRsmooth_chanlist = [DATSTRUCT.FRsmooth_chanlist; FRsmooth_chanlist];
    
    DATSTRUCT.LFP_dat = [DATSTRUCT.LFP_dat; LFP_dat'];
    DATSTRUCT.LFP_chans = [DATSTRUCT.LFP_chans; LFP_chans];
    DATSTRUCT.LFP_bregions = [DATSTRUCT.LFP_bregions; {LFP_bregions}];
    DATSTRUCT.LFP_t = [DATSTRUCT.LFP_t; LFP_t'];
    
end


%% ==== [PLOT] FOR EACH MOTIF, PLOT ALL LFP AND FR ACROSS EXPERIMENTS...
% ONE PLOT FOR EACH UNQIUE MOTIF

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


[indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT.bnum, DATSTRUCT.motifID});

for i=1:length(indsgrpU)
    
    indsthis = indsgrp==indsgrpU(i);
    
    % ============= PLOT FR SMOOTH
    frmat = DATSTRUCT.FRsmooth_dat(indsthis);
    frmat = cell2mat(cellfun(@(x)x', frmat, 'UniformOutput', 0));
    neurlist = DATSTRUCT.FRsmooth_neurID(indsthis);
    neurlist = cell2mat(cellfun(@(x)x(:), neurlist, 'UniformOutput', 0));
    
    motifID = unique(DATSTRUCT.motifID(indsthis));
    bnum = unique(DATSTRUCT.bnum(indsthis));
    bname = SummaryStruct.birds(bnum).birdname;
    
    [~, motiflist_out] = ...
        lt_neural_QUICK_MotifID(bname);
    motifname_common = motiflist_out{motifID};
    
    
    bnum = unique(DATSTRUCT.bnum(indsthis));
    fr_bregion = {SummaryStruct.birds(bnum).neurons(neurlist).NOTE_Location};
    
    frx = DATSTRUCT.FRsmooth_t(indsthis);
    frx = frx{1};
    
    
    % ================= EXTRACT LFP
    lfpmat = DATSTRUCT.LFP_dat(indsthis);
    xlen = min(cellfun(@(x)size(x,1), lfpmat)); % since some are not same length...
    lfpmat = cellfun(@(x)x(1:xlen, :), lfpmat, 'UniformOutput', 0);
    lfpmat = cell2mat(cellfun(@(x)x', lfpmat, 'UniformOutput', 0));
    
    lfpx = DATSTRUCT.LFP_t(indsthis);
    lfpx = lfpx{1}(1:xlen);
    
    lfp_bregion = DATSTRUCT.LFP_bregions(indsthis);
    lfp_bregion = [lfp_bregion{:}]';
    
    
    % ############################## PLOT [LMAN]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    regionthis = 'LMAN';
    pcol = [0.3 0.7 0.3];
    
    title([bname  '-' motifname_common '-' regionthis]);
    
    % ========================== FR SMOOTH
    indsthis = strcmp(fr_bregion, regionthis);
    frthis = frmat(indsthis,:);
    
    plot(frx, frthis', 'Color', pcol);
    frmean = mean(frthis);
    frsem = lt_sem(frthis);
    %     if sum(indsthis)>1
    %         %         shadedErrorBar(frx, frmean, frsem, {'Color', pcol}, 1);
    %         %         plot(frx, frmean, 'Color', pcol, 'LineWidth', 2);
    %     end
    % ---- PLOT LOWER, FOR OVERLAY
    if sum(indsthis)>1
        shadedErrorBar(frx, frmean-4, frsem, {'Color', pcol}, 1);
    else
        plot(frx, frmean-4, 'Color', pcol, 'LineWidth', 2);
    end
    
    
    % ============================= LFP
    indsthis = strcmp(lfp_bregion, regionthis);
    lfpthis = lfpmat(indsthis,:);
    
    % --- flip sign, and add to mov eup plot
    lfpthis = -lfpthis + 4;
    
    plot(lfpx, lfpthis', 'Color', 'k');
    ymean = mean(lfpthis);
    ysem= lt_sem(lfpthis);
    %     if sum(indsthis)>1
    %         %         shadedErrorBar(lfpx, ymean, ysem, {'Color', 'k'}, 1);
    %         %         plot(lfpx, ymean, 'Color', 'k', 'LineWidth', 2);
    %     end
    % --- PLOT LOWER, FOR OVERLAY
    if sum(indsthis)>1
        shadedErrorBar(lfpx, ymean-8, ysem, {'Color', 'k'}, 1);
    else
        plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);
    end
    % plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);
    
    
        % ================ OVERLYA SYL PATCHES
    indsthis = indsgrp==indsgrpU(i);
    ons = DATSTRUCT.sylonsets(indsthis);
    offs = DATSTRUCT.syloffsets(indsthis);
    
    % -- only keep last ons at t=0 and one preceding
    tmp = cell2mat(cellfun(@(x)find(abs(x)<0.001), ons, 'UniformOutput', 0));
    
    onstmp = [];
    for j=1:length(ons)
        try
       onstmp = [onstmp; ons{j}(tmp(j)-1:tmp(j))];
        catch err
            onstmp = [onstmp; ons{j}(tmp(j))];
        end
    end
    
    offstmp = [];
    for j=1:length(ons)
        try
       offstmp = [offstmp; offs{j}(tmp(j)-1:tmp(j))];
        catch err
            offstmp = [offstmp; offs{j}(tmp(j))];
        end
    end
    
    ons = mean(onstmp,1);
    offs = mean(offstmp, 1);
    
    YLIM = ylim;
    lt_neural_QUICK_PlotSylPatches(ons, offs, [YLIM(2)-0.5 YLIM(2)]);

    
    % =========== format
    axis tight;
    xlim(xtoplot);
    lt_plot_zeroline;
    lt_plot_zeroline_vert
    
    
    
    
    % ############################## PLOT [LMAN]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    regionthis = 'RA';
    pcol = [0.8 0.2 0.2];
    
    title([bname  '-' motifname_common '-' regionthis]);
    
    % ========================== FR SMOOTH
    indsthis = strcmp(fr_bregion, regionthis);
    frthis = frmat(indsthis,:);
    
    plot(frx, frthis', 'Color', pcol);
    frmean = mean(frthis);
    frsem = lt_sem(frthis);
    %     if sum(indsthis)>1
    %         %         shadedErrorBar(frx, frmean, frsem, {'Color', pcol}, 1);
    %         %         plot(frx, frmean, 'Color', pcol, 'LineWidth', 2);
    %     end
    % ---- PLOT LOWER, FOR OVERLAY
    if sum(indsthis)>1
        shadedErrorBar(frx, frmean-4, frsem, {'Color', pcol}, 1);
    else
        plot(frx, frmean-4, 'Color', pcol, 'LineWidth', 2);
    end
    
    
    % ============================= LFP
    indsthis = strcmp(lfp_bregion, regionthis);
    lfpthis = lfpmat(indsthis,:);
    
    % --- flip sign, and add to mov eup plot
    lfpthis = -lfpthis + 4;
    
    plot(lfpx, lfpthis', 'Color', 'k');
    ymean = mean(lfpthis);
    ysem= lt_sem(lfpthis);
    %     if sum(indsthis)>1
    %         %         shadedErrorBar(lfpx, ymean, ysem, {'Color', 'k'}, 1);
    %         %         plot(lfpx, ymean, 'Color', 'k', 'LineWidth', 2);
    %     end
    % --- PLOT LOWER, FOR OVERLAY
    if sum(indsthis)>1
        shadedErrorBar(lfpx, ymean-8, ysem, {'Color', 'k'}, 1);
    else
        plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);
    end
    % plot(lfpx, ymean-8, 'Color', 'k', 'LineWidth', 2);
    
    
        % ================ OVERLYA SYL PATCHES
    indsthis = indsgrp==indsgrpU(i);
    ons = DATSTRUCT.sylonsets(indsthis);
    offs = DATSTRUCT.syloffsets(indsthis);
    
    % -- only keep last ons at t=0 and one preceding
    tmp = cell2mat(cellfun(@(x)find(abs(x)<0.001), ons, 'UniformOutput', 0));
    
    onstmp = [];
    for j=1:length(ons)
        try
       onstmp = [onstmp; ons{j}(tmp(j)-1:tmp(j))];
        catch err
            onstmp = [onstmp; ons{j}(tmp(j))];
        end
    end
    
    offstmp = [];
    for j=1:length(ons)
        try
       offstmp = [offstmp; offs{j}(tmp(j)-1:tmp(j))];
        catch err
            offstmp = [offstmp; offs{j}(tmp(j))];
        end
    end
    
    ons = mean(onstmp,1);
    offs = mean(offstmp, 1);
    
    YLIM = ylim;
    lt_neural_QUICK_PlotSylPatches(ons, offs, [YLIM(2)-0.5 YLIM(2)]);
    
    % =========== format
    axis tight;
%     xlim([xtoplot(1)-0.02 xtoplot(2)+0.02]);
    xlim(xtoplot);
    lt_plot_zeroline;
    lt_plot_zeroline_vert
end



%% =========== EXPAND DATSTRUCT to (1) LFP and (2) FR structures
% SO EACH ITEM IS ONE LFP TRACE (MEAN FOR A GIVEN SITE AND MOTIF)
tmp = DATSTRUCT.LFP_dat;
xlen = min(cellfun(@(x)size(x,1), tmp)); % since some are not same length...

DATSTRUCT_LFP.bnum = [];
DATSTRUCT_LFP.motifID = [];
DATSTRUCT_LFP.recsession = [];
DATSTRUCT_LFP.LFP_dat = [];
DATSTRUCT_LFP.LFP_chan = [];
DATSTRUCT_LFP.LFP_bregion = {};


DATSTRUCT_SPK.bnum = [];
DATSTRUCT_SPK.motifID = [];
DATSTRUCT_SPK.recsession = [];
DATSTRUCT_SPK.dat = [];
DATSTRUCT_SPK.chans = [];
DATSTRUCT_SPK.bregion = {};
DATSTRUCT_SPK.nID = [];


RecDays = lt_tools_grp2idx({DATSTRUCT.bnum, DATSTRUCT.enum, DATSTRUCT.switch});




for i=1:length(DATSTRUCT.bnum)
    
    motifID = DATSTRUCT.motifID(i);
    bnum = DATSTRUCT.bnum(i);
    recDay = RecDays(i);
    
    
    % ############################################## SPIKING
    frmat = DATSTRUCT.FRsmooth_dat{i}';
    frchans = DATSTRUCT.FRsmooth_chanlist{i}';
    frneurID= DATSTRUCT.FRsmooth_neurID{i};
    if size(frneurID,2)>1
        frneurID= frneurID';
    end
%     frneurID= DATSTRUCT.FRsmooth_
    frbregions = {SummaryStruct.birds(bnum).neurons(frneurID).NOTE_Location}';
    

    % ========== SAVE
    DATSTRUCT_SPK.bnum = [DATSTRUCT_SPK.bnum; bnum*ones(size(frmat,1),1)];
    DATSTRUCT_SPK.motifID = [DATSTRUCT_SPK.motifID; motifID*ones(size(frmat,1),1)];
    DATSTRUCT_SPK.recsession = [DATSTRUCT_SPK.recsession ; recDay*ones(size(frmat,1),1)];
    DATSTRUCT_SPK.dat = [DATSTRUCT_SPK.dat; frmat];
    DATSTRUCT_SPK.chans = [DATSTRUCT_SPK.chans; frchans];
    DATSTRUCT_SPK.bregion = [DATSTRUCT_SPK.bregion; frbregions];
    DATSTRUCT_SPK.nID = [DATSTRUCT_SPK.nID; frneurID];

    
    
    % ############################################## LFP
    lfpmat = DATSTRUCT.LFP_dat(i);
    
    
    % ================= EXPAND LFP TO ALL CHANS
    lfpmat = cellfun(@(x)x(1:xlen, :), lfpmat, 'UniformOutput', 0);
    lfpmat = cell2mat(cellfun(@(x)x', lfpmat, 'UniformOutput', 0));
    
    % ======= get bregion
    lfp_bregion = DATSTRUCT.LFP_bregions(i);
    lfp_bregion = [lfp_bregion{:}]';
    
    % === chan
    lfp_chans = DATSTRUCT.LFP_chans{i}';
    
    % ======================== SAVE
    DATSTRUCT_LFP.bnum = [DATSTRUCT_LFP.bnum; bnum*ones(size(lfp_chans))];
    DATSTRUCT_LFP.recsession = [DATSTRUCT_LFP.recsession; recDay*ones(size(lfp_chans))];
    DATSTRUCT_LFP.motifID = [DATSTRUCT_LFP.motifID; motifID*ones(size(lfp_chans))];
    DATSTRUCT_LFP.LFP_dat = [DATSTRUCT_LFP.LFP_dat; single(lfpmat)];
    DATSTRUCT_LFP.LFP_chan = [DATSTRUCT_LFP.LFP_chan; lfp_chans];
    DATSTRUCT_LFP.LFP_bregion = [DATSTRUCT_LFP.LFP_bregion; lfp_bregion];
    
end



Params.fr.x = DATSTRUCT.FRsmooth_t{1};
Params.lfp.x = DATSTRUCT.LFP_t{1}(1:length(DATSTRUCT_LFP.LFP_dat(1,:)));


%% ########## [SPIKE-LFP] GO THRU ALL PAIRS OF SPIKE-LFP AND COMPUTE XCOV
maxlag =  (xtoplot(2)-xtoplot(1))*(2/3); % seconds
binsize = 0.001; % defgault, for fr (1ms bins)

% =============== DOWNSAMPLE THE LFP DATA TO 1000HZ

% === prepare output
cc = length(DATSTRUCT_LFP.bnum)*length(DATSTRUCT_SPK.bnum);
DATSTRUCT_LfpSpk_Allpair.bnum = nan(cc,2);
DATSTRUCT_LfpSpk_Allpair.chans = nan(cc,2);
DATSTRUCT_LfpSpk_Allpair.recSessions = nan(cc,2);
DATSTRUCT_LfpSpk_Allpair.motifID = nan(cc,2);
DATSTRUCT_LfpSpk_Allpair.bregion = cell(cc,2);
DATSTRUCT_LfpSpk_Allpair.xcov = cell(cc,1);
DATSTRUCT_LfpSpk_Allpair.dattype = cell(cc,2);

cc = 0;
for i=1:length(DATSTRUCT_LFP.bnum)
    disp([num2str(i) '/' num2str(length(DATSTRUCT_LFP.bnum))]);
    for ii=1:length(DATSTRUCT_SPK.bnum)
        cc=cc+1;
        
        lfpthis = DATSTRUCT_LFP.LFP_dat(i,:);
        lfpx = Params.lfp.x;
        
        frthis = DATSTRUCT_SPK.dat(ii,:);
        frx = Params.fr.x;
        
        % ====== sanity chiec (interpolation looks correcrt)
        if (0)
            lt_figure; hold on;
            tmp = interp1(lfpx, lfpthis, frx);
            plot(lfpx, lfpthis, 'k');
            plot(frx, tmp, 'r')
        end
                
        % ===== interpolate/downsample the lfp (1500hz) to 1000hz(fr);
        lfpthis = interp1(lfpx, lfpthis, frx);
        
        % ===== make sure ther eis no nan
        indsnan = isnan(lfpthis);
        if sum(indsnan)>5
            keyboard
        end
        lfpthis(indsnan) = [];
        frthis(indsnan) = [];
        
        % ============= GET XCOV
        [xc, lags] = xcov(lfpthis, frthis, round(maxlag/binsize), normtype);
        
        % ============= GET STATS
        bb = [DATSTRUCT_LFP.bnum(i) DATSTRUCT_SPK.bnum(ii)];
        mm = [DATSTRUCT_LFP.motifID(i) DATSTRUCT_SPK.motifID(ii)];
        breg = {DATSTRUCT_LFP.LFP_bregion{i} DATSTRUCT_SPK.bregion{ii}};
        recsessions = [DATSTRUCT_LFP.recsession(i) DATSTRUCT_SPK.recsession(ii)];
        chans = [DATSTRUCT_LFP.LFP_chan(i) DATSTRUCT_SPK.chans(ii)];
        dattype = {'lfp', 'spk'};
        
        % =========== FLIP SO IT bregiosn are in alpha order
        [~, indsort] = sort(breg);
        bb = bb(indsort);
        mm = mm(indsort);
        breg = breg(indsort);
        recsessions = recsessions(indsort);
        chans = chans(indsort);
        dattype = dattype(indsort);
        if all(indsort==[1 2])
            %- do nothing
        elseif all(indsort==[2 1])
            % then flip
            xc = flipud(xc);
        else
            disp('HMMM, why');
            pause;
        end
        
        % ============ OUTPUT
        DATSTRUCT_LfpSpk_Allpair.recSessions(cc,:) = single(recsessions);
        DATSTRUCT_LfpSpk_Allpair.bnum(cc,:) = single(bb);
        DATSTRUCT_LfpSpk_Allpair.chans(cc,:) = single(chans);
        DATSTRUCT_LfpSpk_Allpair.motifID(cc,:) = single(mm);
        DATSTRUCT_LfpSpk_Allpair.bregion(cc,:) = breg;
        DATSTRUCT_LfpSpk_Allpair.xcov{cc} = single(xc');
        DATSTRUCT_LfpSpk_Allpair.dattype(cc,:) = dattype';
        
    end
end

assert(~any(cellfun(@(x)any(isnan(x)), DATSTRUCT_LfpSpk_Allpair.xcov)), ' wjhy some xcov give me nan?');

Params.LfpSpkPairs.xlags = lags*binsize;

% NOTE: should all be order of lfp, spk, except for cases across brain
% regions, which are ordered by brain region (alphabetically).


%% ################## [PLOTS] - LFP-SPK PAIRS



figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ####################### LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif/session/chan');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==DATSTRUCT_LfpSpk_Allpair.recSessions(:,2) ...% same session
    & DATSTRUCT_LfpSpk_Allpair.chans(:,1)==DATSTRUCT_LfpSpk_Allpair.chans(:,2); % same session


ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymat');
assert(length(indtmp)==size(ymat,1));
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', lags, 1:size(ymat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(ymat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif/session');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==DATSTRUCT_LfpSpk_Allpair.recSessions(:,2); % same session

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
% dattype1 = 'spk';
% dattype2 = 'lfp';
title('same bird/motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymat');
assert(length(indtmp)==size(ymat,1));
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', lags, 1:size(ymat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(ymat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
% dattype1 = 'spk';
% dattype2 = 'lfp';
title('same bird/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymat');
assert(length(indtmp)==size(ymat,1));
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', lags, 1:size(ymat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(ymat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
% dattype1 = 'spk';
% dattype2 = 'lfp';
title('diff bird');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)~=DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


% #################### RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif/session/chan');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==DATSTRUCT_LfpSpk_Allpair.recSessions(:,2) ...% same session
    & DATSTRUCT_LfpSpk_Allpair.chans(:,1)==DATSTRUCT_LfpSpk_Allpair.chans(:,2); % same session


ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif/session');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==DATSTRUCT_LfpSpk_Allpair.recSessions(:,2); % same session

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird,session/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==DATSTRUCT_LfpSpk_Allpair.recSessions(:,2); % same session

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
% dattype1 = 'spk';
% dattype2 = 'lfp';
title('same bird/motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymat');
assert(length(indtmp)==size(ymat,1));
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', lags, 1:size(ymat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(ymat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
% dattype1 = 'spk';
% dattype2 = 'lfp';
title('same bird/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymat');
assert(length(indtmp)==size(ymat,1));
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', lags, 1:size(ymat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(ymat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
% dattype1 = 'spk';
% dattype2 = 'lfp';
title('diff bird');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)~=DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
if strcmp(bregion1, bregion2)
    ylabel('xcov [LFP <-> SPK]');
else
    ylabel('xcov [LFP SPK order unclear (sorted by bregion)]');
end
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


% ============== FORMAT
linkaxes(hsplots, 'xy');
xlim(XLIM); ylim([-0.9 0.9]);
axis tight;

%% ========= [LFP - SPIKES, SUBTRACTING CONTROLS]
% i.e. each LFP-spike pair for same motif (motif, chan, session combo) picks 1 random case from
% different motif (same LFP, same bregion pairs) and subtracts that...



figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ===================================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same motif subtrac "control"[diff motif]');
ylabel('xcov(lfp vs. spk)');
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';

% NOTE: assuming bregions 1 and 2 are the same, and dattypes are lfp, spk.

% ===== go thru all lfp and for each one extract a matched control
indsthis = find(strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2)); % 

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
% --- collect one control (diff motif) for each case
ymat_control = nan(size(ymat));
for ii=1:length(indsthis)
    i = indsthis(ii);
    disp(ii);
    bnum = DATSTRUCT_LfpSpk_Allpair.bnum(i,1);
    mot = DATSTRUCT_LfpSpk_Allpair.motifID(i,1);
    chan1 = DATSTRUCT_LfpSpk_Allpair.chans(i,1);
    recsess = DATSTRUCT_LfpSpk_Allpair.recSessions(i,1);
    
    %^ ==== find cases that have different second motif.
    indstmp = find(DATSTRUCT_LfpSpk_Allpair.bnum(:,1) == bnum & DATSTRUCT_LfpSpk_Allpair.bnum(:,2) == bnum & ...
        DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==mot & DATSTRUCT_LfpSpk_Allpair.motifID(:,2)~=mot & ...
        DATSTRUCT_LfpSpk_Allpair.chans(:,1)==chan1 & ...
        DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==recsess & ...
        strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) ...
        & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2));
    
    % === pick a random indstmp to compare to
    ymat_control(ii,:) = DATSTRUCT_LfpSpk_Allpair.xcov{indstmp(randi(length(indstmp)))};
end
ymatdiff = ymat - ymat_control;
x = Params.LfpSpkPairs.xlags;

ymean = mean(ymatdiff);
ystd = std(ymatdiff);
ysem = lt_sem(ymatdiff);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;
ylim([-0.8 0.8]);
axis tight;

% -- heat maps
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymatdiff');
[~, indsort] = sort(indtmp);
ymatdiff = ymatdiff(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymatdiff', x, 1:size(ymatdiff,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;





% ===================================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('same motif subtrac "control"[diff motif]');
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';

% NOTE: assuming bregions 1 and 2 are the same, and dattypes are lfp, spk.

% ===== go thru all lfp and for each one extract a matched control
indsthis = find(strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2) ...
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2)); % 

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
% --- collect one control (diff motif) for each case
ymat_control = nan(size(ymat));
for ii=1:length(indsthis)
    i = indsthis(ii);
    disp(ii);
    bnum = DATSTRUCT_LfpSpk_Allpair.bnum(i,1);
    mot = DATSTRUCT_LfpSpk_Allpair.motifID(i,1);
    chan1 = DATSTRUCT_LfpSpk_Allpair.chans(i,1);
    recsess = DATSTRUCT_LfpSpk_Allpair.recSessions(i,1);
    
    %^ ==== find cases that have different second motif.
    indstmp = find(DATSTRUCT_LfpSpk_Allpair.bnum(:,1) == bnum & DATSTRUCT_LfpSpk_Allpair.bnum(:,2) == bnum & ...
        DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==mot & DATSTRUCT_LfpSpk_Allpair.motifID(:,2)~=mot & ...
        DATSTRUCT_LfpSpk_Allpair.chans(:,1)==chan1 & ...
        DATSTRUCT_LfpSpk_Allpair.recSessions(:,1)==recsess & ...
        strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) ...
        & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2));
    
    % === pick a random indstmp to compare to
    ymat_control(ii,:) = DATSTRUCT_LfpSpk_Allpair.xcov{indstmp(randi(length(indstmp)))};
end
ymatdiff = ymat - ymat_control;
x = Params.LfpSpkPairs.xlags;

ymean = mean(ymatdiff);
ystd = std(ymatdiff);
ysem = lt_sem(ymatdiff);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;
xlim(XLIM); ylim([-0.9 0.9]);
axis tight;

% -- heat maps
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(ymatdiff');
[~, indsort] = sort(indtmp);
ymatdiff = ymatdiff(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymatdiff', x, 1:size(ymatdiff,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% ============== FORMAT
linkaxes(hsplots, 'xy');


%% ========= across regions
% figcount=1;
% subplotrows=5;
% subplotcols=4;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];


% ======== LMAN SPK VS RA LFP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
dattype1 = 'spk';
dattype2 = 'lfp';
title('same bird/motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
ylabel(['xcov [' dattype1 '<->' dattype2]);

xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


% ======== LMAN SPK VS RA LFP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
dattype1 = 'spk';
dattype2 = 'lfp';
title('same bird/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
ylabel(['xcov [' dattype1 '<->' dattype2]);

xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


% ======== LMAN SPK VS RA LFP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
ylabel(['xcov [' dattype1 '<->' dattype2]);

xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;


% ======== LMAN SPK VS RA LFP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
ylabel(['xcov [' dattype1 '<->' dattype2]);

xlabel([bregion1 '-' bregion2]);
lt_plot_zeroline_vert;




% ============== FORMAT
linkaxes(hsplots, 'xy');
xlim(XLIM); ylim([-0.9 0.9]);
axis tight;



%% ================= [PLOT - LFP VS. SPIKES] - HEAT MAPS
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



% ==================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% ==================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% ==================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)==DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% ==================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('same bird/diff motif');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)==DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LfpSpk_Allpair.motifID(:,1)~=DATSTRUCT_LfpSpk_Allpair.motifID(:,2) ... % same motif
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
dattype1 = 'lfp';
dattype2 = 'spk';
title('diff bird');
indsthis = ...
    strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,1), bregion1) & strcmp(DATSTRUCT_LfpSpk_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_LfpSpk_Allpair.bnum(:,1)~=DATSTRUCT_LfpSpk_Allpair.bnum(:,2) ... % same bird
    & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,1), dattype1) & strcmp(DATSTRUCT_LfpSpk_Allpair.dattype(:,2), dattype2);

ymat = cell2mat(DATSTRUCT_LfpSpk_Allpair.xcov(indsthis));
x = Params.LfpSpkPairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);

xlabel([bregion1 '-' bregion2]);
ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;

%% ############## [LFP VS. MEAN SPIKING...]
maxlag =  (xtoplot(2)-xtoplot(1))*(2/3); % seconds
binsize = 0.001; % defgault, for fr (1ms bins)

% =============== DOWNSAMPLE THE LFP DATA TO 1000HZ

% === prepare output
% cc = length(DATSTRUCT_LFP.bnum)*length(DATSTRUCT_SPK.bnum);
% % DATSTRUCT_LfpSpk_Allpair.bnum = nan(cc,2);
% % DATSTRUCT_LfpSpk_Allpair.chans = nan(cc,2);
% % DATSTRUCT_LfpSpk_Allpair.recSessions = nan(cc,2);
% % DATSTRUCT_LfpSpk_Allpair.motifID = nan(cc,2);
% % DATSTRUCT_LfpSpk_Allpair.bregion = cell(cc,2);
% % DATSTRUCT_LfpSpk_Allpair.xcov = cell(cc,1);
% % DATSTRUCT_LfpSpk_Allpair.dattype = cell(cc,2);
%
% cc = 0;

Xcov_vsMeanFr_LMAN = [];
Xcov_vsMeanFr_RA = [];
Xcov_vsMeanFr_LMAN_CONTROL = {};
Xcov_vsMeanFr_RA_CONTROL = {};
maxmotifID = max(DATSTRUCT_LFP.motifID);
for i=1:length(DATSTRUCT_LFP.bnum)
    disp([num2str(i) '/' num2str(length(DATSTRUCT_LFP.bnum))]);
    
    % ######## COLLECT ALL SPIKING THAT IS SAME BIRD AND MOTIF
    bnum = DATSTRUCT_LFP.bnum(i);
    motifID = DATSTRUCT_LFP.motifID(i);
    lfpthis = DATSTRUCT_LFP.LFP_dat(i,:);
    lfpx = Params.lfp.x;
    frx = Params.fr.x;
    
    % ===== interpolate/downsample the lfp (1500hz) to 1000hz(fr);
    lfpthis = interp1(lfpx, lfpthis, frx);
    if sum(isnan(lfpthis))>5
        keyboard
    end
    
    
    % ================= LMAN SPIKES
    indsthis = DATSTRUCT_SPK.bnum==bnum & strcmp(DATSTRUCT_SPK.bregion, 'LMAN') ...
        & DATSTRUCT_SPK.motifID==motifID;
    
    frthis = mean(DATSTRUCT_SPK.dat(indsthis,:),1); % --- take mean over all chans
    
    % --- GET XCOV
    [xc_LMAN, lags] = xcov(lfpthis(~isnan(lfpthis)), frthis(~isnan(lfpthis)), ...
        round(maxlag/binsize), normtype);
    
    % --- CONTROL - COMPARE TO MEAN LMAN ACTIVITY ACROSS ALL OTHER MOTIFS [SAME BIRD]
    % ====== go thru all the other motifs
    xcall_LMAN = [];
    for mm=1:maxmotifID
        indsthis = DATSTRUCT_SPK.bnum==bnum & strcmp(DATSTRUCT_SPK.bregion, 'LMAN') ...
            & DATSTRUCT_SPK.motifID~=motifID & DATSTRUCT_SPK.motifID==mm;
        if ~any(indsthis)
            continue
        end
        
        frthis = mean(DATSTRUCT_SPK.dat(indsthis,:),1); % --- take mean over all chans

        % --- GET XCOV
        [xc_control, lags] = xcov(lfpthis(~isnan(lfpthis)), frthis(~isnan(lfpthis)), ...
            round(maxlag/binsize), normtype);
        xcall_LMAN = [xcall_LMAN; xc_control'];
    end
    
    
    
    
    
    
    % ================= RA SPIKES
    indsthis = DATSTRUCT_SPK.bnum==bnum & strcmp(DATSTRUCT_SPK.bregion, 'RA') ...
        & DATSTRUCT_SPK.motifID==motifID;
    
    frthis = mean(DATSTRUCT_SPK.dat(indsthis,:),1); % --- take mean over all chans
    
    % --- GET XCOV
    [xc_RA, lags] = xcov(lfpthis(~isnan(lfpthis)), frthis(~isnan(lfpthis)), ...
        round(maxlag/binsize), normtype);

    % ============ RA
    xcall_RA = [];
    for mm=1:maxmotifID
        indsthis = DATSTRUCT_SPK.bnum==bnum & strcmp(DATSTRUCT_SPK.bregion, 'RA') ...
            & DATSTRUCT_SPK.motifID~=motifID & DATSTRUCT_SPK.motifID==mm;
        if ~any(indsthis)
            continue
        end
        
        frthis = mean(DATSTRUCT_SPK.dat(indsthis,:),1); % --- take mean over all chans

        % --- GET XCOV
        [xc_control, lags] = xcov(lfpthis(~isnan(lfpthis)), frthis(~isnan(lfpthis)), ...
            round(maxlag/binsize), normtype);
        xcall_RA = [xcall_RA; xc_control'];
    end

    
    
    
%     % --- CONTROL - COMPARE TO MEAN RA ACTIVITY ACROSS ALL OTHER MOTIFS [SAME BIRD]
%     indsthis = DATSTRUCT_SPK.bnum==bnum & strcmp(DATSTRUCT_SPK.bregion, 'RA') ...
%         & DATSTRUCT_SPK.motifID~=motifID;
%     
%     frthis = mean(DATSTRUCT_SPK.dat(indsthis,:),1); % --- take mean over all chans
%     
%     % --- GET XCOV
%     [xc_RA_control, lags] = xcov(lfpthis(~isnan(lfpthis)), frthis(~isnan(lfpthis)), ...
%         round(maxlag/binsize), normtype);
%     
    
    
    % ================== OUTPUT
    Xcov_vsMeanFr_LMAN = [Xcov_vsMeanFr_LMAN; xc_LMAN'];
    Xcov_vsMeanFr_RA = [Xcov_vsMeanFr_RA; xc_RA'];
    Xcov_vsMeanFr_LMAN_CONTROL = [Xcov_vsMeanFr_LMAN_CONTROL; xcall_LMAN];
    Xcov_vsMeanFr_RA_CONTROL = [Xcov_vsMeanFr_RA_CONTROL; xcall_RA];
%     Xcov_vsMeanFr_LMAN = [Xcov_vsMeanFr_LMAN; xc_LMAN'];
%     Xcov_vsMeanFr_RA = [Xcov_vsMeanFr_RA; xc_RA'];
%     Xcov_vsMeanFr_LMAN_CONTROL = [Xcov_vsMeanFr_LMAN_CONTROL; xc_LMAN_control'];
%     Xcov_vsMeanFr_RA_CONTROL = [Xcov_vsMeanFr_RA_CONTROL; xc_RA_control'];
end
assert(size(Xcov_vsMeanFr_LMAN,1)==size(DATSTRUCT_LFP.bnum,1))
assert(size(Xcov_vsMeanFr_RA,1)==size(DATSTRUCT_LFP.bnum,1))

DATSTRUCT_LFP.Xcov_vsMeanFr_LMAN = Xcov_vsMeanFr_LMAN;
DATSTRUCT_LFP.Xcov_vsMeanFr_RA = Xcov_vsMeanFr_RA;
DATSTRUCT_LFP.Xcov_vsMeanFr_LMAN_CONTROL = Xcov_vsMeanFr_LMAN_CONTROL;
DATSTRUCT_LFP.Xcov_vsMeanFr_RA_CONTROL = Xcov_vsMeanFr_RA_CONTROL;

Params.lfp.vsFRmeanOverUnits.xlags = lags*binsize;

%% ============== [PLOTS - LFP VS. FR(MEAN ACROSS SITES)];

% [indsgrp, indsgrpU] = lt_tools_grp2idx({DATSTRUCT_LFP.bnum, DATSTRUCT_LFP.motifID});
%
% for i=1:length(indsgrpU)
%    indsthis = indsgrp==indsgrpU(i);
%
%    % === LMAN
%    xcovmat = DATSTRUCT_LFP.Xcov_vsMeanFr_LMAN(indsthis,:);
%
%    % === RA
% end

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ========== LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'LMAN';
bregionSPK = 'LMAN';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];


xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);
% -- plot distrubtions of times of peak
[maxvals, indtmp] = max(xcovmat');
lagtimes = lags(indtmp);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('lags of peaks');
xcenters = linspace(min(lagtimes), max(lagtimes), 20);
lt_plot_histogram(lagtimes, xcenters, 1, 0);
line([mean(lagtimes) mean(lagtimes)], ylim, 'Color', 'b');
% xlim([-1 1]);

% ========== LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'LMAN';
bregionSPK = 'LMAN';
title('samebird/diffmotif');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];

ymat = cell2mat(DATSTRUCT_LFP.(ftoget)(indsthis));
xcovmean = mean(ymat);
xcovstd = std(ymat);
xcovsem = lt_sem(ymat);
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = ymat;
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);

% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);
% -- plot distrubtions of times of peak
[maxvals, indtmp] = max(xcovmat');
lagtimes = lags(indtmp);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('lags of peaks');
xcenters = linspace(min(lagtimes), max(lagtimes), 20);
lt_plot_histogram(lagtimes, xcenters, 1, 0);
line([mean(lagtimes) mean(lagtimes)], ylim, 'Color', 'b');
% xlim([-1 1]);




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'LMAN';
bregionSPK = 'LMAN';
title('samemotif - diffmotif');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget1 = ['Xcov_vsMeanFr_' bregionSPK];
ftoget2 = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];

xcovmat = DATSTRUCT_LFP.(ftoget1)(indsthis,:) - cell2mat(cellfun(@(x)mean(x,1), DATSTRUCT_LFP.(ftoget2)(indsthis,:), 'UniformOutput', 0));
xcovmean = mean(xcovmat);
xcovstd = std(xcovmat);
xcovsem = lt_sem(xcovmat);
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);




% ========== LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'RA';
bregionSPK = 'RA';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);

lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);
% -- plot distrubtions of times of peak
lagtimes = lags(indtmp);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('lags of peaks');
xcenters = linspace(min(lagtimes), max(lagtimes), 20);
lt_plot_histogram(lagtimes, xcenters, 1, 0);
line([mean(lagtimes) mean(lagtimes)], ylim, 'Color', 'b');
% xlim([-1 1]);


% ========== LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'RA';
bregionSPK = 'RA';
title('samebird/diffmotif');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];

ymat = cell2mat(DATSTRUCT_LFP.(ftoget)(indsthis));
xcovmean = mean(ymat);
xcovstd = std(ymat);
xcovsem = lt_sem(ymat);
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = ymat;
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);

% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);
% -- plot distrubtions of times of peak
lagtimes = lags(indtmp);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('lags of peaks');
xcenters = linspace(min(lagtimes), max(lagtimes), 20);
lt_plot_histogram(lagtimes, xcenters, 1, 0);
line([mean(lagtimes) mean(lagtimes)], ylim, 'Color', 'b');
% xlim([-1 1]);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'RA';
bregionSPK = 'RA';
title('samemotif - diffmotif');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget1 = ['Xcov_vsMeanFr_' bregionSPK];
ftoget2 = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];

xcovmat = DATSTRUCT_LFP.(ftoget1)(indsthis,:) - cell2mat(cellfun(@(x)mean(x,1), DATSTRUCT_LFP.(ftoget2)(indsthis,:), 'UniformOutput', 0));
xcovmean = mean(xcovmat);
xcovstd = std(xcovmat);
xcovsem = lt_sem(xcovmat);
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);




% ========== LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'LMAN';
bregionSPK = 'RA';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);

lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'LMAN';
bregionSPK = 'RA';
title('samebird/diffmotif');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];

ymat = cell2mat(DATSTRUCT_LFP.(ftoget)(indsthis));
xcovmean = mean(ymat);
xcovstd = std(ymat);
xcovsem = lt_sem(ymat);
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = ymat;
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);

% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);




% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% hsplots = [hsplots; hsplot];
% bregionLFP = 'LMAN';
% bregionSPK = 'RA';
% title('samemotif - diffmotif');
% ylabel('xcov(coeff)');
% 
% xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
% indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
% ftoget1 = ['Xcov_vsMeanFr_' bregionSPK];
% ftoget2 = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];
% 
% xcovmat = DATSTRUCT_LFP.(ftoget1)(indsthis,:) - DATSTRUCT_LFP.(ftoget2)(indsthis,:);
% xcovmean = mean(xcovmat);
% xcovstd = std(xcovmat);
% xcovsem = lt_sem(xcovmat);
% lags = Params.lfp.vsFRmeanOverUnits.xlags;
% 
% shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
% shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
% axis tight;
% lt_plot_zeroline_vert;
% 
% % --- plot heat map
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% [~, indtmp] = max(xcovmat');
% assert(length(indtmp)==size(xcovmat,1));
% [~, indsort] = sort(indtmp);
% xcovmat = xcovmat(indsort,:); % sort in order of peak.
% 
% lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% % --- plot distribution of peak magnitudes
% [maxvals, indtmp] = max(xcovmat');
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('max corr vals');
% xcenters = -0.95:0.1:0.95;
% lt_plot_histogram(maxvals, xcenters, 1, 0);
% line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
% xlim([-1 1]);



% ========== LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'RA';
bregionSPK = 'LMAN';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
ylim([-0.8 0.8]);


lt_plot_zeroline_vert;
% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregionLFP = 'RA';
bregionSPK = 'LMAN';
title('samebird/diffmotif');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];

ymat = cell2mat(DATSTRUCT_LFP.(ftoget)(indsthis));
xcovmean = mean(ymat);
xcovstd = std(ymat);
xcovsem = lt_sem(ymat);
lags = Params.lfp.vsFRmeanOverUnits.xlags;

shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
axis tight;
lt_plot_zeroline_vert;

% --- plot heat map
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xcovmat = ymat;
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);

% --- plot distribution of peak magnitudes
[maxvals, indtmp] = max(xcovmat');
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
xlabel('max corr vals');
xcenters = -0.95:0.1:0.95;
lt_plot_histogram(maxvals, xcenters, 1, 0);
line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
xlim([-1 1]);




% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% hsplots = [hsplots; hsplot];
% bregionLFP = 'RA';
% bregionSPK = 'LMAN';
% title('samemotif - diffmotif');
% ylabel('xcov(coeff)');
% 
% xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
% indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
% ftoget1 = ['Xcov_vsMeanFr_' bregionSPK];
% ftoget2 = ['Xcov_vsMeanFr_' bregionSPK '_CONTROL'];
% 
% xcovmat = DATSTRUCT_LFP.(ftoget1)(indsthis,:) - DATSTRUCT_LFP.(ftoget2)(indsthis,:);
% xcovmean = mean(xcovmat);
% xcovstd = std(xcovmat);
% xcovsem = lt_sem(xcovmat);
% lags = Params.lfp.vsFRmeanOverUnits.xlags;
% 
% shadedErrorBar(lags, xcovmean, xcovstd, {'Color', [0.5 0.5 0.5]}, 1);
% shadedErrorBar(lags, xcovmean, xcovsem, {'Color', [0.5 0.5 0.5]}, 1);
% axis tight;
% lt_plot_zeroline_vert;
% 
% % --- plot heat map
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% [~, indtmp] = max(xcovmat');
% assert(length(indtmp)==size(xcovmat,1));
% [~, indsort] = sort(indtmp);
% xcovmat = xcovmat(indsort,:); % sort in order of peak.
% 
% lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);
% % --- plot distribution of peak magnitudes
% [maxvals, indtmp] = max(xcovmat');
% [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% xlabel('max corr vals');
% xcenters = -0.95:0.1:0.95;
% lt_plot_histogram(maxvals, xcenters, 1, 0);
% line([mean(maxvals) mean(maxvals)], ylim, 'Color', 'b');
% xlim([-1 1]);



% ======== format
linkaxes(hsplots, 'xy');


%% ====== [PLOTS - HEAT MAPS OF LFP VS. MEAN UNITS]
if (0) % IGNROE, subsumed by above
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


% ============================= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];

bregionLFP = 'LMAN';
bregionSPK = 'LMAN';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
% 
% xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

% -- sort by time of peak
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);



% ============================= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];

bregionLFP = 'RA';
bregionSPK = 'RA';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
% 
% xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

% -- sort by time of peak
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);



% ============================= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];

bregionLFP = 'LMAN';
bregionSPK = 'RA';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
% 
% xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

% -- sort by time of peak
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);



% ============================= 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];

bregionLFP = 'RA';
bregionSPK = 'LMAN';
title('same motif[singleLFP vs. meanFR');
ylabel('xcov(coeff)');

xlabel(['LFP(' bregionLFP  ') <--> SPK(mean,units)[' bregionSPK ']']);
indsthis = strcmp(DATSTRUCT_LFP.LFP_bregion, bregionLFP);
ftoget = ['Xcov_vsMeanFr_' bregionSPK];

xcovmat = DATSTRUCT_LFP.(ftoget)(indsthis,:);
% 
% xcovmean = mean(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovstd = std(DATSTRUCT_LFP.(ftoget)(indsthis,:));
% xcovsem = lt_sem(DATSTRUCT_LFP.(ftoget)(indsthis,:));
lags = Params.lfp.vsFRmeanOverUnits.xlags;

% -- sort by time of peak
[~, indtmp] = max(xcovmat');
assert(length(indtmp)==size(xcovmat,1));
[~, indsort] = sort(indtmp);
xcovmat = xcovmat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(xcovmat', lags, 1:size(xcovmat,1), 1, [], clim);




end

%% ###################################################
%% ####### [SPIKE-SPIKE XCOV]

%% ========= [LFP XCOV] - GO THRU ALL PAIRS OF LFP AND COMPUTE XCOV
binsize = 0.001;
maxlag =  (xtoplot(2)-xtoplot(1))*(2/3); % seconds

% ======= how many rows expected?
cc=0;
for i=1:length(DATSTRUCT_SPK.bnum)
    for ii=i+1:length(DATSTRUCT_SPK.bnum)
        cc = cc+1;
    end
end


DATSTRUCT_SPK_Allpair.bnum = nan(cc,2);
% DATSTRUCT_SPK_Allpair.chans = nan(cc,2);
DATSTRUCT_SPK_Allpair.chans = nan(cc,2);
DATSTRUCT_SPK_Allpair.neurID = nan(cc,2);
DATSTRUCT_SPK_Allpair.recSessions = nan(cc,2);
DATSTRUCT_SPK_Allpair.motifID = nan(cc,2);
DATSTRUCT_SPK_Allpair.bregion = cell(cc,2);
DATSTRUCT_SPK_Allpair.xcov = cell(cc,1);
cc=0;
for i=1:length(DATSTRUCT_SPK.bnum)
    disp(i);
    for ii=i+1:length(DATSTRUCT_SPK.bnum)
        cc = cc+1;
        
        spk1 = DATSTRUCT_SPK.dat(i,:);
        spk2 = DATSTRUCT_SPK.dat(ii,:);
        
        [xc, lags] = xcov(spk1', spk2', round(maxlag/binsize), normtype);
        %         [xc, lags] = xcov(lfp1', lfp2', round(maxlag/lfpbinsize), 'unbiased');
        
        % ==== feature about this pair
        bb = [DATSTRUCT_SPK.bnum(i) DATSTRUCT_SPK.bnum(ii)];
        mm = [DATSTRUCT_SPK.motifID(i) DATSTRUCT_SPK.motifID(ii)];
        breg = {DATSTRUCT_SPK.bregion{i} DATSTRUCT_SPK.bregion{ii}};
        recsessions = [DATSTRUCT_SPK.recsession(i) DATSTRUCT_SPK.recsession(ii)];
        chans = [DATSTRUCT_SPK.chans(i) DATSTRUCT_SPK.chans(ii)];
        neur = [DATSTRUCT_SPK.nID(i) DATSTRUCT_SPK.nID(ii)];
        
        % =========== FLIP SO IT bregiosn are in alpha order
        [~, indsort] = sort(breg);
        bb = bb(indsort);
        mm = mm(indsort);
        breg = breg(indsort);
        recsessions = recsessions(indsort);
        chans = chans(indsort);
        neur = neur(indsort);
        if all(indsort==[1 2])
            %- do nothing
        elseif all(indsort==[2 1])
            % then flip
            xc = flipud(xc);
        else
            disp('HMMM, why');
            pause;
        end
        
        % ============ OUTPUT
        DATSTRUCT_SPK_Allpair.recSessions(cc,:) = recsessions;
        DATSTRUCT_SPK_Allpair.bnum(cc,:) = bb;
        DATSTRUCT_SPK_Allpair.chans(cc,:) = chans;
        DATSTRUCT_SPK_Allpair.neurID(cc,:) = neur;
        DATSTRUCT_SPK_Allpair.motifID(cc,:) = mm;
        DATSTRUCT_SPK_Allpair.bregion(cc,:) = breg;
        DATSTRUCT_SPK_Allpair.xcov{cc} = xc';
        
    end
end

Params.spkpairs.xlags = lags*binsize;

%
% %% ========= [XCORRELATION BETWEEN ALL PAIRS OF TRIALS]
%
% lfpmatall = DATSTRUCT_LFP.LFP_dat;
% lfpbinsize = lfpx(2)-lfpx(1);
% maxlag = 0.08; % seconds
%
% % maxlag =
% [xcovmat, lags] = xcov(lfpmatall', round(maxlag/lfpbinsize), normtype);


assert(~any(cellfun(@(x)any(isnan(x)), DATSTRUCT_SPK_Allpair.xcov)), ' wjhy some xcov give me nan?');



%% =========== [PLOT] - SPIKE PAIRS.

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird/motif/session');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.recSessions(:,1)==DATSTRUCT_SPK_Allpair.recSessions(:,2); %


ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)[spk-spk]');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird/session, diff motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.recSessions(:,1)==DATSTRUCT_SPK_Allpair.recSessions(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird/same motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird/diff motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;







% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/motif/session');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.recSessions(:,1)==DATSTRUCT_SPK_Allpair.recSessions(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/session, diff motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.recSessions(:,1)==DATSTRUCT_SPK_Allpair.recSessions(:,2); %


ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/same motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/diff motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;





[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('same bird/motif/session');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.recSessions(:,1)==DATSTRUCT_SPK_Allpair.recSessions(:,2); %


ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('same bird/session, diff motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.recSessions(:,1)==DATSTRUCT_SPK_Allpair.recSessions(:,2); %


ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('same bird/same motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('same bird/diff motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2); %
ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
xlim(XLIM); ylim([-0.9 0.9]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;





% =========== format plots
linkaxes(hsplots, 'xy');



%% ========= [PLOTS] - MORE, for spike-spike


figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird/motif, diff neuron');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)[spk-spk]');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird/neuron, diffmotif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)==DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('same bird, diffneuron/motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)~=DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'LMAN';
title('diff bird');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)~=DATSTRUCT_SPK_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;







% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/motif, diff neuron');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)[spk-spk]');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/neuron, diffmotif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)==DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird, diffneuron/motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)~=DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('diff bird');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)~=DATSTRUCT_SPK_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/motif, diff neuron');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)[spk-spk]');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird/neuron, diffmotif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)==DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('same bird, diffneuron/motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)~=DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'RA';
bregion2 = 'RA';
title('diff bird');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)~=DATSTRUCT_SPK_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;










% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('same bird/motif, diff neuron');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)==DATSTRUCT_SPK_Allpair.motifID(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)[spk-spk]');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('same bird, diffneuron/motif');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)==DATSTRUCT_SPK_Allpair.bnum(:,2) ... %
    & DATSTRUCT_SPK_Allpair.motifID(:,1)~=DATSTRUCT_SPK_Allpair.motifID(:,2) ...
    & DATSTRUCT_SPK_Allpair.neurID(:,1)~=DATSTRUCT_SPK_Allpair.neurID(:,2); 

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
axis tight;
ylim([-0.8 0.8]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% #################################### 
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
bregion1 = 'LMAN';
bregion2 = 'RA';
title('diff bird');
indsthis = strcmp(DATSTRUCT_SPK_Allpair.bregion(:,1), bregion1) & ...
    strcmp(DATSTRUCT_SPK_Allpair.bregion(:,2), bregion2) & ...
    DATSTRUCT_SPK_Allpair.bnum(:,1)~=DATSTRUCT_SPK_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_SPK_Allpair.xcov(indsthis));
x = Params.spkpairs.xlags;
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);

xlabel([bregion1 '-' bregion2]);
ylabel('xcov (std)');
xlim(XLIM); ylim([-0.9 0.9]);
lt_plot_zeroline_vert;

% ============= PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');

xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% ==============
linkaxes(hsplots, 'xy');


%% ################################################################
%% ========= [LFP XCOV] - GO THRU ALL PAIRS OF LFP AND COMPUTE XCOV
lfpbinsize = Params.lfp.x(2)-Params.lfp.x(1);
maxlag =  (xtoplot(2)-xtoplot(1))*(2/3); % seconds

% ======= how many rows expected?
cc=0;
for i=1:length(DATSTRUCT_LFP.bnum)
    for ii=i+1:length(DATSTRUCT_LFP.bnum)
        cc = cc+1;
    end
end


DATSTRUCT_LFP_Allpair.bnum = nan(cc,2);
DATSTRUCT_LFP_Allpair.chans = nan(cc,2);
DATSTRUCT_LFP_Allpair.recSessions = nan(cc,2);
DATSTRUCT_LFP_Allpair.motifID = nan(cc,2);
DATSTRUCT_LFP_Allpair.bregion = cell(cc,2);
DATSTRUCT_LFP_Allpair.xcov = cell(cc,1);
cc=0;
for i=1:length(DATSTRUCT_LFP.bnum)
    disp(i);
    for ii=i+1:length(DATSTRUCT_LFP.bnum)
        cc = cc+1;
        
        lfp1 = DATSTRUCT_LFP.LFP_dat(i,:);
        lfp2 = DATSTRUCT_LFP.LFP_dat(ii,:);
        
        [xc, lags] = xcov(lfp1', lfp2', round(maxlag/lfpbinsize), normtype);
        %         [xc, lags] = xcov(lfp1', lfp2', round(maxlag/lfpbinsize), 'unbiased');
        
        % ==== feature about this pair
        bb = [DATSTRUCT_LFP.bnum(i) DATSTRUCT_LFP.bnum(ii)];
        mm = [DATSTRUCT_LFP.motifID(i) DATSTRUCT_LFP.motifID(ii)];
        breg = {DATSTRUCT_LFP.LFP_bregion{i} DATSTRUCT_LFP.LFP_bregion{ii}};
        recsessions = [DATSTRUCT_LFP.recsession(i) DATSTRUCT_LFP.recsession(ii)];
        chans = [DATSTRUCT_LFP.LFP_chan(i) DATSTRUCT_LFP.LFP_chan(ii)];
        
        % =========== FLIP SO IT bregiosn are in alpha order
        [~, indsort] = sort(breg);
        bb = bb(indsort);
        mm = mm(indsort);
        breg = breg(indsort);
        recsessions = recsessions(indsort);
        chans = chans(indsort);
        if all(indsort==[1 2])
            %- do nothing
        elseif all(indsort==[2 1])
            % then flip
            xc = flipud(xc);
        else
            disp('HMMM, why');
            pause;
        end
        
        % ============ OUTPUT
        DATSTRUCT_LFP_Allpair.recSessions(cc,:) = recsessions;
        DATSTRUCT_LFP_Allpair.bnum(cc,:) = bb;
        DATSTRUCT_LFP_Allpair.chans(cc,:) = chans;
        DATSTRUCT_LFP_Allpair.motifID(cc,:) = mm;
        DATSTRUCT_LFP_Allpair.bregion(cc,:) = breg;
        DATSTRUCT_LFP_Allpair.xcov{cc} = xc';
        
    end
end

Params.lfppairs.xlags = lags*lfpbinsize;

%
% %% ========= [XCORRELATION BETWEEN ALL PAIRS OF TRIALS]
%
% lfpmatall = DATSTRUCT_LFP.LFP_dat;
% lfpbinsize = lfpx(2)-lfpx(1);
% maxlag = 0.08; % seconds
%
% % maxlag =
% [xcovmat, lags] = xcov(lfpmatall', round(maxlag/lfpbinsize), normtype);


assert(~any(cellfun(@(x)any(isnan(x)), DATSTRUCT_LFP_Allpair.xcov)), ' wjhy some xcov give me nan?');

%% =========== [PLOT] XCOV


figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ############## LMAN VS. RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ystd = std(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/diff motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2); % diff motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('diff bird');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)~=DATSTRUCT_LFP_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% ############################################## LMAN VS. LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN  [diff recordings]');
ylabel('xcov (std)');
title('same bird/same motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN  [diff recordings]');
ylabel('xcov (std)');
title('same bird/diff motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN  [diff recordings]');
ylabel('xcov (std)');
title('diff bird');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)~=DATSTRUCT_LFP_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




% ############## RA VS. RA (DIFF MOTIF)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA  [diff recordings]');
ylabel('xcov (std)');
title('same bird/same motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA  [diff recordings]');
ylabel('xcov (std)');
title('same bird/diff motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2); %

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA  [diff recordings]');
ylabel('xcov (std)');
title('diff bird');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)~=DATSTRUCT_LFP_Allpair.bnum(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
xlim(XLIM); ylim([-0.9 0.9]);

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% ==============
linkaxes(hsplots, 'xy');

%% ############################### SAME VS. DIFFERENT REC SESSION
figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/same motif/same session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/session, diffmotif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/same motif/DIFF session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)~=DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/diff motif/diff session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)~=DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/same motif/same session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;

% ================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA -- RA');
ylabel('xcov (std)');
title('same bird/session, diffmotif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/same motif/DIFF session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)~=DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/diff motif/diff session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)~=DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/same motif/same session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


% ================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN -- RA');
ylabel('xcov (std)');
title('same bird/session, diffmotif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/same motif/DIFF session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)~=DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/diff motif/diff session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)~=DATSTRUCT_LFP_Allpair.recSessions(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline_vert;
xlim(XLIM); ylim([-0.9 0.9]);

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% ==============
linkaxes(hsplots, 'xy');

%% %  ########## DIFF CHANS (SAME MOTIF) VS. SAME CHAN (DIFF MOTIFS)

figcount=1;
subplotrows=6;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same motif/same session/diff chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)~=DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('diff motif/same session/same chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)==DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('diff motif/same session/diff chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)~=DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;


%  ########## DIFF CHANS (SAME MOTIF) VS. SAME CHAN (DIFF MOTIFS)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same motif/same session/diff chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)~=DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('diff motif/same session/same chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)==DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA -- RA');
ylabel('xcov (std)');
title('diff motif/same session/diff chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)~=DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;
axis tight;
ylim([-0.8 0.8]);

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;






%  ########## DIFF CHANS (SAME MOTIF) VS. SAME CHAN (DIFF MOTIFS)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN -- RA ');
ylabel('xcov (std)');
title('same motif/same session/diff chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)~=DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;





[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN -- RA');
ylabel('xcov (std)');
title('diff motif/same session/diff chan');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... %
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2) ...
    & DATSTRUCT_LFP_Allpair.chans(:,1)~=DATSTRUCT_LFP_Allpair.chans(:,2);

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;
% plot(x, ymat', '-k');
ymean = mean(ymat);
ysem = lt_sem(ymat);
shadedErrorBar(x, ymean, ystd, {'Color', [0.5 0.5 0.5]}, 1);
shadedErrorBar(x, ymean, ysem, {'Color', 'k'}, 1);
n = size(ymat,1);
lt_plot_annotation(1, ['N=' num2str(n)], 'b');
lt_plot_zeroline;
xlim(XLIM); ylim([-0.9 0.9]);

% --------------- PLOT HEAT MAP
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.
lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim, '', '', '', 'East');
% xlabel([bregion1 '-' bregion2]);
% ylabel(['xcov [' dattype1 '<->' dattype2]);
lt_plot_zeroline_vert;



% ==============
linkaxes(hsplots, 'xy');




%% ========== [RELATIVE TIMING OF LMAN AND RA LFP]
if (0) % SINCE ALREADY PLOTTING ABOVE, ADJACENT TO MEANS.
figcount=1;
subplotrows=3;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/diffmotif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/motif/session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...% same motif
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- LMAN ');
ylabel('xcov (std)');
title('same bird/session, diff motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'LMAN') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...% same motif
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);


% ===========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/diffmotif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/motif/session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...% same motif
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('RA  -- RA ');
ylabel('xcov (std)');
title('same bird/session, diff motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'RA') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...% same motif
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);


%$ =========================
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/diffmotif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/motif/session');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)==DATSTRUCT_LFP_Allpair.motifID(:,2) ...% same motif
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
xlabel('LMAN  -- RA ');
ylabel('xcov (std)');
title('same bird/session, diff motif');
indsthis = strcmp(DATSTRUCT_LFP_Allpair.bregion(:,1), 'LMAN') & ...
    strcmp(DATSTRUCT_LFP_Allpair.bregion(:,2), 'RA') & ...
    DATSTRUCT_LFP_Allpair.bnum(:,1)==DATSTRUCT_LFP_Allpair.bnum(:,2) ... % same bird
    & DATSTRUCT_LFP_Allpair.motifID(:,1)~=DATSTRUCT_LFP_Allpair.motifID(:,2) ...% same motif
    & DATSTRUCT_LFP_Allpair.recSessions(:,1)==DATSTRUCT_LFP_Allpair.recSessions(:,2); % same motif

ymat = cell2mat(DATSTRUCT_LFP_Allpair.xcov(indsthis));
x = Params.lfppairs.xlags;

% -- sort by time of peak
[~, indtmp] = max(ymat');
[~, indsort] = sort(indtmp);
ymat = ymat(indsort,:); % sort in order of peak.

lt_neural_Coher_Plot(ymat', x, 1:size(ymat,1), 1, [], clim);
end


