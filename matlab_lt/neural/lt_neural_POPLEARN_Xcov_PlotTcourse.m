function lt_neural_POPLEARN_Xcov_PlotTcourse(OUTSTRUCT_XCOV, SwitchStruct, ...
    onlygoodexpt, PARAMS, dattype, lagwindows, clim, ffbinsedges_indstoplot, ...
    plotindivswitch, XLIM, YLIMGRAM)
%% lt 2/6/19 - takes lag window, and plots timecourse, relative to onset of syl

%%
if onlygoodexpt==1
    % ===== filter outstruct
    [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
end


%% group based on syl type

if strcmp(dattype, 'switch')
    fieldtoget = 'XcovgramBase';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_base] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'XcovgramWN';
    [~, ~, ~, ~, allbnum, allenum, allswnum, allDat_wn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
elseif strcmp(dattype, 'chan')
    fieldtoget = 'XcovgramBase';
    [allbnum, allenum, allswnum, allDat_base] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
    
    fieldtoget = 'XcovgramWN';
    [allbnum, allenum, allswnum, allDat_wn] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);
end
allDat_Diff = allDat_wn - allDat_base;


%% ==================== [PLOT] SUMMARY OVER ALL

figcount=1;
subplotrows=5;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG [BASE]');
ymat = squeeze(allDat_base(:,:,1,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');


% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG [WN]');
ymat = squeeze(allDat_wn(:,:,1,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');


% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG [WN - BASE]');
ymat = squeeze(allDat_Diff(:,:,1,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');


% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('SAME [WN - BASE]');
ymat = squeeze(allDat_Diff(:,:,2,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');

% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('DIFF [WN - BASE]');
ymat = squeeze(allDat_Diff(:,:,3,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');


% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG(WN - BASE)-SAME(WN-BASE)');
ymat = squeeze(allDat_Diff(:,:,1,:) - allDat_Diff(:,:,2,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');

% ==========
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('TARG(WN - BASE)-DIFF(WN-BASE)');
ymat = squeeze(allDat_Diff(:,:,1,:) - allDat_Diff(:,:,3,:));
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
xlim(XLIM); ylim(YLIMGRAM);

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');



% ###########################
linkaxes(hsplots, 'xy');

%% =================== [PLOT] ONE FOR EACH SWITCH
if plotindivswitch==1
    [indsgrp, indsgrpU] = lt_tools_grp2idx({allbnum, allenum, allswnum});
    
    
    figcount=1;
    subplotrows=5;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    
    for i=1:length(indsgrpU)
        
        indsthis = indsgrp == indsgrpU(i);
        hsplots = [];
        
        bnum = unique(allbnum(indsthis));
        enum = unique(allenum(indsthis));
        sw = unique(allswnum(indsthis));
        
        bname = SwitchStruct.bird(bnum).birdname;
        ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
        
        
        % ==========
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('TARG [BASE]');
        ylabel([bname '-' ename  '-sw' num2str(sw)]);
        ymat = squeeze(allDat_base(:,:,1,indsthis));
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
            ffbinsedges_indstoplot);
        
        % ==========
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('TARG [WN]');
        ymat = squeeze(allDat_wn(:,:,1,indsthis));
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
            ffbinsedges_indstoplot);
        
        % ==========
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('TARG [WN - BASE]');
        ylabel([bname '-' ename  '-sw' num2str(sw)]);
        ymat = squeeze(allDat_Diff(:,:,1,indsthis));
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
            ffbinsedges_indstoplot);
        
        % ==========
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('TARG(WN - BASE)-DIFF(WN-BASE)');
        ymat = squeeze(allDat_Diff(:,:,1,indsthis) - allDat_Diff(:,:,3,indsthis));
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
            ffbinsedges_indstoplot);
        
        
        % ===========
        linkaxes(hsplots, 'xy');
    end
    
end



