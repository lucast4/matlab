function lt_neural_POPLEARN_Xcov_PlotTcourse(OUTSTRUCT_XCOV, SwitchStruct, ...
    onlygoodexpt, PARAMS, dattype, lagwindows, clim, ffbinsedges_indstoplot, ...
    plotindivswitch, XLIM, YLIMGRAM, plotlines, ffsplitparams, ...
    doOnlyPlotAdjacentSyls)
%% lt 2/6/19 - takes lag window, and plots timecourse, relative to onset of syl

%%
if onlygoodexpt==1
    % ===== filter outstruct
    [OUTSTRUCT_XCOV, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
end

%% ===== ONLY KEEP SYLLABLES ADJACENT TO TARGET
% - finds each targ syl (for a given channel pair). then iterates over all
% otheres to find the preceding syl. Then replaces outstruct so that it
% only has the adjacent syl, and those are named artifically as target.

NUMBACK = 1;
if doOnlyPlotAdjacentSyls==1
    % --- go thru each expt and find the one syl that is next to target
    inds = find(OUTSTRUCT_XCOV.istarg==1);
    indstokeep = [];
    for i=inds'

    % --- find what bird/syl this corresponds to
    nn = OUTSTRUCT_XCOV.neurpairnum(i);
    bb = OUTSTRUCT_XCOV.bnum(i);
    ee = OUTSTRUCT_XCOV.enum(i);
    ss = OUTSTRUCT_XCOV.swnum(i);
    mname = OUTSTRUCT_XCOV.motifname{i};

    bname = SwitchStruct.bird(bb).birdname;

    [~, ~, motifposition_targ] = ...
    lt_neural_QUICK_MotifID(bname, mname);

    % -- find the adjacent syllable
    indstocheck = find(OUTSTRUCT_XCOV.bnum==bb & OUTSTRUCT_XCOV.neurpairnum==nn & ...
        OUTSTRUCT_XCOV.enum==ee & OUTSTRUCT_XCOV.switch==ss);

    disp('---');
    disp([bname '-' num2str(ee) '-' num2str(ss) '-n' num2str(nn)])
    %     disp(motifposition_targ)
    for ii=indstocheck'

        % --- get the posotion of this syl
        mname_2 = OUTSTRUCT_XCOV.motifname{ii};

        [~, ~, motifposition_other] = lt_neural_QUICK_MotifID(bname, mname_2);
    %         disp(motifposition_other)
        if motifposition_other(1)==motifposition_targ(1) & ...
                motifposition_other(2)==(motifposition_targ(2)-NUMBACK)
                % --- if it is preceding the target, then keep
            disp([num2str(ii) ', ' mname_2 ' before ' mname]);    
            indstokeep = [indstokeep ii];
        end


    end

    end

    % =============== ONLY KEEP THESE INDS - RENAME EVERYTHING AS TARGET
    OUTSTRUCT_XCOV = lt_structure_subsample_all_fields(OUTSTRUCT_XCOV, indstokeep, 1);
    OUTSTRUCT_XCOV.istarg = ones(length(OUTSTRUCT_XCOV.istarg),1);
    OUTSTRUCT_XCOV.issame = ones(length(OUTSTRUCT_XCOV.issame),1);
end

%% ============== IF WANT TO PLOT ADAPTIVE PITCH TRIALS (OR NONADAPTIVE)
if ffsplitparams.dosplit==1
    epochtoplot = ffsplitparams.epochtoplot;
    adaptivebin = ffsplitparams.adaptivebin;
    
    % ===== COLLECT ALL FFSPLIT DATA
    ybase = OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Base;
    ywn = cellfun(@(x)nanmean(x(:,:, epochtoplot),3), OUTSTRUCT_XCOV.XcovgramWN_FFsplits_Epochs, 'UniformOutput', 0);
    
    % --- flip those are driving learning down
    indstoflip = OUTSTRUCT_XCOV.learndirTarg==-1;
    ybase(indstoflip, :) = fliplr(ybase(indstoflip, :));
    ywn(indstoflip, :) = fliplr(ywn(indstoflip, :));
    
    if adaptivebin<3
        % --- pull out the desired bin
        ybase = ybase(:, adaptivebin);
        ywn = ywn(:, adaptivebin);
        
    elseif adaptivebin==3
        % --- adaptive minus nonadaptive
        ybase = cellfun(@(x,y)y-x, ybase(:,1), ybase(:,2), 'UniformOutput', 0);
        ywn = cellfun(@(x,y)y-x, ywn(:,1), ywn(:,2), 'UniformOutput', 0);
    end
        % === replace with these
        OUTSTRUCT_XCOV.XcovgramBase = ybase;
        OUTSTRUCT_XCOV.XcovgramWN = ywn;
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
subplotrows=2;
subplotcols=3;
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

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 1, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');
ylim([-1 1]);

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

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 1, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');
ylim([-1 1]);


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

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 1, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');
ylim([-1 1]);

% -- plot each bird
pcolbirds = lt_make_plot_colors(max(allbnum), 0, 0);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
for i=1:max(allbnum)
   indsthis = allbnum==i;
   if ~any(indsthis)
       continue
   end
   
   ymat = squeeze(allDat_Diff(:,:,1,indsthis));

lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot, pcolbirds{i});

lt_plot_text(0, ((1/i))*clim(end), ['bnum' num2str(i)], pcolbirds{i}, 8);

end
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);


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

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 1, lagwindows, '', ...
    ffbinsedges_indstoplot);
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);
lt_plot_annotation(1, ['N=' num2str(sum(~isnan(squeeze(ymat(1,1,:)))))], 'm');
ylim([-1 1]);

% -- plot each bird
pcolbirds = lt_make_plot_colors(max(allbnum), 0, 0);
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
for i=1:max(allbnum)
   indsthis = allbnum==i;
   if ~any(indsthis)
       continue
   end
   
    ymat = squeeze(allDat_Diff(:,:,1,indsthis) - allDat_Diff(:,:,3,indsthis));

lt_neural_Coher_Plot(ymat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 2, '-', clim, 1, 0, lagwindows, '', ...
    ffbinsedges_indstoplot, pcolbirds{i});

lt_plot_text(0, ((1/i))*clim(end), ['bnum' num2str(i)], pcolbirds{i}, 8);

end
lt_plot_zeroline;
hsplots = [hsplots hsplot];
xlim(XLIM);


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



