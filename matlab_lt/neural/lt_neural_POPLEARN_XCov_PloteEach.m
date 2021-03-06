function lt_neural_POPLEARN_XCov_PlotAll(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    plotNotminShuff, onlygoodexpt, clim)

%%
if onlygoodexpt==1
    OUTSTRUCT_XCOV = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
end


%% 1/9/19 - lt plots each expt, targ, same, diff, xcov (minus shift predictor

% ==========
[indsgrp, indsgrpUni] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.swnum});

figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(indsgrpUni)
    hsplots = [];
    
    % ========== expt details
    bnum = unique(OUTSTRUCT_XCOV.bnum(indsgrp==indsgrpUni(i)));
    enum = unique(OUTSTRUCT_XCOV.enum(indsgrp==indsgrpUni(i)));
    swnum = unique(OUTSTRUCT_XCOV.swnum(indsgrp==indsgrpUni(i)));
    
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    
    % ################################ TARGET
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==1;
    ptit = 'TARG';
    
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
    end
    
    % 1) BASE AND WN SEPARATE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    plot(PARAMS.Xcov_ccLags, covbase', '-k');
    plot(PARAMS.Xcov_ccLags, covWN', '-r');
    
    % 1) BASE AND WN SEPARATE [means]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    y = mean(covbase,1);
    ysem = lt_sem(covbase);
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, y, ysem, {'Color', 'k'}, 1);
        y = mean(covWN,1);
        ysem = lt_sem(covWN);
        shadedErrorBar(PARAMS.Xcov_ccLags, y, ysem, {'Color', 'r'}, 1);
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert
    
    % 1) WN MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(ptit);
    ylabel('xcov (k=base; rd=WN)');
    y = covWN-covbase;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    plot(PARAMS.Xcov_ccLags, y', '-b');
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, ymean, ysem, {'Color', 'k'}, 1);
    end
    lt_plot_zeroline;
    
    
    % ################################ SAME
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==1;
    ptit = 'SAME';
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
    end
    
    % 1) BASE AND WN SEPARATE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    if any(indstmp)
        plot(PARAMS.Xcov_ccLags, covbase', '-k');
        plot(PARAMS.Xcov_ccLags, covWN', '-r');
    end
    % 1) BASE AND WN SEPARATE [means]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    y = mean(covbase,1);
    ysem = lt_sem(covbase);
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, y, ysem, {'Color', 'k'}, 1);
        y = mean(covWN,1);
        ysem = lt_sem(covWN);
        shadedErrorBar(PARAMS.Xcov_ccLags, y, ysem, {'Color', 'r'}, 1);
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert
    
    % 1) WN MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(ptit);
    ylabel('xcov (k=base; rd=WN)');
    y = covWN-covbase;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    if any(indstmp)
        plot(PARAMS.Xcov_ccLags, y', '-b');
    end
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, ymean, ysem, {'Color', 'k'}, 1);
    end
    lt_plot_zeroline;
    
    
    % ################################ DIFF
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==0;
    ptit = 'DIFF';
    
    if plotNotminShuff==1
        covbase = OUTSTRUCT_XCOV.XcovBase_NoMinShuff(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN_NoMinShuff(indstmp,:);
    elseif plotNotminShuff==0
        covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
        covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
    end
    
    % 1) BASE AND WN SEPARATE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    plot(PARAMS.Xcov_ccLags, covbase', '-k');
    plot(PARAMS.Xcov_ccLags, covWN', '-r');
    
    % 1) BASE AND WN SEPARATE [means]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    y = mean(covbase,1);
    ysem = lt_sem(covbase);
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, y, ysem, {'Color', 'k'}, 1);
        y = mean(covWN,1);
        ysem = lt_sem(covWN);
        shadedErrorBar(PARAMS.Xcov_ccLags, y, ysem, {'Color', 'r'}, 1);
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert
    
    
    % 1) WN MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(ptit);
    ylabel('xcov (k=base; rd=WN)');
    y = covWN-covbase;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    plot(PARAMS.Xcov_ccLags, y', '-b');
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, ymean, ysem, {'Color', 'k'}, 1);
    end
    lt_plot_zeroline;
    
    linkaxes(hsplots);
    axis tight;
    
    
    %% =========== XGRAM PLOTS
    % ################################ TARGET
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==1;
    ptit = 'TARG';
    
    % -- BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[BASE]'])
    xgram = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % -- WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[WN]'])
    xgram = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % -- WN minus Base
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[WN-base]'])
    xgram1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
    xgram2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
    xgram = xgram1-xgram2;
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    
    
    
    % ################################ SAME
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==1;
    ptit = 'SAME';
    
    % -- BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[BASE]'])
    xgram = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % -- WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[WN]'])
    xgram = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % -- WN minus Base
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[WN-base]'])
    xgram1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
    xgram2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
    xgram = xgram1-xgram2;
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % ################################ DIFF
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==0;
    ptit = 'DIFF';
    
    % -- BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[BASE]'])
    xgram = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % -- WN
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[WN]'])
    xgram = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
    % -- WN minus Base
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([ptit '[WN-base]'])
    xgram1 = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramWN(indstmp));
    xgram2 = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.XcovgramBase(indstmp));
    xgram = xgram1-xgram2;
    lt_neural_Coher_Plot(xgram, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], clim);
    axis tight;
    lt_plot_zeroline;
    
end


