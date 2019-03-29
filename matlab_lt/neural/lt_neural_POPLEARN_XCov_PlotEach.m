function lt_neural_POPLEARN_XCov_PlotAll(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    plotNotminShuff, onlygoodexpt, clim, btoplot, etoplot, stoplot)

%%
if onlygoodexpt==1
    OUTSTRUCT_XCOV = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, ...
        SwitchStruct, 'xcov_spikes');
end


%% 1/9/19 - lt plots each expt, targ, same, diff, xcov (minus shift predictor

% ==========
[indsgrp, indsgrpUni] = lt_tools_grp2idx({OUTSTRUCT_XCOV.bnum, ...
    OUTSTRUCT_XCOV.enum, OUTSTRUCT_XCOV.swnum, OUTSTRUCT_XCOV.motifnum});

figcount=1;
subplotrows=7;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for i=1:length(indsgrpUni)
    
    indstmp = find(indsgrp==indsgrpUni(i));

    % ========== expt details
    bnum = unique(OUTSTRUCT_XCOV.bnum(indstmp));
    enum = unique(OUTSTRUCT_XCOV.enum(indstmp));
    swnum = unique(OUTSTRUCT_XCOV.swnum(indstmp));
    mname = unique(OUTSTRUCT_XCOV.motifname(indstmp)); mname = mname{1};
    
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    if bnum~=btoplot | enum ~=etoplot |swnum ~=stoplot
        continue
    end
    
    % ================    
    covbase = OUTSTRUCT_XCOV.XcovBase(indstmp,:);
    covWN = OUTSTRUCT_XCOV.XcovWN(indstmp,:);
    
    
    % 1) BASE AND WN SEPARATE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=base; rd=WN)');
    plot(PARAMS.Xcov_ccLags, covbase', '-k');
    plot(PARAMS.Xcov_ccLags, covWN', '-r');
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % 1) BASE AND WN SEPARATE [means]
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(mname);
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
%     title(ptit);
    ylabel('xcov (k=base; rd=WN)');
    y = covWN-covbase;
    ymean = mean(y,1);
    ysem = lt_sem(y);
    plot(PARAMS.Xcov_ccLags, y', '-b');
    if length(ysem)>1
        shadedErrorBar(PARAMS.Xcov_ccLags, ymean, ysem, {'Color', 'k'}, 1);
    end
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
end

linkaxes(hsplots, 'xy');
axis tight;



