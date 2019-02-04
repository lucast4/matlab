function lt_neural_POPLEARN_XCov_PlotAllGram(OUTSTRUCT_XCOV, SwitchStruct, ...
    PARAMS, clim)


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

    covbase = OUTSTRUCT_XCOV.XcovgramBase(indstmp);
    covWN = OUTSTRUCT_XCOV.XcovgramWN(indstmp);
    % --- mean across all pairs
    covbase = lt_neural_Coher_Cell2Mat(covbase);
    covWN = lt_neural_Coher_Cell2Mat(covWN);
    
    % 1) BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[BASE]']);
    xlabel('xbin (center');
    ylabel('lag (neg = LMAN lead');
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covbase, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);
    
    % 1) WN 
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[WN]']);
    xlabel('xbin (center');
    ylabel(ptit);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covWN, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);

    % wn MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[WN-base]']);
    xlabel('xbin (center');
    ylabel(ptit);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covWN-covbase, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);
    
    
    
    
    % ################################ SAME
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==1;
    ptit = 'SAME';

    covbase = OUTSTRUCT_XCOV.XcovgramBase(indstmp);
    covWN = OUTSTRUCT_XCOV.XcovgramWN(indstmp);
    % --- mean across all pairs
    covbase = lt_neural_Coher_Cell2Mat(covbase);
    covWN = lt_neural_Coher_Cell2Mat(covWN);
    
    % 1) BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[BASE]']);
    xlabel('xbin (center');
    ylabel('lag (neg = LMAN lead');
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covbase, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);
    
    % 1) WN 
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[WN]']);
    xlabel('xbin (center');
    ylabel(ptit);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covWN, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);

    % wn MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[WN-base]']);
    xlabel('xbin (center');
    ylabel(ptit);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covWN-covbase, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);
    
    
    
    % ################################ DIFF
    indstmp = indsgrp==indsgrpUni(i) & OUTSTRUCT_XCOV.istarg==0 & OUTSTRUCT_XCOV.issame==0;
    ptit = 'DIFF';

    covbase = OUTSTRUCT_XCOV.XcovgramBase(indstmp);
    covWN = OUTSTRUCT_XCOV.XcovgramWN(indstmp);
    % --- mean across all pairs
    covbase = lt_neural_Coher_Cell2Mat(covbase);
    covWN = lt_neural_Coher_Cell2Mat(covWN);
    
    % 1) BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[BASE]']);
    xlabel('xbin (center');
    ylabel('lag (neg = LMAN lead');
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covbase, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);
    
    % 1) WN 
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[WN]']);
    xlabel('xbin (center');
    ylabel(ptit);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covWN, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);

    % wn MINUS BASE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum) '[WN-base]']);
    xlabel('xbin (center');
    ylabel(ptit);
    hsplots = [hsplots hsplot];
    lt_neural_Coher_Plot(covWN-covbase, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, [], '', clim);
    
end