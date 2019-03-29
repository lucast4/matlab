function lt_neural_POPLEARN_XCov_PlotAll(OUTSTRUCT_XCOV, SwitchStruct, PARAMS, ...
    plotNotminShuff, onlygoodexpt, clim, btoplot, etoplot, stoplot, ...
    epochstoplot)

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
subplotcols=4;
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
    covbase = OUTSTRUCT_XCOV.Xcovslice_ffsplits_base(indstmp);    
    covWN = OUTSTRUCT_XCOV.Xcovslice_ffsplits_epochs(indstmp);
    
    % --- % extract fflo and ffhi versions
    covbase_lo = cellfun(@(x)x{1}, covbase, 'UniformOutput', 0);
    covbase_hi = cellfun(@(x)x{2}, covbase, 'UniformOutput', 0);
    
    covWN_lo = cellfun(@(x)x{1}, covWN, 'UniformOutput', 0);
    covWN_hi = cellfun(@(x)x{2}, covWN, 'UniformOutput', 0);
    
    % --- pick out epochs of interest
    covWN_lo = cellfun(@(x)mean(x(:,:, epochstoplot),3), covWN_lo, 'UniformOutput', 0);
    covWN_hi = cellfun(@(x)mean(x(:,:, epochstoplot),3), covWN_hi, 'UniformOutput', 0);
    
    % ========================
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    ylabel('xcov (k=lopitch; rd=Hipitch)');
    xlabel('BASE');
    
    ylo = cell2mat(covbase_lo);
    yhi = cell2mat(covbase_hi);
    plot(PARAMS.Xcov_ccLags, ylo', '-k');
    plot(PARAMS.Xcov_ccLags, yhi', '-r');
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ========================
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(mname);
    hsplots = [hsplots hsplot];
    xlabel('BASE');

    shadedErrorBar(PARAMS.Xcov_ccLags, mean(ylo,1), lt_sem(ylo), {'Color', 'k'}, 1);
    shadedErrorBar(PARAMS.Xcov_ccLags, mean(yhi,1), lt_sem(yhi), {'Color', 'r'}, 1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ========================
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(swnum)]);
    hsplots = [hsplots hsplot];
    xlabel('WN');
    
    ylo = cell2mat(covWN_lo);
    yhi = cell2mat(covWN_hi);
    plot(PARAMS.Xcov_ccLags, ylo', '-k');
    plot(PARAMS.Xcov_ccLags, yhi', '-r');
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ========================
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(mname);
    hsplots = [hsplots hsplot];
    xlabel('WN');

    shadedErrorBar(PARAMS.Xcov_ccLags, mean(ylo,1), lt_sem(ylo), {'Color', 'k'}, 1);
    shadedErrorBar(PARAMS.Xcov_ccLags, mean(yhi,1), lt_sem(yhi), {'Color', 'r'}, 1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;

end

linkaxes(hsplots, 'xy');
axis tight;



