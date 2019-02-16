function lt_neural_POPLEARN_Cons_PlotOver(OUTSTRUCT_units, SwitchStruct, ...
    SummaryStruct, birdtoplot, onlygoodexpt)
%% lt 2/5/19 - Overview plots for consistency of smoothed FR (x-trials)

% birdtoplot = []; % empty for all. 

%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter outstruct
    expttype = 'xcov_spikes';
    [OUTSTRUCT_units] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_units, SwitchStruct, expttype);
    
end


%% =========== [PLOTS] One for each switch

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT_units.bnum, OUTSTRUCT_units.enum, OUTSTRUCT_units.switch, ...
    OUTSTRUCT_units.neurID});

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(indsgrpU)
    
    indsthis = indsgrp==indsgrpU(i);
    
    rho_baseWn = OUTSTRUCT_units.xtrialFrRho_BaseWn(indsthis, :);
    istarg = OUTSTRUCT_units.istarg(indsthis);
    issame = OUTSTRUCT_units.issame(indsthis);
    motifname = OUTSTRUCT_units.motifname(indsthis);
    
    % -- details
    bnum = unique(OUTSTRUCT_units.bnum(indsthis));
    enum = unique(OUTSTRUCT_units.enum(indsthis));
    sw = unique(OUTSTRUCT_units.switch(indsthis));
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    nID = unique(OUTSTRUCT_units.neurID(indsthis));
    bregion = SummaryStruct.birds(bnum).neurons(nID).NOTE_Location;
    
    if ~isempty(birdtoplot)
       if ~ismember(birdtoplot, bnum)
           continue
       end
    end
    % ====== PLOT
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([bname '-' ename '-sw' num2str(sw) '-nID' num2str(nID) '[' bregion ']']);
    xlabel('BASE (FR corr, x-trials)');
    ylabel('WN');
    
    rhomean = cellfun(@nanmean, rho_baseWn);
    rhosem = cellfun(@lt_sem, rho_baseWn);
    
    % ------ target
    indstmp = istarg==1;
    pcol = 'r';
    
    x = rhomean(indstmp, 1);
    xerr = rhosem(indstmp,1);
    y = rhomean(indstmp, 2);
    yerr = rhosem(indstmp, 2);
    lt_plot(x, y, {'Errors', yerr, 'Xerrors', xerr, 'Color', pcol});
    
    
    % ------ same
    indstmp = istarg==0 & issame==1;
    pcol = 'b';
    
    x = rhomean(indstmp, 1);
    xerr = rhosem(indstmp,1);
    y = rhomean(indstmp, 2);
    yerr = rhosem(indstmp, 2);
    lt_plot(x, y, {'Errors', yerr, 'Xerrors', xerr, 'Color', pcol});
    
    
        % ------ diff
    indstmp = istarg==0 & issame==0;
    pcol = 'k';
    
    x = rhomean(indstmp, 1);
    xerr = rhosem(indstmp,1);
    y = rhomean(indstmp, 2);
    yerr = rhosem(indstmp, 2);
    lt_plot(x, y, {'Errors', yerr, 'Xerrors', xerr, 'Color', pcol});
    
    % ----------------------------------- FORMAT PLOT
    lt_plot_makesquare_plot45line(gca, 'k', []);
  
end

linkaxes(hsplots, 'xy');
