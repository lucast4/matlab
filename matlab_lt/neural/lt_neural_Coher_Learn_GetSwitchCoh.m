function Yallswitch = lt_neural_Coher_Learn_GetSwitchCoh(SwitchStruct, OUTSTRUCT,...
    useabs, plotON, PARAMS, vssametype)
%% lt 10/12/18 - extracts mean coherence by switch (one matrix for pre andone for post)
% vssametype = 1; then instead of nontarg comapres to same type

% ============== PARAMS
% useabs = 1; % if 1, then absolute values (wn minus diff)
% plotON =0;

tbins = PARAMS.tbins;
ffbins = PARAMS.ffbins;

%% ===============
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});
indsgrp_unique = unique(indsgrp)';


% ===============
figcount=1;
subplotrows=3;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

clim = [-0.4 0.4];

% ==== to plot mean across all switches
Yallswitch = {}; % switch x (pre, post)
for j=indsgrp_unique
    
    bname = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).birdname;
    ename = SwitchStruct.bird(unique(OUTSTRUCT.bnum(indsgrp==j))).exptnum(unique(OUTSTRUCT.enum(indsgrp==j))).exptname;
    swnum = unique(OUTSTRUCT.switch(indsgrp==j));
    
    % ======================= TRAIN MINUS BASELINE COHERENCE
    Yall = {};
    % --- targets
    indsthis = OUTSTRUCT.istarg==1 & indsgrp==j;
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    if useabs==1
        cohmat = abs(cohmat);
    end
    %     Yall = [Yall cohmat];
    cohmean = nanmean(cohmat,3);
    Yall = [Yall cohmean];
    
    % --- nontarg
    if vssametype==1
        indsthis = OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1 & indsgrp==j;
    else
        indsthis = OUTSTRUCT.istarg==0 & indsgrp==j;
    end
    cohmat = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohMean_WNminusBase(indsthis));
    if useabs==1
        cohmat = abs(cohmat);
    end
    %     Yall = [Yall cohmat];
    cohmean = nanmean(cohmat,3);
    Yall = [Yall cohmean];
    if isempty(cohmean)
        continue
    end
    % --------------------- targ minus nontarg
    cohdiff = nanmean(Yall{1},3) - nanmean(Yall{2},3);
    
    % ================= PLOT
    if plotON==1
        % ----- 1) TARG (WN minsu base)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('targ(WN-base)[all chan pairs]');
        title([bname '-' ename '-sw' num2str(swnum)]);
        lt_neural_Coher_Plot(Yall{1}, tbins, ffbins, 1, '-', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('targ(WN-base)');
        lt_neural_Coher_Plot(Yall{1}, tbins, ffbins, 2, '-', clim);
        lt_plot_zeroline;
        
        % ----- 2) NONTARG (WN minsu base)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('nontarg(WN-base)');
        lt_neural_Coher_Plot(Yall{2}, tbins, ffbins, 1, '-', clim);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        ylabel('nontarg(WN-base)');
        lt_neural_Coher_Plot(Yall{2}, tbins, ffbins, 2, '-', clim);
        lt_plot_zeroline;
        
        % ----- 3) TARG-NONTARG (WN minsu base)
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohdiff, tbins, ffbins, 1, '-', clim);
        ylabel('targ(WN-base) - nontarg(WN-base)');
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        lt_neural_Coher_Plot(cohdiff, tbins, ffbins, 2, '-', clim);
        lt_plot_zeroline;
    end
    
    % ============================= COLLECT
    Yallswitch = [Yallswitch; Yall];
end