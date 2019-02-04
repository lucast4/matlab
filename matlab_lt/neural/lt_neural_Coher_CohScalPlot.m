function lt_neural_Coher_CohScalPlot(OUTSTRUCT, SwitchStruct, cohdiff_norm_to_global, ...
    rho_norm_to_global, centerdata, onlyfirstswitch, useFFchangeInsteadOfLearnDir, ...
    onlygoodexpt)

%% lt 11/11/18 - plots summary across and within experiments


% cohdiff_norm_to_global = 1; % if 1, then finds (for a given channel)
% rho_norm_to_global = 0;
% centerdata = 0; % then for each switch, centers data (i.e. across channels)
% % NOTE: this is useful if want to ask about, within a tgiven expt, is there
% % correlation.
%
% % change across syllables.
% onlyfirstswitch = 1; % if 0, then all siwtches.
%
% useFFchangeInsteadOfLearnDir=0; % then sign(ff(WN)-ff(base)); if 0, then learndir
%


%% ======= 2) FILTER DATA
% == filter by specific types of switches.

if onlygoodexpt==1
    OUTSTRUCT = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct);
end

%% +++++++++++++++++++++++++++++++++++++++
if cohdiff_norm_to_global==1
    % ------- for each channel, get a global across motifs.
    Cohdiffmean = nan(length(OUTSTRUCT.bnum),1);
    Rhomean = nan(length(OUTSTRUCT.bnum),1);
    [inds_grp, inds_unique] = lt_tools_grp2idx(...
        {OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
    for j=inds_unique'
        
        indsthis = inds_grp==j;
        assert(length(unique(OUTSTRUCT.motifname(indsthis)))== ...
            length(OUTSTRUCT.motifname(indsthis)), 'should all be unique motifs');
        
        meancohdiff = mean(OUTSTRUCT.cohscal_WN_minus_base(indsthis));
        meanrho = nanmean(OUTSTRUCT.cohscal_ff_rho_rhoSE_base(indsthis,1));
        % ======= save outoput
        Cohdiffmean(indsthis) = meancohdiff;
        Rhomean(indsthis) = meanrho;
    end
    assert(~any(isnan(Cohdiffmean)));
end

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ======== ONE PLOT FOR EACH EXPERIMENT
[inds_grp, inds_unique] = lt_tools_grp2idx(...
    {OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});

RhoAll = [];
CohDiffAll = [];
LearnDirAll = [];
GrpIndAll = [];
for j=inds_unique'
    
    indsthis = inds_grp==j & OUTSTRUCT.istarg==1;
    
    % ======== for each channel pair, collect correlation and change in coh
    % -------2 ) Collect metadata
    learndir = unique(OUTSTRUCT.learndirTarg(indsthis));
    if useFFchangeInsteadOfLearnDir==1
        learndir = sign(unique(OUTSTRUCT.ff_WNminusBase(indsthis)));
    end
    
    if learndir==1
        pcol = 'b';
    elseif learndir==-1
        pcol = 'r';
    else
        pcol = [0.7 0.7 0.7];
    end
    
    bnum = unique(OUTSTRUCT.bnum(indsthis));
    enum = unique(OUTSTRUCT.enum(indsthis));
    swnum = unique(OUTSTRUCT.switch(indsthis));
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    if onlyfirstswitch==1 & swnum~=1
        continue
    end
    
    % ----- 1) Colelct data
    rho_base = OUTSTRUCT.cohscal_ff_rho_rhoSE_base(indsthis,:);
    assert(~any(isnan(rho_base(:))), 'shoud not be nan since have ff since is targ');
    
    coh_WNminBase = OUTSTRUCT.cohscal_WN_minus_base(indsthis);
    if cohdiff_norm_to_global==1
        % --- get global mean of coh diff
        cohglobmean = Cohdiffmean(indsthis);
        coh_WNminBase = coh_WNminBase-cohglobmean;
    end
    if rho_norm_to_global==1
        % --- norm to global (within channel, across motifs);
        rhoglobmean = Rhomean(indsthis);
        rho_base(:,1) = rho_base(:,1)-rhoglobmean;
    end
    
    if centerdata==1
        rho_base(:,1) = rho_base(:,1) - mean(rho_base(:,1));
        coh_WNminBase = coh_WNminBase - mean(coh_WNminBase);
    end
    
    % -- plot
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title({[bname '-' ename '-sw' num2str(swnum)], ['learndir ' num2str(learndir)]});
    xlabel('rho (coh vs. ff, BASE)');
    ylabel('change in coh(WN - base)');
    
    lt_plot(rho_base(:,1), coh_WNminBase, {'Xerrors', rho_base(:,2), 'Color', pcol});
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % =============== COLLECT
    RhoAll = [RhoAll; rho_base(:,1)];
    CohDiffAll = [CohDiffAll; coh_WNminBase];
    LearnDirAll = [LearnDirAll; ones(size(coh_WNminBase))*learndir];
    GrpIndAll = [GrpIndAll; j*ones(size(coh_WNminBase))];
end

linkaxes(hsplots, 'xy');

%% ====================== limits

YLIM = [-max(abs(CohDiffAll))-0.01 max(abs(CohDiffAll))+0.01];
XLIM = [-max(abs(RhoAll))-0.01 max(abs(RhoAll))+0.01];
    
%% ====================== COMBINE ALL DATAPOINTS INTO ONE PLOT
% DATAPOINT = channelpair;

% -------------- 1) POSITIVE LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('All dat (pos learn)');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)');

x = RhoAll(LearnDirAll==1);
y = CohDiffAll(LearnDirAll==1);
plot(x,y, 'ok');
xlim(XLIM); ylim(YLIM);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% -------------- 1) NEGATIVE LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('All dat (neg learn)');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)');

x = RhoAll(LearnDirAll== -1);
y = CohDiffAll(LearnDirAll== -1);
plot(x,y, 'ok');
xlim(XLIM); ylim(YLIM);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---------------- 2) COMBINE ONTO SAME PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('All dat');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)[flip if neg learn');

x = RhoAll;
y = CohDiffAll;
y(LearnDirAll==-1) = -y(LearnDirAll==-1);
plot(x,y, 'ok');
lt_plot_zeroline;
lt_plot_zeroline_vert;
lt_regress(y, x, 0, 0, 0, 1, 'r',1);
xlim(XLIM); ylim(YLIM);


% ---------------- 2) COMBINE ONTO SAME PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('All dat');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)[no flip]');

x = RhoAll;
y = CohDiffAll;
plot(x,y, 'ok');
lt_plot_zeroline;
lt_plot_zeroline_vert;
lt_regress(y, x, 0, 0, 0, 1, 'r',1);
xlim(XLIM); ylim(YLIM);

%% ====== one plot per SWITCH
RhoAll = grpstats(RhoAll, GrpIndAll, {'mean'});
CohDiffAll = grpstats(CohDiffAll, GrpIndAll, {'mean'});
LearnDirAll = grpstats(LearnDirAll, GrpIndAll, {'mean'});

% -------------- 1) POSITIVE LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('DAT=switch (pos learn)');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)');

x = RhoAll(LearnDirAll==1);
y = CohDiffAll(LearnDirAll==1);
plot(x,y, 'ok');
xlim(XLIM); ylim(YLIM);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% -------------- 1) NEGATIVE LEARNING
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('DAT=switch (neg learn)');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)');

x = RhoAll(LearnDirAll== -1);
y = CohDiffAll(LearnDirAll== -1);
plot(x,y, 'ok');
xlim(XLIM); ylim(YLIM);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---------------- 2) COMBINE ONTO SAME PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('DAT=switch');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)[flip if neg learn');

x = RhoAll;
y = CohDiffAll;
y(LearnDirAll==-1) = -y(LearnDirAll==-1);
plot(x,y, 'ok');
lt_plot_zeroline;
lt_plot_zeroline_vert;
lt_regress(y, x, 0, 0, 0, 1, 'r',1);
xlim(XLIM); ylim(YLIM);

% ---------------- 2) COMBINE ONTO SAME PLOT
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots; hsplot];
title('DAT=switch');
xlabel('rho (coh vs. ff, BASE)');
ylabel('change in coh(WN - base)[no flip]');

x = RhoAll;
y = CohDiffAll;
plot(x,y, 'ok');
lt_plot_zeroline;
lt_plot_zeroline_vert;
lt_regress(y, x, 0, 0, 0, 1, 'r',1);
xlim(XLIM); ylim(YLIM);
