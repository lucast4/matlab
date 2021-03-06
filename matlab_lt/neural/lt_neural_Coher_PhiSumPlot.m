function lt_neural_Coher_PhiSumPlot(OUTSTRUCT, SwitchStruct, ...
    useFFchangeInsteadOfLearnDir, onlygood)

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


if onlygood==1
    OUTSTRUCT = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct);
end

%% +++++++++++++++++++++++++++++++++++++++
% if cohdiff_norm_to_global==1
%     % ------- for each channel, get a global
%     Cohdiffmean = nan(length(OUTSTRUCT.bnum),1);
%     Rhomean = nan(length(OUTSTRUCT.bnum),1);
%     [inds_grp, inds_unique] = lt_tools_grp2idx(...
%         {OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
%     for j=inds_unique'
%         
%         indsthis = inds_grp==j;
%         assert(length(unique(OUTSTRUCT.motifname(indsthis)))== ...
%             length(OUTSTRUCT.motifname(indsthis)), 'should all be unique motifs');
%         
%         meancohdiff = mean(OUTSTRUCT.cohscal_WN_minus_base(indsthis));
%         meanrho = nanmean(OUTSTRUCT.cohscal_ff_rho_rhoSE_base(indsthis,1));
%         % ======= save outoput
%         Cohdiffmean(indsthis) = meancohdiff;
%         Rhomean(indsthis) = meanrho;
%     end
%     assert(~any(isnan(Cohdiffmean)));
% end

%%
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ======== ONE PLOT FOR EACH EXPERIMENT
[inds_grp, inds_unique, X_cell] = lt_tools_grp2idx(...
    {OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});

% RhoAll = [];
% CohDiffAll = [];
AllLearnDir = [];
% GrpIndAll = [];
AllPhiMean_BaseAngle = [];
AllGrpInd = [];
AllPhiPLV_BaseWN = [];

AllExpt_PhiMean_BaseAngle = [];
AllExpt_GrpInd = [];
AllExpt_AcrossChanVar = [];
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

    % ----- 1) Colelct data
    phimean_BaseWN = OUTSTRUCT.Phi_mean_BaseWN(indsthis,:); % [base, WN]
    phiPLV_BaseWN = OUTSTRUCT.Phi_PLV_BaseWN(indsthis,:);
    
    % -- plot
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots; hsplot];
    title({[bname '-' ename '-sw' num2str(swnum)], ['learndir ' num2str(learndir)]});
    
    circ_plot(phimean_BaseWN(:,1),'pretty','ko',true,'linewidth',2,'color','k');
    tmp1 = circ_stats(phimean_BaseWN(:,1));
    lt_plot_text(0,0, ['base:' num2str(tmp1.mean/pi) '*pi'], 'm');
    hold on;
    circ_plot(phimean_BaseWN(:,2),'pretty','ro',true,'linewidth',2,'color','r');
    
    % ---- what is change in phi?
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('diff in phi (WN minis base)');
    phiDiff = phimean_BaseWN(:,2) - phimean_BaseWN(:,1);
    circ_plot(phiDiff,'pretty','bo',true,'linewidth',2,'color','b');
    
    
    
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)');
%     
%     lt_plot(rho_base(:,1), coh_WNminBase, {'Xerrors', rho_base(:,2), 'Color', pcol});
%     xlim([-1 1]);
%     ylim([-1 1]);
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
    
    % =============== COLLECT
    AllPhiPLV_BaseWN = [AllPhiPLV_BaseWN; phiPLV_BaseWN];
    AllPhiMean_BaseAngle = [AllPhiMean_BaseAngle; phimean_BaseWN];
    AllGrpInd = [AllGrpInd; ones(size(phimean_BaseWN,1),1)*j];
    AllLearnDir = [AllLearnDir; ones(size(phimean_BaseWN,1),1)*learndir];
    
    tmp1 = circ_stats(phimean_BaseWN(:,1));
    tmp2 = circ_stats(phimean_BaseWN(:,2));
    AllExpt_PhiMean_BaseAngle = [AllExpt_PhiMean_BaseAngle; [tmp1.mean tmp2.mean]];
    AllExpt_GrpInd = [AllExpt_GrpInd; j];
    
    % ============= across channel stats
    if size(phimean_BaseWN,1)==1
    AllExpt_AcrossChanVar = [AllExpt_AcrossChanVar; [nan nan]];
    else
    AllExpt_AcrossChanVar = [AllExpt_AcrossChanVar; [tmp1.var tmp2.var]];
    end
    
end

% linkaxes(hsplots, 'xy');

%% =========== learn dir across experiments

AllExpt_LearnDir = grpstats(AllLearnDir, AllGrpInd);
assert(all(unique(AllGrpInd)==AllExpt_GrpInd));

AllExpt_PhiPLV_BaseWN = ...
    [grpstats(AllPhiPLV_BaseWN(:,1), AllGrpInd), grpstats(AllPhiPLV_BaseWN(:,2), AllGrpInd)];

%% ================ compute change in angle




%% ================ plot summary across expts

lt_figure; hold on;

lt_subplot(3,2,1); hold on;
title('base [neg=LMAN leads]');
circ_plot(AllExpt_PhiMean_BaseAngle(:,1),'pretty','ko',true,'linewidth',2,'color','k');
hold on;
[pval] = circ_rtest(AllExpt_PhiMean_BaseAngle(:,1));
lt_plot_pvalue(pval, 'rayleigh test', 1);

lt_subplot(3,2,2); hold on;
title('WN');
circ_plot(AllExpt_PhiMean_BaseAngle(:,2),'pretty','ro',true,'linewidth',2,'color','r');
hold on ;
[pval] = circ_rtest(AllExpt_PhiMean_BaseAngle(:,2));
lt_plot_pvalue(pval, 'rayleigh test', 1);

lt_subplot(3,2,3); hold on;
title('base(bk) - WN(rd) [dat=sw]');
circ_plot(AllExpt_PhiMean_BaseAngle(:,1),'pretty','ko',true,'linewidth',2,'color','k');
hold on;
circ_plot(AllExpt_PhiMean_BaseAngle(:,2),'pretty','ro',true,'linewidth',2,'color','r');
hold on;
[pval] = circ_wwtest(AllExpt_PhiMean_BaseAngle(:,1), AllExpt_PhiMean_BaseAngle(:,2));
lt_plot_pvalue(pval, 'ww test [diff means]');

lt_subplot(3,2,4); hold on;
title('PLV (dat = mean across chans');
xlabel('base, WN');
ylabel('mean PLV(over trials)');
x = [1 2];
y = AllExpt_PhiPLV_BaseWN;
plot(x, y', '-ok');
xlim([0 3]);
lt_plot(x+0.2, nanmean(y,1), {'Errors', lt_sem(y)});
[~, p] = ttest(y(:,1), y(:,2));
lt_plot_pvalue(p, 'ttest', 1);


lt_subplot(3,2,5); hold on;
title('learn dir -1');
y = AllExpt_PhiMean_BaseAngle(AllExpt_LearnDir==-1, :);
circ_plot(y(:,1),'pretty','ko',true,'linewidth',2,'color','k');
hold on;
circ_plot(y(:,2),'pretty','ro',true,'linewidth',2,'color','r');
hold on;
[pval] = circ_wwtest(y(:,1), y(:,2));
lt_plot_pvalue(pval, 'ww test [diff means]');

lt_subplot(3,2,6); hold on;
title('learn dir +1');
y = AllExpt_PhiMean_BaseAngle(AllExpt_LearnDir==1, :);
circ_plot(y(:,1),'pretty','ko',true,'linewidth',2,'color','k');
hold on;
circ_plot(y(:,2),'pretty','ro',true,'linewidth',2,'color','r');
hold on;
[pval] = circ_wwtest(y(:,1), y(:,2));
lt_plot_pvalue(pval, 'ww test [diff means]');


%% ================ DO CHANNELS BECOME MORE similar?
% i.e. for each switch, one val for variance across chans (variance =
% 1-PLV);

lt_figure; hold on;

% ===
lt_subplot(3,2,1); hold on;
xlabel('base - WN');
ylabel('phi var (i.e. 1-PLV) across chans');

x = [1 2];
y = AllExpt_AcrossChanVar;
plot(x, y', '-ok');
xlim([0 3]);
lt_plot(x+0.2, nanmean(y,1), {'Errors', lt_sem(y)});

[~, p] = ttest(y(:,1), y(:,2));
lt_plot_pvalue(p, 'ttest', 1);

%% ================ one value of mean for each expt


circ_stats(phimean_BaseWN(:,1))

lt_neural_QUICK_PhaseLockVal(phimean_BaseWN(:,1))
%% wrap angles to 2pi
AllPhiMean_BaseAngle = wrapTo2Pi(AllPhiMean_BaseAngle);


%% ======================


%% ====================== PLOT SUMMARY [SHOW EACH EXPT, AND COMBNIE]

lt_figure; hold on; 
xlabel('bird - expt - swnum');
ylabel('phi (plotted 2x)');
title('sq value = learndir');
x = [AllGrpInd; AllGrpInd];
y = [AllPhiMean_BaseAngle; AllPhiMean_BaseAngle+2*pi];

ldir = grpstats(AllLearnDir, AllGrpInd);

plot(x-0.2, y(:,1), 'ok');
plot(x+0.2, y(:,2), 'or');
plot([x-0.2 x+0.2]', y', '-k');
plot(unique(AllGrpInd), ldir, 'sm');
lt_plot_zeroline;
set(gca, 'XTick', 1:max(AllGrpInd), 'XTickLabel', unique(X_cell));
rotateXLabels(gca, 90);


%% ====================== PLOT EACH EXPT SEPARATELY ON CIRCLE


%% ====================== IS THERE A CONSISTENT DIRECTION CHANGE?
lt_figure; hold on;

% ==== all channel pairs
lt_subplot(3,2,1); hold on;
xlabel('BASE -- WN');
ylabel('phi');

x = [1 2];
y = [AllPhiMean_BaseAngle; AllPhiMean_BaseAngle+2*pi];

plot(x, y', '-k');




%% ====================== COMBINE ALL DATAPOINTS INTO ONE PLOT
% DATAPOINT = channelpair;
% 
% % -------------- 1) POSITIVE LEARNING
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title('All dat (pos learn)');
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)');
% 
%     x = RhoAll(LearnDirAll==1);
%     y = CohDiffAll(LearnDirAll==1);
%     plot(x,y, 'ok');
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
%     
% % -------------- 1) NEGATIVE LEARNING
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title('All dat (neg learn)');
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)');
% 
%     x = RhoAll(LearnDirAll== -1);
%     y = CohDiffAll(LearnDirAll== -1);
%     plot(x,y, 'ok');
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
% 
% % ---------------- 2) COMBINE ONTO SAME PLOT
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title('All dat');
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)[flip if neg learn');
% 
%     x = RhoAll;
%     y = CohDiffAll;
%     y(LearnDirAll==-1) = -y(LearnDirAll==-1);
%     plot(x,y, 'ok');
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
%     lt_regress(y, x, 0, 0, 0, 1, 'r',1);
% 
%% ====== one plot per SWITCH
% RhoAll = grpstats(RhoAll, GrpIndAll, {'mean'});
% CohDiffAll = grpstats(CohDiffAll, GrpIndAll, {'mean'});
% LearnDirAll = grpstats(LearnDirAll, GrpIndAll, {'mean'});
%     
% % -------------- 1) POSITIVE LEARNING
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title('DAT=switch (pos learn)');
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)');
% 
%     x = RhoAll(LearnDirAll==1);
%     y = CohDiffAll(LearnDirAll==1);
%     plot(x,y, 'ok');
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
%     
% % -------------- 1) NEGATIVE LEARNING
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title('DAT=switch (neg learn)');
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)');
% 
%     x = RhoAll(LearnDirAll== -1);
%     y = CohDiffAll(LearnDirAll== -1);
%     plot(x,y, 'ok');
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
% 
% % ---------------- 2) COMBINE ONTO SAME PLOT
%     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%     hsplots = [hsplots; hsplot];
%     title('DAT=switch');
%     xlabel('rho (coh vs. ff, BASE)');
%     ylabel('change in coh(WN - base)[flip if neg learn');
% 
%     x = RhoAll;
%     y = CohDiffAll;
%     y(LearnDirAll==-1) = -y(LearnDirAll==-1);
%     plot(x,y, 'ok');
%     lt_plot_zeroline;
%     lt_plot_zeroline_vert;
%     lt_regress(y, x, 0, 0, 0, 1, 'r',1);
