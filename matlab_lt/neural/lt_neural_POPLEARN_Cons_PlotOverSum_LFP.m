function lt_neural_POPLEARN_Cons_PlotOverSum_LFP(OUTSTRUCT_lfp, OUTSTRUCT, ...
    SwitchStruct, SummaryStruct, birdtoplot, onlygoodexpt, plotlevel, ...
    corrwind, expttoplot)
%% lt 2/5/19 - Overview plots for consistency of smoothed FR (x-trials)

% birdtoplot = []; % empty for all.

%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter outstruct
    expttype = 'xcov_spikes';
    [OUTSTRUCT_lfp] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_lfp, SwitchStruct, expttype);

% [OUTSTRUCT_lfp] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_lfp, SwitchStruct);
end

%% ========= filter by bird/expt
if ~isempty(birdtoplot)
    indstokeep = ismember(OUTSTRUCT_lfp.bnum, birdtoplot);
    OUTSTRUCT_lfp = lt_structure_subsample_all_fields(OUTSTRUCT_lfp, indstokeep, 1);
end

if ~isempty(expttoplot)
    indstokeep = ismember(OUTSTRUCT_lfp.enum, expttoplot);
    OUTSTRUCT_lfp = lt_structure_subsample_all_fields(OUTSTRUCT_lfp, indstokeep, 1);
end

%% ========= GET MEAN RHO FOR EACH CASE

OUTSTRUCT_lfp.xtrialFrRho_Means = cellfun(@nanmean, OUTSTRUCT_lfp.xtrialFrRho_BaseWn);
OUTSTRUCT_lfp.chanpair = OUTSTRUCT_lfp.neurID; % just ad hoc for the code that expects this (GrpStats)

%% ========= SPLIT DATASET BASED ON BRAIN REGION

inds_RA = strcmp(OUTSTRUCT_lfp.bregion, 'RA');
OUTSTRUCT_units_RA = lt_structure_subsample_all_fields(OUTSTRUCT_lfp, inds_RA, 1);

inds_LMAN = strcmp(OUTSTRUCT_lfp.bregion, 'LMAN');
OUTSTRUCT_units_LMAN = lt_structure_subsample_all_fields(OUTSTRUCT_lfp, inds_LMAN, 1);

%% COLLECT DATA

% ========= LAERNING
if strcmp(plotlevel, 'chanpair')
    [allbnum_RA, allenum_RA, allswnum_RA, allDat_RA] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_units_RA, 'xtrialFrRho_Means');
    [allbnum_LMAN, allenum_LMAN, allswnum_LMAN, allDat_LMAN] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_units_LMAN, 'xtrialFrRho_Means');
elseif strcmp(plotlevel, 'switch')
    [~,~,~,~, allbnum_RA, allenum_RA, allswnum_RA, allDat_RA] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_units_RA, 'xtrialFrRho_Means');
    [~,~,~,~, allbnum_LMAN, allenum_LMAN, allswnum_LMAN, allDat_LMAN] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_units_LMAN, 'xtrialFrRho_Means');    
end



%% ========== COLLECT COHERENCE SCALAR BASE AND WN
disp('NEED TO MODIFY CODE TO MATCH EXACTLY THE TRIALS USED TO GET COHERENCE SCALAR AND LFP CORRELATIONS');
disp('DO that by referring back to OUTSTRUCT_lfp here')
pause;

% CohScal_BaseWN = nan(length(OUTSTRUCT.bnum), 2, 2);
CohScal_BaseWN = cell(length(OUTSTRUCT.bnum), 1);

for i=1:length(OUTSTRUCT.bnum)
    
    indsbase = OUTSTRUCT.indsbase_epoch{i};
    indswn = OUTSTRUCT.indsWN_epoch{i};
    
    cohbase = mean(OUTSTRUCT.cohscal{i}(indsbase));
    cohwn = mean(OUTSTRUCT.cohscal{i}(indswn));
    
    % =============== JUST AS AHACK, DO THIS (i.e. have 2x2 instead of 1x2 outoput.
    % SO THAT FOLLOWING CODE WORKS (PREBIOUSLY
    % WRITTEN FOR XCOV SCALARA, WHICH HAS 2 WINDOWS OF INTERST).
    tmp = nan(2,2);
    tmp(1,1) = cohbase;
    tmp(1,2) = cohwn;
    
    tmp(2,1) = cohbase;
    tmp(2,2) = cohwn;
    
    CohScal_BaseWN{i} = tmp;
end
OUTSTRUCT.CohScal_BaseWN = CohScal_BaseWN;


%% ========== COLLECT DATA FOR PAIRS OF UNITS

% === FOR EVERY NEURON PAIR, GET FR RHO FOR BASE AND WN, FOR UNIT 1 AND
% UNIT 2

FrRho_LMANRA_BaseWN = cell(length(OUTSTRUCT.bnum),1);

for i=1:length(OUTSTRUCT.bnum)
   bb = OUTSTRUCT.bnum(i);
   ee = OUTSTRUCT.enum(i);
   ss = OUTSTRUCT.switch(i);
   mm = OUTSTRUCT.motifnum(i);
   chanpair = OUTSTRUCT.chanpair(i,:);
   
   % ========== skip if this doesn['t have unit dat
   indthis = OUTSTRUCT_lfp.bnum==bb & OUTSTRUCT_lfp.enum==ee & ...
       OUTSTRUCT_lfp.switch==ss & OUTSTRUCT_lfp.motifnum==mm;
   if ~any(indthis)
       continue
   end
   
   % ============= ASSUME THAT FIRS NERUON IS LMAN, SECOND IS RA
   indstmp = OUTSTRUCT_lfp.bnum==bb & OUTSTRUCT_lfp.enum==ee & ...
       OUTSTRUCT_lfp.switch==ss & OUTSTRUCT_lfp.motifnum==mm & OUTSTRUCT_lfp.chanthis==chanpair(1);
   assert(strcmp(OUTSTRUCT_lfp.bregion{indstmp}, 'LMAN'));
   indstmp = OUTSTRUCT_lfp.bnum==bb & OUTSTRUCT_lfp.enum==ee & ...
       OUTSTRUCT_lfp.switch==ss & OUTSTRUCT_lfp.motifnum==mm & OUTSTRUCT_lfp.chanthis==chanpair(2);
   assert(strcmp(OUTSTRUCT_lfp.bregion{indstmp}, 'RA'));
   
   
   rhomeansall = nan(2,2); % neuron(2) x base/wn(2)
   
   % ================ neuron 1
   nthis = chanpair(1);
   
   indthis = OUTSTRUCT_lfp.bnum==bb & OUTSTRUCT_lfp.enum==ee & ...
       OUTSTRUCT_lfp.switch==ss & OUTSTRUCT_lfp.motifnum==mm ...
       & OUTSTRUCT_lfp.chanthis==nthis;
    assert(sum(indthis)==1);

    rho_means = OUTSTRUCT_lfp.xtrialFrRho_Means(indthis,:);
    rhomeansall(1,:) = rho_means;
    
   % ================ neuron 1
   nthis = chanpair(2);
   
   indthis = OUTSTRUCT_lfp.bnum==bb & OUTSTRUCT_lfp.enum==ee & ...
       OUTSTRUCT_lfp.switch==ss & OUTSTRUCT_lfp.motifnum==mm ...
       & OUTSTRUCT_lfp.chanthis==nthis;
    assert(sum(indthis)==1);

    rho_means = OUTSTRUCT_lfp.xtrialFrRho_Means(indthis,:);
    rhomeansall(2,:) = rho_means;
    
    % ====== OUTPUT
    FrRho_LMANRA_BaseWN{i} = rhomeansall;
   
end

OUTSTRUCT.FrRho_LMANRA_BaseWN = FrRho_LMANRA_BaseWN;


% ======================= ONLY KEEP CASES THAT HAVE DATA FOR UNITS
indstokeep = ~cellfun('isempty', OUTSTRUCT.FrRho_LMANRA_BaseWN);
OUTSTRUCT = lt_structure_RmvEmptyField(OUTSTRUCT);
OUTSTRUCT = lt_structure_subsample_all_fields(OUTSTRUCT, indstokeep, 1);

% ======================= GET GRP STATS
if strcmp(plotlevel, 'chanpair')
    [allbnum_pairs, allenum_pairs, allswnum_pairs, allFrRho_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'FrRho_LMANRA_BaseWN');
    [allbnum_pairs, allenum_pairs, allswnum_pairs, allXcovScal_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'CohScal_BaseWN');
elseif strcmp(plotlevel, 'switch')
    [~,~,~,~,allbnum_pairs, allenum_pairs, allswnum_pairs, allFrRho_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'FrRho_LMANRA_BaseWN');
    [~,~,~,~,allbnum_pairs, allenum_pairs, allswnum_pairs, allXcovScal_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT, 'CohScal_BaseWN');
end


%% =========== [PLOTS - comparing pairs to units] 
% =========== e.g. greater increae in correlation predicts more consistent
% LMAN/RA activity across trials?

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ###################################### TARG
% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [LMAN]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;
sylind = 1;

bregionind = 1; % 1= lmnan; 2 = ra

x = squeeze(allXcovScal_pairs(xcovwind, 2, sylind, :) - allXcovScal_pairs(xcovwind, 1, sylind, :));
y = squeeze(allFrRho_pairs(bregionind, 2, sylind, :) - allFrRho_pairs(bregionind, 1, sylind, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [RA]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;
sylind = 1;

bregionind = 2; % 1= lmnan; 2 = ra

x = squeeze(allXcovScal_pairs(xcovwind, 2, sylind, :) - allXcovScal_pairs(xcovwind, 1, sylind, :));
y = squeeze(allFrRho_pairs(bregionind, 2, sylind, :) - allFrRho_pairs(bregionind, 1, sylind, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= LMAN/RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [mean of LMAN/RA]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;
sylind = 1;

x = squeeze(allXcovScal_pairs(xcovwind, 2, sylind, :) - allXcovScal_pairs(xcovwind, 1, sylind, :));
y1 = squeeze(allFrRho_pairs(1, 2, sylind, :) - allFrRho_pairs(1, 1, sylind, :));
y2 = squeeze(allFrRho_pairs(2, 2, sylind, :) - allFrRho_pairs(2, 1, sylind, :));
y = mean([y1 y2],2);
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ###################################### TARG
% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [LMAN]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;
sylind = 1;

bregionind = 1; % 1= lmnan; 2 = ra

x = squeeze(allXcovScal_pairs(xcovwind, 2, sylind, :) - allXcovScal_pairs(xcovwind, 1, sylind, :));
y = squeeze(allFrRho_pairs(bregionind, 2, sylind, :) - allFrRho_pairs(bregionind, 1, sylind, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [RA]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;
sylind = 1;

bregionind = 2; % 1= lmnan; 2 = ra

x = squeeze(allXcovScal_pairs(xcovwind, 2, sylind, :) - allXcovScal_pairs(xcovwind, 1, sylind, :));
y = squeeze(allFrRho_pairs(bregionind, 2, sylind, :) - allFrRho_pairs(bregionind, 1, sylind, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= LMAN/RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [mean of LMAN/RA]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;
sylind = 1;

x = squeeze(allXcovScal_pairs(xcovwind, 2, sylind, :) - allXcovScal_pairs(xcovwind, 1, sylind, :));
y1 = squeeze(allFrRho_pairs(1, 2, sylind, :) - allFrRho_pairs(1, 1, sylind, :));
y2 = squeeze(allFrRho_pairs(2, 2, sylind, :) - allFrRho_pairs(2, 1, sylind, :));
y = mean([y1 y2],2);
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ###################################### TARG [mean of wind1 and 2]
% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [LMAN]');
xlabel('change in neuron-neuron corr (mean (wind1/2))');
ylabel('change in trial-trial corr (Wn - base');
sylind = 1;
bregionind = 1; % 1= lmnan; 2 = ra

x1 = squeeze(allXcovScal_pairs(1, 2, sylind, :) - allXcovScal_pairs(1, 1, sylind, :));
x2 = squeeze(allXcovScal_pairs(2, 2, sylind, :) - allXcovScal_pairs(2, 1, sylind, :));
x = mean([x1 x2],2);
y = squeeze(allFrRho_pairs(bregionind, 2, sylind, :) - allFrRho_pairs(bregionind, 1, sylind, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [RA]');
xlabel('change in neuron-neuron corr (mean (wind1/2))');
ylabel('change in trial-trial corr (Wn - base');
sylind = 1;

bregionind = 2; % 1= lmnan; 2 = ra

x1 = squeeze(allXcovScal_pairs(1, 2, sylind, :) - allXcovScal_pairs(1, 1, sylind, :));
x2 = squeeze(allXcovScal_pairs(2, 2, sylind, :) - allXcovScal_pairs(2, 1, sylind, :));
x = mean([x1 x2],2);
y = squeeze(allFrRho_pairs(bregionind, 2, sylind, :) - allFrRho_pairs(bregionind, 1, sylind, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= LMAN/RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG [mean of LMAN/RA]');
xlabel('change in neuron-neuron corr (mean (wind1/2))');
ylabel('change in trial-trial corr (Wn - base');
sylind = 1;

x1 = squeeze(allXcovScal_pairs(1, 2, sylind, :) - allXcovScal_pairs(1, 1, sylind, :));
x2 = squeeze(allXcovScal_pairs(2, 2, sylind, :) - allXcovScal_pairs(2, 1, sylind, :));
x = mean([x1 x2],2);
y1 = squeeze(allFrRho_pairs(1, 2, sylind, :) - allFrRho_pairs(1, 1, sylind, :));
y2 = squeeze(allFrRho_pairs(2, 2, sylind, :) - allFrRho_pairs(2, 1, sylind, :));
y = mean([y1 y2],2);
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



%% ================= NORMALIZE TO DIFFERENT TYPE SYKLS ...

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ###################################### TARG
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG-DIFF [LMAN]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;

bregionind = 1; % 1= lmnan; 2 = ra

x1 = squeeze(allXcovScal_pairs(xcovwind, 2, 1, :) - allXcovScal_pairs(xcovwind, 1, 1, :));
x2 = squeeze(allXcovScal_pairs(xcovwind, 2, 3, :) - allXcovScal_pairs(xcovwind, 1, 3, :));
x = x1-x2;

y1 = squeeze(allFrRho_pairs(bregionind, 2, 1, :) - allFrRho_pairs(bregionind, 1, 1, :));
y2 = squeeze(allFrRho_pairs(bregionind, 2, 3, :) - allFrRho_pairs(bregionind, 1, 3, :));
y = y1-y2;

% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG-DIFF [LMAN]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;

bregionind = 1; % 1= lmnan; 2 = ra

x1 = squeeze(allXcovScal_pairs(xcovwind, 2, 1, :) - allXcovScal_pairs(xcovwind, 1, 1, :));
x2 = squeeze(allXcovScal_pairs(xcovwind, 2, 3, :) - allXcovScal_pairs(xcovwind, 1, 3, :));
x = x1-x2;

y1 = squeeze(allFrRho_pairs(bregionind, 2, 1, :) - allFrRho_pairs(bregionind, 1, 1, :));
y2 = squeeze(allFrRho_pairs(bregionind, 2, 3, :) - allFrRho_pairs(bregionind, 1, 3, :));
y = y1-y2;

% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG-DIFF [RA]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;

bregionind = 2; % 1= lmnan; 2 = ra

x1 = squeeze(allXcovScal_pairs(xcovwind, 2, 1, :) - allXcovScal_pairs(xcovwind, 1, 1, :));
x2 = squeeze(allXcovScal_pairs(xcovwind, 2, 3, :) - allXcovScal_pairs(xcovwind, 1, 3, :));
x = x1-x2;

y1 = squeeze(allFrRho_pairs(bregionind, 2, 1, :) - allFrRho_pairs(bregionind, 1, 1, :));
y2 = squeeze(allFrRho_pairs(bregionind, 2, 3, :) - allFrRho_pairs(bregionind, 1, 3, :));
y = y1-y2;

% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);




% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG-DIFF [RA]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;

bregionind = 2; % 1= lmnan; 2 = ra

x1 = squeeze(allXcovScal_pairs(xcovwind, 2, 1, :) - allXcovScal_pairs(xcovwind, 1, 1, :));
x2 = squeeze(allXcovScal_pairs(xcovwind, 2, 3, :) - allXcovScal_pairs(xcovwind, 1, 3, :));
x = x1-x2;

y1 = squeeze(allFrRho_pairs(bregionind, 2, 1, :) - allFrRho_pairs(bregionind, 1, 1, :));
y2 = squeeze(allFrRho_pairs(bregionind, 2, 3, :) - allFrRho_pairs(bregionind, 1, 3, :));
y = y1-y2;

% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);


% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG-DIFF [LMAN]');
xlabel('change in neuron-neuron corr (mean(wind1,2))');
ylabel('change in trial-trial corr (Wn - base');
% xcovwind = 2;

bregionind = 1; % 1= lmnan; 2 = ra

x1 = squeeze(mean(allXcovScal_pairs([1 2], 2, 1, :),1) - mean(allXcovScal_pairs([1 2], 1, 1, :),1));
x2 = squeeze(mean(allXcovScal_pairs([1 2], 2, 3, :),1) - mean(allXcovScal_pairs([1 2], 1, 3, :),1));
% x2 = squeeze(allXcovScal_pairs(xcovwind, 2, 3, :) - allXcovScal_pairs(xcovwind, 1, 3, :));
x = x1-x2;

y1 = squeeze(allFrRho_pairs(bregionind, 2, 1, :) - allFrRho_pairs(bregionind, 1, 1, :));
y2 = squeeze(allFrRho_pairs(bregionind, 2, 3, :) - allFrRho_pairs(bregionind, 1, 3, :));
y = y1-y2;

% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG-DIFF [RA]');
xlabel('change in neuron-neuron corr (mean(wind1,2))');
ylabel('change in trial-trial corr (Wn - base');
% xcovwind = 2;

bregionind = 2; % 1= lmnan; 2 = ra

x1 = squeeze(mean(allXcovScal_pairs([1 2], 2, 1, :),1) - mean(allXcovScal_pairs([1 2], 1, 1, :),1));
x2 = squeeze(mean(allXcovScal_pairs([1 2], 2, 3, :),1) - mean(allXcovScal_pairs([1 2], 1, 3, :),1));
% x2 = squeeze(allXcovScal_pairs(xcovwind, 2, 3, :) - allXcovScal_pairs(xcovwind, 1, 3, :));
x = x1-x2;

y1 = squeeze(allFrRho_pairs(bregionind, 2, 1, :) - allFrRho_pairs(bregionind, 1, 1, :));
y2 = squeeze(allFrRho_pairs(bregionind, 2, 3, :) - allFrRho_pairs(bregionind, 1, 3, :));
y = y1-y2;

% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);

%% ================ COMBINE ALL SYLS/CASES
AllSyl_allXcovScal_pairs = lt_neural_Coher_Cell2Mat(OUTSTRUCT.CohScal_BaseWN);
AllSyl_allFrRho_pairs = lt_neural_Coher_Cell2Mat(OUTSTRUCT.FrRho_LMANRA_BaseWN);


% ###################################### TARG
% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('ALL SYLS [LMAN]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;

bregionind = 1; % 1= lmnan; 2 = ra

x = squeeze(AllSyl_allXcovScal_pairs(xcovwind, 2, :) - AllSyl_allXcovScal_pairs(xcovwind, 1, :));
y = squeeze(AllSyl_allFrRho_pairs(bregionind, 2, :) - AllSyl_allFrRho_pairs(bregionind, 1, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ========================= LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('ALL SYLS [LMAN]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;

bregionind = 1; % 1= lmnan; 2 = ra

x = squeeze(AllSyl_allXcovScal_pairs(xcovwind, 2, :) - AllSyl_allXcovScal_pairs(xcovwind, 1, :));
y = squeeze(AllSyl_allFrRho_pairs(bregionind, 2, :) - AllSyl_allFrRho_pairs(bregionind, 1, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ========================= RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('ALL SYLS [RA]');
xlabel('change in neuron-neuron corr (window 1)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 1;

bregionind = 2; % 1= lmnan; 2 = ra

x = squeeze(AllSyl_allXcovScal_pairs(xcovwind, 2, :) - AllSyl_allXcovScal_pairs(xcovwind, 1, :));
y = squeeze(AllSyl_allFrRho_pairs(bregionind, 2, :) - AllSyl_allFrRho_pairs(bregionind, 1, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);



% ========================= RA
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('ALL SYLS [RA]');
xlabel('change in neuron-neuron corr (window 2)');
ylabel('change in trial-trial corr (Wn - base');
xcovwind = 2;

bregionind = 2; % 1= lmnan; 2 = ra

x = squeeze(AllSyl_allXcovScal_pairs(xcovwind, 2, :) - AllSyl_allXcovScal_pairs(xcovwind, 1, :));
y = squeeze(AllSyl_allFrRho_pairs(bregionind, 2, :) - AllSyl_allFrRho_pairs(bregionind, 1, :));
% plot(x,y, 'ok');
lt_regress(y, x, 1, 0, 1, 1);

%% =========== [PLOTS] 

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG, RA');
xlabel('BASE (FR corr, x-trials)');
ylabel('WN');
lt_plot_annotation(1, ['corr wind: ' num2str(corrwind)], 'b');

x = squeeze(allDat_RA(1, 1, 1, :));
y = squeeze(allDat_RA(1, 2, 1, :));
lt_regress(y, x, 1, 0, 1, 1, 'r');
lt_plot_makesquare_plot45line(gca, 'k');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG, LMAN');
xlabel('BASE (FR corr, x-trials)');
ylabel('WN');

x = squeeze(allDat_LMAN(1, 1, 1, :));
y = squeeze(allDat_LMAN(1, 2, 1, :));
lt_regress(y, x, 1, 0, 1, 1, 'b');
lt_plot_makesquare_plot45line(gca, 'k');

% ======================== SAME
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('SAME, RA');
xlabel('BASE (FR corr, x-trials)');
ylabel('WN');

x = squeeze(allDat_RA(1, 1, 2, :));
y = squeeze(allDat_RA(1, 2, 2, :));
if any(~isnan(y))
lt_regress(y, x, 1, 0, 1, 1, 'r');
lt_plot_makesquare_plot45line(gca, 'k');
end

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('SAME, LMAN');
xlabel('BASE (FR corr, x-trials)');
ylabel('WN');

x = squeeze(allDat_LMAN(1, 1, 2, :));
y = squeeze(allDat_LMAN(1, 2, 2, :));
if any(~isnan(y))    lt_regress(y, x, 1, 0, 1, 1, 'b');
lt_plot_makesquare_plot45line(gca, 'k');
end

% ======================== DIFF
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('DIFF, RA');
xlabel('BASE (FR corr, x-trials)');
ylabel('WN');

x = squeeze(allDat_RA(1, 1, 3, :));
y = squeeze(allDat_RA(1, 2, 3, :));
lt_regress(y, x, 1, 0, 1, 1, 'r');
lt_plot_makesquare_plot45line(gca, 'k');


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('DIFF, LMAN');
xlabel('BASE (FR corr, x-trials)');
ylabel('WN');

x = squeeze(allDat_LMAN(1, 1, 3, :));
y = squeeze(allDat_LMAN(1, 2, 3, :));
lt_regress(y, x, 1, 0, 1, 1, 'b');
lt_plot_makesquare_plot45line(gca, 'k');



% ====== link
linkaxes(hsplots, 'xy');


%% ========== [PLOT] BASE-WN, LINE PLOTS
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG, RA');
ylabel('FR corr, x-trials');
xlabel('BASE-WN');

y1 = squeeze(allDat_RA(1, 1, 1, :));
y2 = squeeze(allDat_RA(1, 2, 1, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if length(y)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('TARG, LMAN');
ylabel('FR corr, x-trials');
xlabel('BASE-WN');

y1 = squeeze(allDat_LMAN(1, 1, 1, :));
y2 = squeeze(allDat_LMAN(1, 2, 1, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if length(y)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('SAME, RA');
ylabel('FR corr, x-trials');
xlabel('BASE-WN');

y1 = squeeze(allDat_RA(1, 1, 2, :));
y2 = squeeze(allDat_RA(1, 2, 2, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if length(y)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('SAME, LMAN');
ylabel('FR corr, x-trials');
xlabel('BASE-WN');

y1 = squeeze(allDat_LMAN(1, 1, 2, :));
y2 = squeeze(allDat_LMAN(1, 2, 2, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if length(y)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;




[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('DIFF, RA');
ylabel('FR corr, x-trials');
xlabel('BASE-WN');

y1 = squeeze(allDat_RA(1, 1, 3, :));
y2 = squeeze(allDat_RA(1, 2, 3, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if length(y)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('DIFF, LMAN');
ylabel('FR corr, x-trials');
xlabel('BASE-WN');

y1 = squeeze(allDat_LMAN(1, 1, 3, :));
y2 = squeeze(allDat_LMAN(1, 2, 3, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if length(y)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


%% ===== [PLOT] different format
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('RA, BASE');
ylabel('FR corr, x-trials');
xlabel('TARG - DIFF');

y1 = squeeze(allDat_RA(:, 1, 1, :));
y2 = squeeze(allDat_RA(:, 1, 3, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('RA, WN');
ylabel('FR corr, x-trials');
xlabel('TARG - DIFF');

y1 = squeeze(allDat_RA(:, 2, 1, :));
y2 = squeeze(allDat_RA(:, 2, 3, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;



[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LMAN, BASE');
ylabel('FR corr, x-trials');
xlabel('TARG - DIFF');

y1 = squeeze(allDat_LMAN(:, 1, 1, :));
y2 = squeeze(allDat_LMAN(:, 1, 3, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LMAN, WN');
ylabel('FR corr, x-trials');
xlabel('TARG - DIFF');

y1 = squeeze(allDat_LMAN(:, 2, 1, :));
y2 = squeeze(allDat_LMAN(:, 2, 3, :));
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;

% =========================== DIFFERENCES (TARG - DIFF)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('RA');
ylabel('FR corr, x-trials');
xlabel('corr, targ-diff [BASE, WN]');

ybase = squeeze(allDat_RA(:, 1, 1, :) - allDat_RA(:, 1, 3, :));
ywn = squeeze(allDat_RA(:, 2, 1, :) - allDat_RA(:, 2, 3, :));

Y = [ybase ywn];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LMAN');
ylabel('FR corr, x-trials');
xlabel('corr, targ-diff [BASE, WN]');

ybase = squeeze(allDat_LMAN(:, 1, 1, :) - allDat_LMAN(:, 1, 3, :));
ywn = squeeze(allDat_LMAN(:, 2, 1, :) - allDat_LMAN(:, 2, 3, :));

Y = [ybase ywn];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


% ######################## COMBINE RA AND LMAN
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LMAN+RA, BASE');
ylabel('FR corr, x-trials');
xlabel('TARG - DIFF');

y1 = [squeeze(allDat_RA(:, 1, 1, :)); squeeze(allDat_LMAN(:, 1, 1, :))];
y2 = [squeeze(allDat_RA(:, 1, 3, :)); squeeze(allDat_LMAN(:, 1, 3, :))];
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LMAN+RA, WN');
ylabel('FR corr, x-trials');
xlabel('TARG - DIFF');

y1 = [squeeze(allDat_RA(:, 2, 1, :)); squeeze(allDat_LMAN(:, 2, 1, :))];
y2 = [squeeze(allDat_RA(:, 2, 3, :)); squeeze(allDat_LMAN(:, 2, 3, :))];
Y = [y1 y2];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;


% =========================== DIFFERENCES (TARG - DIFF)
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('LMAN + RA');
ylabel('FR corr, x-trials');
xlabel('corr, targ-diff [BASE, WN]');

ybase = [squeeze(allDat_RA(:, 1, 1, :) - allDat_RA(:, 1, 3, :)); ...
    squeeze(allDat_LMAN(:, 1, 1, :) - allDat_LMAN(:, 1, 3, :))];
ywn = [squeeze(allDat_RA(:, 2, 1, :) - allDat_RA(:, 2, 3, :)); ...
    squeeze(allDat_LMAN(:, 2, 1, :) - allDat_LMAN(:, 2, 3, :))];

Y = [ybase ywn];
x = 1:size(Y,2);
plot(x, Y', '-ok');
if size(Y,1)>1
lt_plot(x+0.15, nanmean(Y), {'Errors', lt_sem(Y), 'Color', 'r', 'LineStyle', '-'});
p = signrank(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'srank');
end
xlim([0 3]);
lt_plot_zeroline;

