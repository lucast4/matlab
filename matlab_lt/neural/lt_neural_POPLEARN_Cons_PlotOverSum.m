function lt_neural_POPLEARN_Cons_PlotOverSum(OUTSTRUCT_units, OUTSTRUCT_XCOV, ...
    SwitchStruct, SummaryStruct, birdtoplot, onlygoodexpt, plotlevel, ...
    corrwind, expttoplot)
%% lt 2/5/19 - Overview plots for consistency of smoothed FR (x-trials)

% birdtoplot = []; % empty for all.

%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter outstruct
    expttype = 'xcov_spikes';
    [OUTSTRUCT_units] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_units, SwitchStruct, expttype);
end

%% ========= filter by bird/expt
if ~isempty(birdtoplot)
    indstokeep = ismember(OUTSTRUCT_units.bnum, birdtoplot);
    OUTSTRUCT_units = lt_structure_subsample_all_fields(OUTSTRUCT_units, indstokeep, 1);
end

if ~isempty(expttoplot)
    indstokeep = ismember(OUTSTRUCT_units.enum, expttoplot);
    OUTSTRUCT_units = lt_structure_subsample_all_fields(OUTSTRUCT_units, indstokeep, 1);
end

%% ========= GET MEAN RHO FOR EACH CASE

OUTSTRUCT_units.xtrialFrRho_Means = cellfun(@nanmean, OUTSTRUCT_units.xtrialFrRho_BaseWn);
OUTSTRUCT_units.chanpair = OUTSTRUCT_units.neurID; % just ad hoc for the code that expects this (GrpStats)

%% ========= SPLIT DATASET BASED ON BRAIN REGION

inds_RA = strcmp(OUTSTRUCT_units.bregion, 'RA');
OUTSTRUCT_units_RA = lt_structure_subsample_all_fields(OUTSTRUCT_units, inds_RA, 1);

inds_LMAN = strcmp(OUTSTRUCT_units.bregion, 'LMAN');
OUTSTRUCT_units_LMAN = lt_structure_subsample_all_fields(OUTSTRUCT_units, inds_LMAN, 1);

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




%% ========== COLLECT DATA FOR PAIRS OF UNITS

% === FOR EVERY NEURON PAIR, GET FR RHO FOR BASE AND WN, FOR UNIT 1 AND
% UNIT 2

FrRho_LMANRA_BaseWN = cell(length(OUTSTRUCT_XCOV.bnum),1);

for i=1:length(OUTSTRUCT_XCOV.bnum)
   bb = OUTSTRUCT_XCOV.bnum(i);
   ee = OUTSTRUCT_XCOV.enum(i);
   ss = OUTSTRUCT_XCOV.switch(i);
   mm = OUTSTRUCT_XCOV.motifnum(i);
   neurpair = OUTSTRUCT_XCOV.neurpair(i,:);
   
   % ========== skip if this doesn['t have unit dat
   indthis = OUTSTRUCT_units.bnum==bb & OUTSTRUCT_units.enum==ee & ...
       OUTSTRUCT_units.switch==ss & OUTSTRUCT_units.motifnum==mm;
   if ~any(indthis)
       continue
   end
   
   % ============= ASSUME THAT FIRS NERUON IS LMAN, SECOND IS RA
   assert(strcmp(SummaryStruct.birds(bb).neurons(neurpair(1)).NOTE_Location, 'LMAN'));
   assert(strcmp(SummaryStruct.birds(bb).neurons(neurpair(2)).NOTE_Location, 'RA'));
   
   
   rhomeansall = nan(2,2); % neuron(2) x base/wn(2)
   % ================ neuron 1
   nthis = neurpair(1);
   
   indthis = OUTSTRUCT_units.bnum==bb & OUTSTRUCT_units.enum==ee & ...
       OUTSTRUCT_units.switch==ss & OUTSTRUCT_units.motifnum==mm ...
       & OUTSTRUCT_units.neurID==nthis;
    assert(sum(indthis)==1);

    rho_means = OUTSTRUCT_units.xtrialFrRho_Means(indthis,:);
    rhomeansall(1,:) = rho_means;
    
   % ================ neuron 1
   nthis = neurpair(2);
   
   indthis = OUTSTRUCT_units.bnum==bb & OUTSTRUCT_units.enum==ee & ...
       OUTSTRUCT_units.switch==ss & OUTSTRUCT_units.motifnum==mm ...
       & OUTSTRUCT_units.neurID==nthis;
    assert(sum(indthis)==1);

    rho_means = OUTSTRUCT_units.xtrialFrRho_Means(indthis,:);
    rhomeansall(2,:) = rho_means;
    
    % ====== OUTPUT
    FrRho_LMANRA_BaseWN{i} = rhomeansall;
   
end

OUTSTRUCT_XCOV.FrRho_LMANRA_BaseWN = FrRho_LMANRA_BaseWN;


% ======================= ONLY KEEP CASES THAT HAVE DATA FOR UNITS
indstokeep = ~cellfun('isempty', OUTSTRUCT_XCOV.FrRho_LMANRA_BaseWN);
OUTSTRUCT_XCOV = lt_structure_subsample_all_fields(OUTSTRUCT_XCOV, indstokeep, 1);

% ======================= GET GRP STATS
if strcmp(plotlevel, 'chanpair')
    [allbnum_pairs, allenum_pairs, allswnum_pairs, allFrRho_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, 'FrRho_LMANRA_BaseWN');
    [allbnum_pairs, allenum_pairs, allswnum_pairs, allXcovScal_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, 'Xcovscal_window_BaseWN');
elseif strcmp(plotlevel, 'switch')
    [~,~,~,~,allbnum_pairs, allenum_pairs, allswnum_pairs, allFrRho_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, 'FrRho_LMANRA_BaseWN');
    [~,~,~,~,allbnum_pairs, allenum_pairs, allswnum_pairs, allXcovScal_pairs] = ...
        lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, 'Xcovscal_window_BaseWN');
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
AllSyl_allXcovScal_pairs = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.Xcovscal_window_BaseWN);
AllSyl_allFrRho_pairs = lt_neural_Coher_Cell2Mat(OUTSTRUCT_XCOV.FrRho_LMANRA_BaseWN);


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

