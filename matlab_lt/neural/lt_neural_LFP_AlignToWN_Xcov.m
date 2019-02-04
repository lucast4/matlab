function lt_neural_LFP_AlignToWN_Xcov(OUTSTRUCT_XCOV, OUTSTRUCT, SwitchStruct, ...
    PARAMS, dozscore, prewind_relWN)
%% lt 1/4/19 - coherence change better aligned to WN onset of syl onset?

% prctiletouse = 50;

% === for plotting
clim = [-0.15 0.15];
% ffbintoget = [20 35];
% dozscore = 0; % 0=no, 1=yes, 2=center onoy
prctiletouse = 2.5;

%%

%% ======= 2) FILTER DATA
% == filter by specific types of switches.
dataset = 'xcov_spikes';
OUTSTRUCT = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct, ...
    dataset);

OUTSTRUCT_XCOV = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT_XCOV, SwitchStruct, ...
    dataset);

%% ======== only keep if LFP


%% ======== 1) FOR EACH SWITCH, EXTRACT TARGET HEAT MAP OF XCOV MATRIX.
% ============
%     fieldtoget = 'CohMean_WNminusBase';
%     [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
%         lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

fieldtoget = 'XcovgramBase';
[~, ~, ~, ~, allbnum, allenum, allswnum, allDat1] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);

fieldtoget = 'XcovgramWN';
[~, ~, ~, ~, allbnum, allenum, allswnum, allDat2] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT_XCOV, fieldtoget);

% Get difference (i.e. change from baseline)
allDatDiff = allDat2 - allDat1;


%% ======== GET TIMING OF WN ONSETS
% ======== for each switch, figure out timing for this WN file
% --- if multiple targets, then will take the mean of value for each
% target.

allwntimes = []; % 2.5 PERCENTILEs
% all_midtimes = []; % m edian
for i=1:length(allbnum)
    
    % === what are target syls?
    indstmp = OUTSTRUCT.bnum==allbnum(i) & OUTSTRUCT.enum==allenum(i) & OUTSTRUCT.switch==allswnum(i) ...
        & OUTSTRUCT.istarg==1;
    
    wntime = unique(OUTSTRUCT.WNonset(indstmp));
    

    % ========= IOUTPUT
    allwntimes = [allwntimes; wntime];
end

%% ======================= sort everything by order of wn time
[~, indsort] = sort(allwntimes);

allwntimes = allwntimes(indsort);
allbnum = allbnum(indsort);
allenum = allenum(indsort);
allswnum = allswnum(indsort);
allDatDiff = allDatDiff(:,:,:, indsort);
% 
% allcohdiff = allcohdiff(:,:,:,indsort);

%% ############################ [PLOTS]
%% ================ 1) plot all xcovgrams, sorted by WN onset
figcount=1;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(allbnum)
    
    bname = SwitchStruct.bird(allbnum(i)).birdname;
    ename = SwitchStruct.bird(allbnum(i)).exptnum(allenum(i)).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(allswnum(i))]);
    ylabel('Targ - diff (WN - base)');
    
    % ==== target minus diff
    cohdiff = allDatDiff(:,:,:,i);
    Y = cohdiff(:,:,1) - cohdiff(:,:,3);
    lt_neural_Coher_Plot(Y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
%     ylim([20 100]);
    
    % === note down time of min WN
    line([allwntimes(i) allwntimes(i)], ylim, 'Color', [0.2 0.8 0.5], 'LineWidth', 2);
    lt_plot_zeroline;    
end

%% ================ 1) plot all xcovgrams, sorted by WN onset
figcount=1;
subplotrows=2;
subplotcols=4;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(allbnum)
    
    bname = SwitchStruct.bird(allbnum(i)).birdname;
    ename = SwitchStruct.bird(allbnum(i)).exptnum(allenum(i)).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(allswnum(i))]);
    ylabel('Targ (WN - base)');
    
    % ==== target minus diff
    cohdiff = allDatDiff(:,:,:,i);
    Y = cohdiff(:,:,1);
    lt_neural_Coher_Plot(Y, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '', clim);
%     ylim([20 100]);
    
    % === note down time of min WN
    line([allwntimes(i) allwntimes(i)], ylim, 'Color', [0.2 0.8 0.5], 'LineWidth', 2);
    lt_plot_zeroline;    
end



%% =============== TAKE AVERAGE BY MEDIAN SPLIT
lt_figure; hold on;

% ============================= short latency
indtmp = allwntimes <= median(allwntimes);

tmin = min(allwntimes(indtmp));
tmax = max(allwntimes(indtmp));
tmed = median(allwntimes(indtmp));

% --- TARG - DIFF
lt_subplot(2,3,1); hold on;
title('short latency (median split)');
cohmat = squeeze(allDatDiff(:,:,1,indtmp) - allDatDiff(:,:,3, indtmp));
ylabel('TARG - DIFF (WN - base)');
lt_neural_Coher_Plot(cohmat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '-', clim, 1, 0);
YLIM = ylim;
line([tmin tmin], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], YLIM, 'LineWidth', 2, 'Color', 'm');
lt_plot_zeroline;

% --- TARG - DIFF
lt_subplot(2,3,2); hold on;
title('short latency (median split)');
cohmat = squeeze(allDatDiff(:,:,1,indtmp) - allDatDiff(:,:,3, indtmp));
ylabel('TARG - DIFF (WN - base)');
lt_neural_Coher_Plot(cohmat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, '-', clim, 1, 0);
YLIM = ylim;
line([tmin tmin], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], YLIM, 'LineWidth', 2, 'Color', 'm');
lt_plot_zeroline;

% --- TARG
lt_subplot(2,3,3); hold on;
title('short latency (median split)');
cohmat = squeeze(allDatDiff(:,:,1,indtmp));
ylabel('TARG (WN - base)');
lt_neural_Coher_Plot(cohmat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '-', clim, 1, 0);
YLIM = ylim;
line([tmin tmin], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], YLIM, 'LineWidth', 2, 'Color', 'm');
lt_plot_zeroline;


% ============================= long latency
indtmp = allwntimes > median(allwntimes);

tmin = min(allwntimes(indtmp));
tmax = max(allwntimes(indtmp));
tmed = median(allwntimes(indtmp));

% --- TARG - DIFF
lt_subplot(2,3,4); hold on;
title('long latency (median split)');
cohmat = squeeze(allDatDiff(:,:,1,indtmp) - allDatDiff(:,:,3, indtmp));
ylabel('TARG - DIFF (WN - base)');
lt_neural_Coher_Plot(cohmat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '-', clim, 1, 0);
YLIM = ylim;
line([tmin tmin], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], YLIM, 'LineWidth', 2, 'Color', 'm');
lt_plot_zeroline;

% --- TARG - DIFF
lt_subplot(2,3,5); hold on;
title('long latency (median split)');
cohmat = squeeze(allDatDiff(:,:,1,indtmp) - allDatDiff(:,:,3, indtmp));
ylabel('TARG - DIFF (WN - base)');
lt_neural_Coher_Plot(cohmat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 3, '-', clim, 1, 0);
YLIM = ylim;
line([tmin tmin], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], YLIM, 'LineWidth', 2, 'Color', 'm');
lt_plot_zeroline;

% --- TARG
lt_subplot(2,3,6); hold on;
title('long latency (median split)');
cohmat = squeeze(allDatDiff(:,:,1,indtmp));
ylabel('TARG (WN - base)');
lt_neural_Coher_Plot(cohmat, PARAMS.xcenters_gram, PARAMS.Xcov_ccLags, 1, '-', clim, 1, 0);
YLIM = ylim;
line([tmin tmin], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], YLIM, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], YLIM, 'LineWidth', 2, 'Color', 'm');
lt_plot_zeroline;
