function lt_neural_LFP_AlignToWN(OUTSTRUCT, SwitchCohStruct, MOTIFSTATS_pop, ...
    SwitchStruct, PARAMS, prctiletouse, dozscore, prewind_relWN, ...
    onlyusegoodtargsyls)
%% lt 1/4/19 - coherence change better aligned to WN onset of syl onset?

% prctiletouse = 50;

% === for plotting
clim = [-0.15 0.15];
ffbintoget = [20 35];
% dozscore = 0; % 0=no, 1=yes, 2=center onoy

%%

%% ======= 2) FILTER DATA
% == filter by specific types of switches.
OUTSTRUCT = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct);

%% ======== 1) FOR EACH SWITCH, EXTRACT TARGET COHERENCE MATRIX.
% ============
%     fieldtoget = 'CohMean_WNminusBase';
%     [~, ~, ~, ~, allbnum1, allenum, allswnum, allDat] = ...
%         lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

fieldtoget = 'CohMean_WN';
[~, ~, ~, ~, allbnum, allenum, allswnum, allDat2] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

fieldtoget = 'CohMean_Base';
[~, ~, ~, ~, allbnum, allenum, allswnum, allDat1] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget);

% ---- subtract WN from base
allcohdiff = nan(size(allDat1));
for i=1:size(allcohdiff,3)
    for ii=1:size(allcohdiff,4)
        
        allcohdiff(:,:,i,ii) = allDat2(:,:,i, ii) - allDat1(:,:,i, ii);
    end
end



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
    
    motiflist = unique(OUTSTRUCT.motifnum(indstmp));
    bname = SwitchStruct.bird(allbnum(i)).birdname;
    ename = SwitchStruct.bird(allbnum(i)).exptnum(allenum(i)).exptname;
    sw = allswnum(i);
    
    % ===
    all_pctiles = [];
    for mm=1:length(motiflist)
        dat = SwitchCohStruct.bird(allbnum(i)).exptnum(allenum(i)).switchlist(allswnum(i)).motifnum(motiflist(mm));
        
        if onlyusegoodtargsyls==1
            sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, sw, dat.motifname);
            if sylbad==1
                continue
            end
        end
        
        wntimes = dat.WNhittimes_min - PARAMS.motif_predur;
        indsgood = [dat.indsbase_epoch dat.indsWN_epoch];
        
        % === only look at WN times within baseline or WN epochs
        wntimes = wntimes(indsgood);
        wntimes = wntimes(~isnan(wntimes)); % only look at hit trials
        
        % ============ GET PERCENTILES FOR TIMES
        all_pctiles = [all_pctiles; prctile(wntimes, prctiletouse)];
    end
    
    % ========= IOUTPUT
    allwntimes = [allwntimes; mean(all_pctiles)];
end

%% ======================= sort everything by order of wn time
[~, indsort] = sort(allwntimes);

allwntimes = allwntimes(indsort);
allbnum = allbnum(indsort);
allenum = allenum(indsort);
allswnum = allswnum(indsort);
allcohdiff = allcohdiff(:,:,:,indsort);
%% ############################ [PLOTS]


%% ================ 1) plot all coherograms, sorted by WN onset
figcount=1;
subplotrows=2;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(allbnum)
    
    bname = SwitchStruct.bird(allbnum(i)).birdname;
    ename = SwitchStruct.bird(allbnum(i)).exptnum(allenum(i)).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([bname '-' ename '-sw' num2str(allswnum(i))]);
    
    % ==== target minus diff
    cohdiff = allcohdiff(:,:,:,i);
    coh = cohdiff(:,:,1) - cohdiff(:,:,3);
    lt_neural_Coher_Plot(coh, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
    ylim([20 100]);
    
    % === note down time of min WN
    line([allwntimes(i) allwntimes(i)], ylim, 'Color', [0.2 0.8 0.5], 'LineWidth', 2)
end



%% ================ 2) smoothed coherence in F bin [PLOT EACH EXPT]
figcount=1;
subplotrows=2;
subplotcols=6;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

for i=1:length(allbnum)
    
    bname = SwitchStruct.bird(allbnum(i)).birdname;
    ename = SwitchStruct.bird(allbnum(i)).exptnum(allenum(i)).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title([bname '-' ename '-sw' num2str(allswnum(i))]);
    
    % ==== target minus diff
    cohdiff = allcohdiff(:,:,:,i);
    coh = cohdiff(:,:,1) - cohdiff(:,:,3);
    %     lt_neural_Coher_Plot(coh, PARAMS.tbins, PARAMS.ffbins, 1, '', clim);
    lt_neural_Coher_Plot(coh, PARAMS.tbins, PARAMS.ffbins, 2, '-', clim, 1, 0, PARAMS.ffbinsedges);
    %     ylim([20 100]);
    
    % === note down time of min WN
    line([allwntimes(i) allwntimes(i)], ylim, 'Color', [0.2 0.8 0.5], 'LineWidth', 2);
    lt_plot_zeroline;
    
end

linkaxes(hsplots, 'x');


%% ================ 2) smoothed coherence in F bin [ONE WATERFAL PLOT]
%  ====== ASSUMES WANT FIRST

indf = PARAMS.ffbins>ffbintoget(1) & PARAMS.ffbins<ffbintoget(2);

allcohTimecourse = [];
lt_figure; hold on;
title('coherence, cross is time of WN');
ylabel('experiment');

ythisall = [];
for i=1:length(allbnum)
    
    bname = SwitchStruct.bird(allbnum(i)).birdname;
    ename = SwitchStruct.bird(allbnum(i)).exptnum(allenum(i)).exptname;
    
    % ==== target minus diff
    cohdiff = allcohdiff(:,:,:,i);
    coh = cohdiff(:,:,1) - cohdiff(:,:,3);
    cohvec = mean(coh(:, indf),2);
    if dozscore==1
        cohvec = (cohvec-mean(cohvec))./std(cohvec);
        ythis  = 3*i-3;
    elseif dozscore==2
        cohvec = (cohvec-mean(cohvec));
        ythis  = 0.2*i-0.2;
    else
        ythis  = 0.2*i-0.2;
    end
    
    % ==== PLOT
    plot(PARAMS.tbins, cohvec+ythis, '-k');
    % --- highlight positiv eand negative parts using line plot
    thresh = 0.05;
    plot(PARAMS.tbins(cohvec>thresh), cohvec(cohvec>thresh)+ythis, 'or', 'LineWidth', 2);
    line([PARAMS.tbins(1) PARAMS.tbins(end)], [ythis ythis], 'Color', [0.3 0.3 0.3], 'LineStyle', ':');
    
    % === plot marker indicating time of onset
    plot(allwntimes(i), ythis, 'xb', 'MarkerSize', 8);
    
    % ---- plot time of premotor window
    tmp = allwntimes(i) + prewind_relWN;
    plot(tmp, ythis, 'sb');
    
    
    % == indicate bird and expt
    lt_plot_text(PARAMS.tbins(end), ythis, [bname '-' ename], 'm', 8);
    
    % ==== collect alL
    allcohTimecourse = [allcohTimecourse; cohvec'];
    ythisall = [ythisall; ythis];
end

axis tight;
lt_plot_zeroline_vert;

% == line for traditional premotior window
line([-0.07 -0.07], ylim, 'Color', 'b', 'LineStyle', '--');
line([-0.03 -0.03], ylim, 'Color', 'b', 'LineStyle', '--');

% ===== plot as heat map
lt_figure; hold on;
imagesc(PARAMS.tbins, ythisall, allcohTimecourse, [-0.1 0.1]);
lt_plot_colormap('centered');
colorbar('EastOutside');

plot(allwntimes, ythisall, 'xk', 'MarkerSize', 8);
axis tight;
lt_plot_zeroline_vert;



%% ========== take mean timecourse for first and second median split of data
lt_figure; hold on;

% ---- short latency
indtmp = allwntimes < median(allwntimes);
cohmat = squeeze(allcohdiff(:,:,1,indtmp) - allcohdiff(:,:,3, indtmp));

tmin = min(allwntimes(indtmp));
tmax = max(allwntimes(indtmp));
tmed = median(allwntimes(indtmp));

lt_subplot(3,2,1); hold on;
title('short latency (median split)');
xlabel('[magenta: min, med, max]');
lt_neural_Coher_Plot(cohmat, PARAMS.tbins, PARAMS.ffbins, 1, '', clim, 0, 0, PARAMS.ffbinsedges);
line([tmin tmin], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], ylim, 'LineWidth', 2, 'Color', 'm');

lt_subplot(3,2,2); hold on;
title('short latency (median split)');
lt_neural_Coher_Plot(cohmat, PARAMS.tbins, PARAMS.ffbins, 2, '-', clim, 1, 0, PARAMS.ffbinsedges);
line([tmin tmin], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], ylim, 'LineWidth', 2, 'Color', 'm');


% ---- short latency
indtmp = allwntimes >= median(allwntimes);
cohmat = squeeze(allcohdiff(:,:,1,indtmp) - allcohdiff(:,:,3, indtmp));

tmin = min(allwntimes(indtmp));
tmax = max(allwntimes(indtmp));
tmed = median(allwntimes(indtmp));

lt_subplot(3,2,3); hold on;
title('long latency (median split)');
lt_neural_Coher_Plot(cohmat, PARAMS.tbins, PARAMS.ffbins, 1, '', clim, 0, 0, PARAMS.ffbinsedges);
line([tmin tmin], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], ylim, 'LineWidth', 2, 'Color', 'm');

lt_subplot(3,2,4); hold on;
title('long latency (median split)');
lt_neural_Coher_Plot(cohmat, PARAMS.tbins, PARAMS.ffbins, 2, '-', clim, 1, 0, PARAMS.ffbinsedges);
line([tmin tmin], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmax tmax], ylim, 'LineWidth', 1, 'Color', 'm');
line([tmed tmed], ylim, 'LineWidth', 2, 'Color', 'm');

