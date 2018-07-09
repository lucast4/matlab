function lt_neural_v2_CTXT_Acoustic_Compare(analyfnames, bregionstoplot, ...
    onlyGoodDecode, onlyIfTwoCtxts)
%% LT 5/23/18 - compares multiple (currently max is 2) datasets,
% interaction between context and FF.

% onlyGoodDecode = 1; %

%% LOAD ALL DATA

savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

ALLBRANCH = struct;

for j=1:length(analyfnames)
    
    sdir = [savedir '/' analyfnames{j} '/Acoustic_Corr'];
    
    ALLBRANCH(j).analyfname = analyfnames{j};
    ALLBRANCH(j).bregion = bregionstoplot{j};
    
    % ------- load all things
    fnames = dir([sdir '/AllBranch_*']);
    
    for jj=1:length(fnames)
        
        tmp = load([sdir '/' fnames(jj).name]);
        
        % ======== stick into output structure
        fieldnamethis = fieldnames(tmp);
        assert(length(fieldnamethis)==1, 'sdfasd');
        
        ALLBRANCH(j).DAT.(fieldnamethis{1}) = tmp.(fieldnamethis{1});
        
    end
end



%% ============= FILTER DATASET

for rr=1:length(ALLBRANCH)
    
    % -------------------------------------
    indsgood = logical(ones(size(ALLBRANCH(rr).DAT.AllBranch_birdnum)));
    
    if onlyGoodDecode==1
        indsgood = indsgood & ALLBRANCH(rr).DAT.AllBranch_DecodeP<0.05;
    elseif onlyGoodDecode==2
        % then only bad decode
        indsgood = indsgood & ALLBRANCH(rr).DAT.AllBranch_DecodeP>=0.05;
    end
    
    if onlyIfTwoCtxts==1
        indsgood = indsgood & ALLBRANCH(rr).DAT.AllBranch_Nctxt==2;
    end
    % -------------------------------------
    
    ALLBRANCH(rr).IndsGood = indsgood;
end

%% ============= SHOW DISTRIBUTIOSN

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



for rr=1:length(ALLBRANCH)
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(ALLBRANCH(rr).bregion);
    xlabel('slope(overall');
    ylabel('slope(interaction)');
    
    % -------------------------------------
    indsgood = ALLBRANCH(rr).IndsGood;
    % -------------------------------------
    
    
    slope_overall = ALLBRANCH(rr).DAT.AllBranch_SlopeOverallCoeff_MedAbs(indsgood);
    slope_overall_sig = ALLBRANCH(rr).DAT.AllBranch_SlopeOverall(indsgood);
    
    slope_diff = ALLBRANCH(rr).DAT.AllBranch_SlopeCoeff_MedAbs(indsgood);
    slope_diff_sig = ALLBRANCH(rr).DAT.AllBranch_SlopeDiff(indsgood);
    
    
    % --- none sig
    indtmp = slope_overall_sig==0 & slope_diff_sig==0;
    pcol = [0.6 0.6 0.6];
    plot(slope_overall(indtmp), slope_diff(indtmp), 'o', 'Color', pcol);
    
    % --- both sig
    indtmp = slope_overall_sig==1 & slope_diff_sig==1;
    pcol = 'm';
    plot(slope_overall(indtmp), slope_diff(indtmp), 'o', 'Color', pcol);
    
    % --- one sig
    indtmp = slope_overall_sig==1 & slope_diff_sig==0;
    pcol = 'b';
    plot(slope_overall(indtmp), slope_diff(indtmp), 'o', 'Color', pcol);
    
    
    % --- one sig
    indtmp = slope_overall_sig==0 & slope_diff_sig==1;
    pcol = 'r';
    plot(slope_overall(indtmp), slope_diff(indtmp), 'o', 'Color', pcol);
    
    % --- fraction with interaction greater than overall
    frac = sum(slope_overall<slope_diff)./length(slope_overall);
    lt_plot_annotation(1, ['frac(intera>overall)=' num2str(frac)], 'r');
    
    lt_plot_makesquare_plot45line(gca, 'k', -0.2);
    
    
    
    % ==================== Distributions (1d)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title('b=overall; r=interaction');
    % ----- SLOPES (OVERALL)
    pcol = 'b';
    [~, xbin] = lt_plot_histogram(slope_overall, '', 1, 0, '', 1, pcol);
    
    % ----- SLOPES (INTERA)
    pcol = 'r';
    lt_plot_histogram(slope_diff, xbin, 1, 0, '', 1, pcol);
    
    
    % ==================== FRACTIONS FOR EACH CLASS (SLOPE, INTERACTION ,BOTH)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(ALLBRANCH(rr).bregion);
    ylabel('num sig cases in each class');
    xlabel('none -- overallOnly -- InteracOnly -- Both');
    
    Y = [];
    % --- none sig
    indtmp = slope_overall_sig==0 & slope_diff_sig==0;
    Y = [Y sum(indtmp)];
    
    % --- one sig
    indtmp = slope_overall_sig==1 & slope_diff_sig==0;
    Y = [Y sum(indtmp)];
    
    % --- one sig
    indtmp = slope_overall_sig==0 & slope_diff_sig==1;
    Y = [Y sum(indtmp)];
    
    % --- both sig
    indtmp = slope_overall_sig==1 & slope_diff_sig==1;
    Y = [Y sum(indtmp)];
    
    % ------------ PLOT
    x =1:4;
    lt_plot_bar(x, Y);
    xlim([0 5]);
    
    
    % ====================== FOR EACH CASE PAIRED PLOT (OVERALL VS. INTERACT)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(ALLBRANCH(rr).bregion);
    ylabel('effect size (abs)');
    xlabel('Overall -- Interact');
    
    x = [1 2];
    Y = {slope_overall, slope_diff};
    
    % % --- none sig
    %     indtmp = slope_overall_sig==0 & slope_diff_sig==0;
    %     pcol = [0.6 0.6 0.6];
    %     plot(x, Y(indtmp,:), '-o', 'Color', pcol);
    %
    %     % --- both sig
    %     indtmp = slope_overall_sig==1 & slope_diff_sig==1;
    %     pcol = 'm';
    %     plot(x+0.1, Y(indtmp,:), '-o', 'Color', pcol);
    %
    %     % --- one sig
    %     indtmp = slope_overall_sig==1 & slope_diff_sig==0;
    %     pcol = 'b';
    %     plot(x+0.2, Y(indtmp,:), '-o', 'Color', pcol);
    %
    %
    %     % --- one sig
    %     indtmp = slope_overall_sig==0 & slope_diff_sig==1;
    %     pcol = 'r';
    %     plot(x+0.3, Y(indtmp,:), '-o', 'Color', pcol);
    
    lt_plot_MultDist(Y, x, 1, 'k', 0, 0, 1);
    plot(x, cell2mat(Y), '-');
    xlim([0 3]);
    
end

%% =============  QQ PLOT comparing 2 brain regions

datfield = 'AllBranch_SlopeCoeff_MedAbs';
shufffield = 'SlopeCoeff_MedAbs';
% datfield = 'AllBranch_SlopeCoeff_MedAbs';
% shufffield = 'SlopeCoeff_MedAbs';

PvalsAll = {};
for rr=1:length(ALLBRANCH)
    
    lt_figure; hold on;
    
    indsgood = ALLBRANCH(rr).IndsGood;
    
    % ================ EXTRACT COEFFIFIENTS
    ydat = ALLBRANCH(rr).DAT.(datfield);
    yshuff = [ALLBRANCH(rr).DAT.AllBranch_SHUFFSTRUCT.cycle.(shufffield)];
    
    ydat = ydat(indsgood,:);
    yshuff = yshuff(indsgood,:);
    
    
    % ================= get pvalues (empirical
    pvalmat = (1+sum(yshuff>=ydat, 2))./(1+size(yshuff,2));
    PvalsAll = [PvalsAll pvalmat];
    
    % ============ PLOT - fraction of cases signifncat
    ind_sigcases = pvalmat<0.05;
    disp([num2str(sum(ind_sigcases)) ' sig cases out of ' num2str(length(ind_sigcases))]);
    lt_subplot(3,2,5); hold on;
    title('num sig cases');
    lt_plot_bar(1, sum(ind_sigcases)/length(ind_sigcases));
    ylim([0 1]); xlim([0 2]);
    
    % =============== PLOT - score of all data rel shuffle
    lt_subplot(3,2,1); hold on;
    title('singal case shuff distr (rand choice)');
    casethis = randi(size(yshuff,1),1);
    lt_plot_histogram(yshuff(casethis,:));
    
    lt_subplot(3,2,2); hold on;
    title('all case zscore rel own shuffles');
    y_zscore = (ydat - mean(yshuff,2))./std(yshuff,0,2);
    
    lt_plot_histogram(y_zscore);
    lt_plot_zeroline_vert;
    
    % ================== PLOT, ALL P VALUES
    lt_subplot(3,2,3); hold on;
    title('all cases, percentiles rel own shuffle');
    xcenters = 0:0.05:1;
    [Ybinned] = lt_plot_histogram(pvalmat, xcenters);
    
    lt_subplot(3,2,4); hold on;
    title('all cases, percentiles rel own shuffle');
    xlabel('prct');
    ylabel('frac cases equal to or below this percent')
    
    x = xcenters;
    binsize = xcenters(2)-xcenters(1);
    x = x+binsize/2;
    x(end) = 1;
    y = cumsum(Ybinned)./sum(Ybinned);
    plot(x,y, '-k');
    line([0 1], [0 1], 'Color', 'b', 'LineStyle', '--');
    lt_subplot(3,2,6); hold on;
    title('all cases, percentiles rel own shuffle');
    xlabel('log10(p)');
    [Ybinned] = lt_plot_histogram(log10(pvalmat+(1/size(yshuff,2))));
    
    
end

% =========== QQ PLOT
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% ================== [for each plot against uniform distribution]
for j=1:length(ALLBRANCH)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    pd = makedist('Uniform');
    qqplot(PvalsAll{j}, pd);
    title(['qq vs. uniform, ' ALLBRANCH(j).bregion]);
    line([0 1], [0 1], 'Color', 'k');
end


% ############################# QQ DIFERENCES BETWEEN BRAIN REGIONS
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
[~, xx, yy] = lt_plot_qqplot(PvalsAll{1}, PvalsAll{2});
xlabel(ALLBRANCH(1).bregion);
ylabel(ALLBRANCH(2).bregion);
title('qq comparing 2 datasets');
line([0 1], [0 1], 'Color', 'k');

% ====== COLLECT THE DIFFERENCE FROM UNITY LINE
qqdiffDAT = mean(yy - xx');


% ######################### GET SHUFFLE DISTRIBUTION OF QQ DIFFERENCES
% for all 1000 shuffles, recalculate p values for a shuffle rendition
% (relative to the rest of the distribution, also throwing in the data
% vector...

% ================= COMBINE DATA AND SHUFFLE FOR EACH BRAIN REGION
DatMatAll = {};
for j=1:length(ALLBRANCH)
    
    ydat = ALLBRANCH(j).DAT.(datfield);
    yshuff = [ALLBRANCH(j).DAT.AllBranch_SHUFFSTRUCT.cycle.(shufffield)];
    
    % -------------------------------------
    indsgood = ALLBRANCH(j).IndsGood;
    % -------------------------------------
    
    DatMatAll{j} = [ydat(indsgood,:) yshuff(indsgood,:)];
end

qqDiffShuffAll = [];
for j=2:size(DatMatAll{1},2) % start from 2 since 1 is actual dat.
    disp(num2str(j));
    %     shuffmat = yshuff;
    %     datmat = repmat(ydat, 1, size(shuffmat,2));
    %
    %     % ================= get pvalues (empirical
    %     pvalmat = (1+sum(shuffmat>=datmat, 2))./size(shuffmat,2);
    %     PvalsAll = [PvalsAll pvalmat];
    %
    
    % ======================== COLLECT PVALS FOR EACH DATASET
    PvalsAll_shuff = {};
    for jj=1:length(ALLBRANCH)
        
        ytmp = DatMatAll{jj}(:,j);
        yshufftmp = DatMatAll{jj}(:,[1:j-1 j+1:end]);
        
        pvalshuff = (1+sum(yshufftmp>=ytmp, 2))./(size(yshufftmp,2)+1);
        PvalsAll_shuff{jj} = pvalshuff;
    end
    
    % ======================= do qq plot to get qq diff
    [~, xx, yy] = lt_plot_qqplot(PvalsAll_shuff{1}, PvalsAll_shuff{2}, [], 0);
    
    qqdiff = mean(yy - xx');
    qqDiffShuffAll = [qqDiffShuffAll; qqdiff];
end

pval_qqdiff = (1 + sum(abs(qqDiffShuffAll)>=abs(qqdiffDAT)))./(length(qqDiffShuffAll)+1);
lt_plot_pvalue(pval_qqdiff, 'qqdiff, 2 datasets', 1);

%% ============= COMPARE EFFECT SIZES

datfield = 'AllBranch_SlopeCoeff_MedAbs';
shufffield = 'SlopeCoeff_MedAbs';
% datfield = 'AllBranch_SlopeDiff';
% shufffield = 'SlopeDiff';


lt_figure; hold on;

% ========== plot each case on its own
lt_subplot(4,2,1); hold on;
xlabel([ALLBRANCH(1).bregion '-' ALLBRANCH(2).bregion]);
ylabel('median(abs(effect))');

Y = {};
Yshuff = {};
for i=1:length(ALLBRANCH)
    ydat = ALLBRANCH(i).DAT.(datfield);
    yshuff = [ALLBRANCH(i).DAT.AllBranch_SHUFFSTRUCT.cycle.(shufffield)];
    
    indsgood = ALLBRANCH(i).IndsGood;
    ydat = ydat(indsgood,:);
    yshuff = yshuff(indsgood,:);
    
    
    % ----- take median across all cases
    if strcmp(datfield, 'AllBranch_SlopeDiff')
        % --- then get fraction
        ydat = sum(ydat)/length(ydat);
        yshuff = sum(yshuff,1)/size(yshuff,1);
    else
        ydat = median(ydat);
        yshuff = median(yshuff,1);
    end
    
    % ---- overlay shuffle distribution
    %    lt_plot_MultDist({yshuff}, i, 0, [0.7 0.7 0.7])
    yshuffCI = prctile(yshuff, [2.5 97.5]);
    patch([i-0.3 i+0.3 i+0.3 i-0.3], [yshuffCI(1) yshuffCI(1) yshuffCI(2) yshuffCI(2)], [0.7 0.7 0.7]);
    set(gca, 'Color', 'none');
    line([i-0.3 i+0.3], [median(yshuff) median(yshuff)], 'Color', [0.4 0.4 0.4]);
    
    % ---- plot dat
    lt_plot(i, ydat);
    
    % ============= COLLECT
    Y{i} = ydat;
    Yshuff{i} = yshuff;
end

% ================ PLOT DIFFERENCE
lt_subplot(4,2,2); hold on;
xlabel('diff in median abs effect');

% --- shuffle (each iteraction calcualte medians and differnce
ydiff_shuff = Yshuff{2} - Yshuff{1};
lt_plot_histogram(ydiff_shuff);

% --- data
ydiff = Y{2} - Y{1};
line([ydiff ydiff], ylim, 'Color', 'r');

% --- pval
p = (1+sum(abs(ydiff_shuff) >= abs(ydiff)))./(1+length(ydiff_shuff));
lt_plot_pvalue(p, 'vs shuff', 1);


%%
% %% ============= [DAT VS. SHUFFLE] COMPARE DISTRIBUTIONS
% lt_figure; hold on;
% shuffcycle = 1; % which one ot plot;
% onlysigdecode = 0; % then filters based on decode [0 1 or 2]
% NctxtToPlot = [1:10]; % can't be empty
%
% if onlysigdecode==1
%     inds = AllBranch_DecodeP<0.05 & ismember(AllBranch_Nctxt, NctxtToPlot);
% elseif onlysigdecode==0
%     inds = logical(ones(size(AllBranch_DecodeP))) & ismember(AllBranch_Nctxt, NctxtToPlot);
% elseif onlysigdecode==2
%     inds = AllBranch_DecodeP>=0.05 & ismember(AllBranch_Nctxt, NctxtToPlot);
% end
%
% % ---- NCASES, slope interaction
% lt_subplot(4,2,1); hold on;
% title('ncases slope interaction');
% y = sum(AllBranch_SlopeDiff(inds));
% tmp = [AllBranch_SHUFFSTRUCT.cycle.SlopeDiff];
% yshuff = sum(tmp(inds,:),1);
%
% lt_plot_histogram(yshuff);
% line([y, y], ylim, 'Color', 'r');
% lt_plot_annotation(1, ['N=' num2str(sum(inds))], 'b');
% % pval
% p = (1+sum(yshuff>=y))./length(yshuff);
% lt_plot_pvalue(p, 'vs. shuff', 1);
%
%
% % ---- [SINGLE SHUFF] median absolute(slope-interaction effect)
% lt_subplot(4,2,2); hold on;
% title('[SINGLE SHUFF]median(across classes) abs(slope-context interact)');
% xlabel('bk=shuff');
% y = AllBranch_SlopeCoeff_MedAbs(inds);
% yshuff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeCoeff_MedAbs(inds);
% [~, xcenters] = lt_plot_histogram(y, '', 1, 1, '', 1, 'r');
% lt_plot_histogram(yshuff, xcenters, 1, 1, '', 1, 'k');
%
% p = signrank(y, yshuff);
% lt_plot_pvalue(p, 'srank', 1);
%
% % ---- [MULT SHUFF, each with median] median absolute(slope-interaction effect)
% lt_subplot(4,2,3); hold on;
% title('[MULT SHUFF]median(across classes) abs(slope-context interact)');
%
% y = median(AllBranch_SlopeCoeff_MedAbs(inds));
% tmp = [AllBranch_SHUFFSTRUCT.cycle.SlopeCoeff_MedAbs];
% yshuff = median(tmp(inds,:), 1);
% lt_plot_histogram(yshuff, '', 1, 1, '', 1, 'k');
% line([y y], ylim, 'Color', 'r');
% % -pvalue
% p = (1+sum(yshuff>=y))./length(yshuff);
% lt_plot_pvalue(p, 'vs. shuff', 1);
% lt_plot_zeroline_vert;
%
%
% % ---- [SINGLE SHUFF] overall slopes
% lt_subplot(4,2,4); hold on;
% title('[SINGLE SHUFF] median(across classes) abs(overall slope)');
% xlabel('bk=shuff');
% y = AllBranch_SlopeOverallCoeff_MedAbs(inds);
% yshuff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeOverallCoeff_MedAbs(inds);
% [~, xcenters] = lt_plot_histogram(y, '', 1, 1, '', 1, 'r');
% lt_plot_histogram(yshuff, xcenters, 1, 1, '', 1, 'k');
%
% p = signrank(y, yshuff);
% lt_plot_pvalue(p, 'srank', 1);
%
%
% % ---- [SINGLE SHUFF] intercept x class interaction
% lt_subplot(4,2,5); hold on;
% title('[SINGLE SHUFF] median(across classes) abs(intercept-context interaction)');
% xlabel('bk=shuff');
% y = AllBranch_IntCoeff_MedAbs(inds);
% yshuff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).IntCoeff_MedAbs(inds);
% [~, xcenters] = lt_plot_histogram(y, '', 1, 1, '', 1, 'r');
% lt_plot_histogram(yshuff, xcenters, 1, 1, '', 1, 'k');
%
% p = signrank(y, yshuff);
% lt_plot_pvalue(p, 'srank', 1);
%
%
%
% %% ============= [DAT VS. SHUFFLE] PLOT
% %% [PLOT] OVERLAY DATA WITH ONE INSTANCE OF SHUFFLE
% lt_figure; hold on;
% count = 1;
% decodelist = [1 2]; % iterate over diff contingencies on decode performance
% plotfrac = 1;
% shuffcycle = 1;
% for useshuff = [0 1]
%     for onlydecode = decodelist
%         % onlydecode = 1;
%         % 0: don't care
%         % 1: only cases with significant decode
%         % 2: only cases with NON significant decode
%         lt_subplot(2,2,count); hold on;
%         if plotfrac==1
%             ylabel('fraction syls');
%         else
%             ylabel('n syls');
%         end
%
%         xlabel('int_only - slope_only - both - neither(slope_main) - neither(no slope_main)');
%         for i=1:numbirds
%
%             if onlydecode==0
%                 inds = AllBranch_birdnum==i;
%                 title('all data (dont care about decode)');
%             elseif onlydecode==1
%                 inds =  AllBranch_birdnum==i & AllBranch_DecodeP<0.05;
%                 title('only if significant decode');
%             elseif onlydecode==2
%                 inds =  AllBranch_birdnum==i & AllBranch_DecodeP>=0.05;
%                 title('only if NONsig decode');
%             end
%             line([i+0.5 i+0.5], ylim);
%
%             if ~any(inds)
%                 continue
%             end
%
%             if useshuff==0
%                 % then use dat
%                 intdiff = AllBranch_IntDiff(inds);
%                 slopediff = AllBranch_SlopeDiff(inds);
%                 slopeoverall = AllBranch_SlopeOverall(inds);
%             else
%                 intdiff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).IntDiff(inds);
%                 slopediff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeDiff(inds);
%                 slopeoverall = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeOverall(inds);
%             end
%
%             Y = [];
%             % --- only int
%             y = sum(intdiff==1 & slopediff==0);
%             Y = [Y y];
%
%             % --- only slope
%             y = sum(intdiff==0 & slopediff==1);
%             Y = [Y y];
%
%             % --- both
%             y = sum(intdiff==1 & slopediff==1);
%             Y = [Y y];
%
%             % --- none (with overall slope effect)
%             y = sum(intdiff==0 & slopediff==0 & slopeoverall==1);
%             Y = [Y y];
%
%             % --- none (WITHOUT overall slope effect)
%             y = sum(intdiff==0 & slopediff==0 & slopeoverall==0);
%             Y = [Y y];
%
%             % --- nromalize to get fraction
%             if plotfrac==1
%                 Y = Y./length(intdiff);
%             end
%
%             % --- get x locations
%             X = i+[-0.25 -0.1 0.05 0.2 0.35];
%
%             % ---------------- SIG EFFECTS
%             lt_plot_bar(X(1:3), Y(1:3), {'BarWidth', 0.6, 'Color', 'r'});
%
%             % ---------------- NOT SIG (with slope efefct
%             lt_plot_bar(X, [nan nan nan Y(4) nan], {'BarWidth', 0.6});
%
%             % --------------- NOT SIG (no slope effect)
%             lt_plot_bar(X, [nan nan nan nan Y(5)], {'BarWidth', 0.6, 'Color', 'none'});
%
%
%         end
%
%         % ###################### across entire dataset
%         if onlydecode==0
%             inds = ones(size(AllBranch_birdnum));
%         elseif onlydecode==1
%             inds =  AllBranch_DecodeP<0.05;
%         elseif onlydecode==2
%             inds =  AllBranch_DecodeP>=0.05;
%         end
%         line([i+0.5 i+0.5], ylim);
%
%             if useshuff==0
%                 % then use dat
%                 intdiff = AllBranch_IntDiff(inds);
%                 slopediff = AllBranch_SlopeDiff(inds);
%                 slopeoverall = AllBranch_SlopeOverall(inds);
%             else
%                 intdiff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).IntDiff(inds);
%                 slopediff = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeDiff(inds);
%                 slopeoverall = AllBranch_SHUFFSTRUCT.cycle(shuffcycle).SlopeOverall(inds);
%             end
%
%         Y = [];
%         % --- only int
%         y = sum(intdiff==1 & slopediff==0);
%         Y = [Y y];
%
%         % --- only slope
%         y = sum(intdiff==0 & slopediff==1);
%         Y = [Y y];
%
%         % --- both
%         y = sum(intdiff==1 & slopediff==1);
%         Y = [Y y];
%
%         % --- none (with overall slope effect)
%         y = sum(intdiff==0 & slopediff==0 & slopeoverall==1);
%         Y = [Y y];
%
%         % --- none (WITHOUT overall slope effect)
%         y = sum(intdiff==0 & slopediff==0 & slopeoverall==0);
%         Y = [Y y];
%
%         % --- nromalize to get fraction
%         if plotfrac==1
%             Y = Y./length(intdiff);
%         end
%
%         % --- get x locations
%         X = numbirds+1+[-0.25 -0.1 0.05 0.2 0.35];
%
%         % ---------------- SIG EFFECTS
%         lt_plot_bar(X(1:3), Y(1:3), {'BarWidth', 0.6, 'Color', 'r'});
%
%         % ---------------- NOT SIG (with slope efefct
%         lt_plot_bar(X, [nan nan nan Y(4) nan], {'BarWidth', 0.6});
%
%         % --------------- NOT SIG (no slope effect)
%         lt_plot_bar(X, [nan nan nan nan Y(5)], {'BarWidth', 0.6, 'Color', 'none'});
%
%
%         % ===================== formatting
%         set(gca, 'XTick', 1:numbirds+1, 'XTickLabel', [{SummaryStruct.birds.birdname} 'ALL']);
%         rotateXLabels(gca, 90);
%         count = count+1;
%         if useshuff==1
%             lt_plot_annotation(1, 'SHUFF', 'b');
%         end
%     end
% end
%
% %% = PLOT "RASTER" OF ALL CASES, INDICATING MAGNITUDE AND SIGNIFICANCE
% effecttosortby = 2;
%
% plotshuffle = 1; % have to have done shuffle already
% sc = 1; % choose which one.
%
% % ================
% lt_figure; hold on;
% ylabel('case #');
% xlabel('Int*ctxt -- Slope*ctxt -- SlopeOverall');
%
% % --- ciollect data
% if plotshuffle==0
% Y = [AllBranch_IntCoeff_MedAbs AllBranch_SlopeCoeff_MedAbs AllBranch_SlopeOverallCoeff_MedAbs];
% Ysig = [AllBranch_IntDiff AllBranch_SlopeDiff AllBranch_SlopeOverall];
% Ydecode = AllBranch_DecodeP<0.05;
% title('[DAT]median abs(cofficeint)');
% elseif plotshuffle==1
% Y = [AllBranch_SHUFFSTRUCT.cycle(sc).IntCoeff_MedAbs AllBranch_SHUFFSTRUCT.cycle(sc).SlopeCoeff_MedAbs ...
%     AllBranch_SHUFFSTRUCT.cycle(sc).SlopeOverallCoeff_MedAbs];
% Ysig = [AllBranch_SHUFFSTRUCT.cycle(sc).IntDiff AllBranch_SHUFFSTRUCT.cycle(sc).SlopeDiff ...
%     AllBranch_SHUFFSTRUCT.cycle(sc).SlopeOverall];
% Ydecode = AllBranch_DecodeP<0.05;
% title('[SHUFF]median abs(cofficeint)');
% end
% % ---- sort data, increasing effect *(choose effect)
% [~, inds] = sort(Y(:,effecttosortby));
% Y = Y(inds,:);
% Ysig = Ysig(inds, :);
% Ydecode = Ydecode(inds,:);
%
%
% % ---- plot
% x = 1:size(Y,2);
% imagesc(Y, [0 2]);
% colorbar;
% colormap('gray');
%
% % ---- mark those that are significnat
% for j=1:size(Ysig,2)
%    plot(j-0.3, find(Ysig(:,j)), 'or');
% end
%
% % ---- mark those with significnat decode (regression)
% plot(0.3, find(Ydecode), 'ob');
%
% % lt_figure; hold on;
% % for j=1:size(Y,2)
% %     lt_plot_stem3(j*ones(size(Y,1),1), 1:size(Y,1)', Y(:, j), 'k', 1);
% % end

