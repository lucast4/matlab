function lt_neural_v2_BRANCH_DecodePlot(DecodeStruct)
anum = 1;
plottext =0;

%% ==================== general params

DatAll = DecodeStruct.analynum(anum).dat;
DATSTRUCT_BYBRANCH = DecodeStruct.analynum(anum).datbybranch;
SummaryStruct = DecodeStruct.analynum(anum).SummaryStruct;

Numbirds = length(DecodeStruct.analynum(anum).datbybranch.bird);
Numbranch = max(DatAll.AllBranchNum);

BrainRegions = unique(DecodeStruct.analynum(anum).dat.AllBrainRegion);
plotcols_brainreg = lt_make_plot_colors(length(BrainRegions), 0,0);

%% ================== [PLOT] each bird, all neur x branch
ColorList = {'sig', 'branchnum'};
for howtocolor = ColorList
    % howtocolor = 'branchnum';
    % sig = by p<0.05 threshold
    % branchnum = by branch id
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for bregion = BrainRegions'
        pcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        for i=1:Numbirds
            birdname = SummaryStruct.birds(i).birdname;
            
            % ================ iterate over number of contexts
            indtmp = DatAll.AllBirdNum==i & strcmp(DatAll.AllBrainRegion, bregion);
            if ~any(indtmp)
                continue
            end
            
            % ------------- what is max number of context across all branches?
            X = unique(DatAll.AllNumCtxts(indtmp)); % collect context values
            
            % ============ alternative output (for plotspread)
            inds = DatAll.AllBirdNum==i & strcmp(DatAll.AllBrainRegion, bregion);
            Ydat_vec = DatAll.AllDecode(inds)';
            Yneg_vec = DatAll.AllDecode_neg_mean(inds);
            X_vec = DatAll.AllNumCtxts(inds);
            Pval_vec = DatAll.AllPdat(inds);
            Ybranchnum = DatAll.AllBranchNum(inds);
            
            functmp1 = @(x) prctile(x, [2.5]);
            functmp2 = @(x) prctile(x, [97.5]);
            [tmp1, tmp2, tmp3] = grpstats(Yneg_vec, X_vec, {'mean', functmp1, functmp2});
            Yshuff_means = tmp1;
            Yshuff_CI = [tmp2 tmp3];
            
            
            
            % =================== PLOT FOR THIS BIRD
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' bregion]);
            xlabel('num ctxts in branch');
            ylabel('decode (F1) [col=p<0.05], gray=empirical neg');
            xlim([min(DatAll.AllNumCtxts)-0.5 max(DatAll.AllNumCtxts)+0.5]);
            
            % ------ overlay theoretical neg
            if (0)
                %         yneg_theor = 1./X;
                for x = X'
                    line([x-0.5 x+0.5], [1/x 1/x], 'Color', [0.2 0.3 0.4], ...
                        'LineStyle', '--')
                end
            end
            
            % ------ overlay empirical neg
            for j=1:length(X)
                x = [X(j)-0.5 X(j)+0.5];
                y = Yshuff_means(j);
                yCI = Yshuff_CI(j,:);
                
                patch([x(1) x(2) x(2) x(1)], [yCI(1) yCI(1) yCI(2) yCI(2)], [0.8 0.8 0.8], ...
                    'EdgeColor', 'none');
                line(x, [y y], 'Color', [0.6 0.6 0.6]);
                
                %            shadedErrorBar(x, [y y], flipud([yCI' yCI']-y), {'Color', [0.7 0.7 0.7]},1);
            end
            
            % ----- plot dat
            %         lt_plot_MultDist(Ydat, X, 0, pcol);
            % -- not significant
            if strcmp(howtocolor, 'sig')
                plotSpread(Ydat_vec, 'distributionIdx', X_vec, 'categoryIdx', Pval_vec<0.05, ...
                    'categoryColors', {'k', pcol});
            elseif strcmp(howtocolor, 'branchnum')
                plotcolstmp = lt_make_plot_colors(length(unique(Ybranchnum)), 0,0);
                plotSpread(Ydat_vec, 'distributionIdx', X_vec, 'categoryIdx', Ybranchnum, ...
                    'categoryColors', plotcolstmp);
                
            end
            %         lt_plot_MultDist(Ydat, X, 0, pcol, 1);
            
            
            % --------- overlay mean and SE
            [ymean, ysem] = grpstats(Ydat_vec, X_vec, {'mean', @lt_sem});
            %         ymean = cellfun(@mean, Ydat);
            %         ysem = cellfun(@lt_sem, Ydat);
            lt_plot(X+0.3, ymean, {'Errors', ysem, 'Color', pcol});
            
            ylim([0 1]);
        end
        
        %% ======================== COMBINE DATA ACROSS ALL BIRDS
        % ============ alternative output (for plotspread)
        inds = strcmp(DatAll.AllBrainRegion, bregion);
        
        X = unique(DatAll.AllNumCtxts(inds)); % collect context values
        Ydat_vec = DatAll.AllDecode(inds)';
        Yneg_vec = DatAll.AllDecode_neg_mean(inds);
        X_vec = DatAll.AllNumCtxts(inds);
        Pval_vec = DatAll.AllPdat(inds);
        Ybranchnum = DatAll.AllBranchNum(inds);
        
        functmp1 = @(x) prctile(x, [2.5]);
        functmp2 = @(x) prctile(x, [97.5]);
        [tmp1, tmp2, tmp3] = grpstats(Yneg_vec, X_vec, {'mean', functmp1, functmp2});
        Yshuff_means = tmp1;
        Yshuff_CI = [tmp2 tmp3];
        
                
        % =================== PLOT FOR THIS BIRD
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['ALL BIRDS']);
        xlabel('num ctxts in branch');
        ylabel('decode (F1) [col=p<0.05], gray=empirical neg');
        xlim([min(DatAll.AllNumCtxts)-0.5 max(DatAll.AllNumCtxts)+0.5]);
        
        % ------ overlay theoretical neg
        if (0)
            %         yneg_theor = 1./X;
            for x = X'
                line([x-0.5 x+0.5], [1/x 1/x], 'Color', [0.2 0.3 0.4], ...
                    'LineStyle', '--')
            end
        end
        
        % ------ overlay empirical neg
        for j=1:length(X)
            x = [X(j)-0.5 X(j)+0.5];
            y = Yshuff_means(j);
            yCI = Yshuff_CI(j,:);
            
            patch([x(1) x(2) x(2) x(1)], [yCI(1) yCI(1) yCI(2) yCI(2)], [0.8 0.8 0.8], ...
                'EdgeColor', 'none');
            line(x, [y y], 'Color', [0.6 0.6 0.6]);
            
            %            shadedErrorBar(x, [y y], flipud([yCI' yCI']-y), {'Color', [0.7 0.7 0.7]},1);
        end
        
        % ----- plot dat
        %         lt_plot_MultDist(Ydat, X, 0, pcol);
        % -- not significant
        if strcmp(howtocolor, 'sig')
            plotSpread(Ydat_vec, 'distributionIdx', X_vec, 'categoryIdx', Pval_vec<0.05, ...
                'categoryColors', {'k', pcol});
        elseif strcmp(howtocolor, 'branchnum')
            plotcolstmp = lt_make_plot_colors(length(unique(Ybranchnum)), 0,0);
            plotSpread(Ydat_vec, 'distributionIdx', X_vec, 'categoryIdx', Ybranchnum, ...
                'categoryColors', plotcolstmp);
            
        end
        %         lt_plot_MultDist(Ydat, X, 0, pcol, 1);
        
        
        % --------- overlay mean and SE
        [ymean, ysem] = grpstats(Ydat_vec, X_vec, {'mean', @lt_sem});
        %         ymean = cellfun(@mean, Ydat);
        %         ysem = cellfun(@lt_sem, Ydat);
        lt_plot(X+0.3, ymean, {'Errors', ysem, 'Color', pcol});
        
        ylim([0 1]);
        
    end
end


%% ============== [SEPARATE BY BRANCH ID];
% =============== FOR EACH BIRD SEPARATE BY BRANCH
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



for bregion = BrainRegions'
    pcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
    
    Nall = [];
    Nsig = [];
    BirdnumAll = [];
    for i=1:Numbirds
        birdname = SummaryStruct.birds(i).birdname;
        
        % ================ iterate over number of contexts
        inds = DatAll.AllBirdNum==i & strcmp(DatAll.AllBrainRegion, bregion);
        if ~any(inds)
            continue
        end
        
        % ------------- what is max number of context across all branches?
        %             X = unique(DatAll.AllNumCtxts(inds)); % collect context values
        
        % ============ alternative output (for plotspread)
        Ydat_vec = DatAll.AllDecode(inds)';
        Yneg_vec = DatAll.AllDecode_neg_mean(inds);
        Nctxt = DatAll.AllNumCtxts(inds);
        Pval_vec = DatAll.AllPdat(inds);
        Ybranchnum = DatAll.AllBranchNum(inds);
        
        %             functmp1 = @(x) prctile(x, [2.5]);
        %             functmp2 = @(x) prctile(x, [97.5]);
        %             [tmp1, tmp2, tmp3] = grpstats(Yneg_vec, X_vec, {'mean', functmp1, functmp2});
        %             Yshuff_means = tmp1;
        %             Yshuff_CI = [tmp2 tmp3];
        %
        
        % ======== FOR EACH BRANCH, DETERMINE THE NUMBER OF SIGNIFICANT
        % CASES
        for j=1:max(Ybranchnum)
            
            indtmp = Ybranchnum==j;
            if ~any(indtmp)
                continue
            end
            
            Nall = [Nall; sum(Pval_vec(indtmp)<=1)];
            Nsig = [Nsig; sum(Pval_vec(indtmp)<=0.05)];
            BirdnumAll = [BirdnumAll; i];
        end
        
        
        % =================== PLOT FOR THIS BIRD
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title([birdname '-' bregion]);
        xlabel('nsites total (per syl)');
        ylabel('nsites (sig)');
        
        indtmp = BirdnumAll==i;
        plot(Nall(indtmp)+0.2*rand(sum(indtmp),1)-0.2, Nsig(indtmp)+0.2*rand(sum(indtmp),1)-0.2, 'o', 'Color',pcol);
        lt_plot_makesquare_plot45line(gca, 'k', [-1]);
    end
    
    % ================= all birds for this brain region
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['ALL BIRDS-' bregion]);
        xlabel('nsites total (per syl)');
        ylabel('nsites (sig)');
        
        plot(Nall+0.2*rand(size(Nall))-0.2, Nsig+0.2*rand(size(Nsig))-0.2, 'o', 'Color',pcol);
        lt_plot_makesquare_plot45line(gca, 'k', [-1]);
        
        % ============= % syls with at least one significant neuron,
        % stratified by number of neurons
        X = unique(Nall);
        Y = [];
        for x=X'
           
            numsylssig = sum(Nsig(Nall == x)>0)/sum(Nall == x);
            
               Y = [Y; numsylssig];         
        end
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['ALL BIRDS-' bregion]);
        xlabel('num neurons for this syl');
        ylabel('fraction syls w/a/l one neuron sig');
        plot(X, Y, 'ok');
        
        
        %% =================== for each value of NumNeurons, downsample any
        % sylalbles that have more than that number of neurons. then will
        % have monotonically decreasing size of dataset
        [X, FracSig, Nsamp] = fn_Nsigdecode_dsamp(Nall, Nsig);
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['ALL BIRDS-' bregion{1} '[w/downsamp]']);
        xlabel('num neurons for this syl');
        ylabel('fraction syls w/a/l one neuron sig');
        plot(X, FracSig, 'ok');
        for j=1:length(X)
            lt_plot_text(X(j), FracSig(j), ['N=' num2str(Nsamp(j))], 'r', 8)
        end
        
        
        % ============================= SAME, BUT SEPARATELY FOR EACH BIRD
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(['EACH BIRD-' bregion{1} '[w/downsamp]']);
        xlabel('num neurons for this syl');
        ylabel('fraction syls w/a/l one neuron sig');
        
        % --------------------- go thru each bird
        pcolsbirds = lt_make_plot_colors(max(BirdnumAll), 0, 0);
        Bnames = {};
        for j=1:max(BirdnumAll)
            indtmp = BirdnumAll==j;
            if ~any(indtmp)
                continue
            end
            
            % ================= do multiple cycles of donwsampling and then
            % take average
            ncycles = 10;
            FracSig = [];
            for n=1:ncycles
            [X, FracSigTmp, Nsamp] = fn_Nsigdecode_dsamp(Nall(indtmp), Nsig(indtmp));
            FracSig = [FracSig FracSigTmp];
            end
%             shadedErrorBar(X, mean(FracSig,2), std(FracSig'), {'Color', pcolsbirds{j}}, 1);
            plot(X, mean(FracSig,2), '-', 'Color', pcolsbirds{j});
            
            Bnames = [Bnames SummaryStruct.birds(j).birdname];
        end
        legend(gca, Bnames);
        
        % ---- overlay all birds combined
        if (0)
        [X, FracSig, Nsamp] = fn_Nsigdecode_dsamp(Nall, Nsig);
        plot(X, FracSig, '-k', 'LineWidth', 2);
        end
        
        ylim([0 1]);
        
        
        % --- fit negative exponential model
        
end


%% ============== [REGRESSION MODEL] decode
% NOTE: THIS IS LMAN ONLY
LMANonly =1;

% ===================== 1) negative shuffle
subtractshuffmean = 0; % for each case, subtract it's own neg control
UseDecodeShuff = 1; % then uses shuf instaed of actual decode...
lme_NEG = lt_neural_v2_CTXT_BRANCH_DecodePlot_Regr(DatAll, subtractshuffmean, ...
    UseDecodeShuff, LMANonly);
disp(lme_NEG);

% ==================== 2) dat - negative.
subtractshuffmean = 1; % for each case, subtract it's own neg control
UseDecodeShuff = 0; % then uses shuf instaed of actual decode...
lme_POSminusNEG = lt_neural_v2_CTXT_BRANCH_DecodePlot_Regr(DatAll, subtractshuffmean, ...
    UseDecodeShuff,LMANonly);
disp(lme_POSminusNEG);


%% ============== MAKE PLOT OF COEFFICIENTS
lt_figure; hold on;

% ************************************************* REGRESSION COEFFICIENTS
lt_subplot(1,2,1); hold on;
title('pop estimate (lme)');
ylabel('dat minus shuff (overlayed on shuff)');

% #################################### 1) PLOT NEGATIVE
lme = lme_NEG;

% --- fixed effects
fixednames = lme.CoefficientNames;
fixedeff = lme.fixedEffects;

% --- get coeff (i.e. add to intercept)
coeff_intercept = fixedeff(strcmp(fixednames, '(Intercept)'));
assert(find(strcmp(fixednames, '(Intercept)')) ==1, 'assuming non-first index are all fixed effects');

fixedeff(2:end) = fixedeff(2:end) + coeff_intercept;

fixedeffCI = lme.coefCI;
fixedeffCI(2:end, :) = fixedeffCI(2:end,:) + coeff_intercept;



% ----------------- plot a patch
for x = 1:length(fixednames)
    xthis = [x-0.3 x+0.3];
    ythis = fixedeff(x);
    yCI = fixedeffCI(x,:);
    patch([xthis(1) xthis(2) xthis(2) xthis(1)], [yCI(1) yCI(1) yCI(2) yCI(2)], [0.8 0.8 0.8], ...
        'EdgeColor', 'none');
    line(xthis, [ythis ythis], 'Color', [0.3 0.3 0.3], 'LineWidth', 2);
end

% ------------- save fixed ffects
fixedCoeff_NEG = fixedeff;



% #################################### 2) PLOT (POSITIVE - NEGATIVE) + POSITIVE)
lme = lme_POSminusNEG;

% --- fixed effects
fixednames = lme.CoefficientNames;
fixedeff = lme.fixedEffects;

% --- get coeff (i.e. add to intercept)
coeff_intercept = fixedeff(strcmp(fixednames, '(Intercept)'));
assert(find(strcmp(fixednames, '(Intercept)')) ==1, 'assuming non-first index are all fixed effects');

fixedeff(2:end) = fixedeff(2:end) + coeff_intercept;

fixedeffCI = lme.coefCI;
fixedeffCI(2:end, :) = fixedeffCI(2:end,:) + coeff_intercept;


% --------------- ADD TO NEGATIVE SHUFFLE
fixedCoeff_POS = fixedeff + fixedCoeff_NEG;
fixedCoeffCI_POS = fixedeffCI + fixedCoeff_NEG;


% ----------------- PLOT DATA
x = 1:length(fixednames);
errorbar(x, fixedCoeff_POS,  fixedCoeff_POS-fixedCoeffCI_POS(:,1), ...
    -fixedCoeff_POS+fixedCoeffCI_POS(:,2), 'LineStyle', 'none', ...
    'Color', 'k');
lt_plot(x, fixedCoeff_POS, {'Color', 'k'});

% --------- plot p values
sigthis = fixedeffCI(:,1)>0;
for x=1:length(sigthis)
    if sigthis(x) ==1
        lt_plot_text(x-0.1, fixedCoeff_POS(x)+0.1, '*', 'r', 10)
    end
    
end

ylim([-0.1 1]);
xlim([0 x(end)+1]);
lt_plot_zeroline;
set(gca, 'XTick', x, 'XTickLabel', fixednames);
rotateXLabels(gca, 45);




% ************************************************* ACTUAL DATA
lt_subplot(1,2,2); hold on;
title('data');
ylabel('dat (solid = sign, vs. own shuff (signrank)');

% ------- for each ctxt class and each bird get mean of dat
CtxtsToPlots = [2 3 4];
assert(LMANonly==1, 'this is coded for lman only..');
pcolors = lt_make_plot_colors(Numbirds, 0,0);

for i =1:Numbirds
    
    Yall = []; % by context
    Ystd = [];
    Pvals = [];
    xthis = [];
    for cc = 1:length(CtxtsToPlots)
        ctxtnum_this = CtxtsToPlots(cc);
        
        inds = DatAll.AllBirdNum==i & DatAll.AllNumCtxts==ctxtnum_this ...
            & strcmp(DatAll.AllBrainRegion, 'LMAN');
        
        if ~any(inds)
            continue
        end
        
        ydecode = DatAll.AllDecode(inds)';
        ymean = mean(ydecode);
        ystd = std(ydecode);
        
        Yall = [Yall ymean];
        Ystd = [Ystd ystd];
        xthis = [xthis cc];
        
        % ---- this bird and context num sig diff from its own shuffle>?
        ydecode_neg = DatAll.AllDecode_neg_mean(inds);
        
        p = signrank(ydecode, ydecode_neg, 'tail', 'right');
        Pvals = [Pvals p];
    end
    
    if isempty(Yall)
        continue
    end
    
    
    %     xthis = 1:length(CtxtsToPlots);
    xthis = xthis + 0.6*rand-0.3;
    %     pcol = [rand rand rand];
    pcol = pcolors{i};
    
    % -- nonsig
    indtmp = Pvals>=0.05;
    if ~isempty(indtmp)
        lt_plot(xthis(indtmp), Yall(indtmp), {'Errors', Ystd(indtmp), 'Color', pcol, ...
            'MarkerFaceColor', 'none'});
    end
    
    % -- sig
    indtmp = Pvals<0.05;
    if ~isempty(indtmp)
        lt_plot(xthis(indtmp), Yall(indtmp), {'Errors', Ystd(indtmp), 'Color', pcol});
    end
    
    % --- connect all
    plot(xthis, Yall, '-', 'Color', pcol);
end

ylim([-0.1 1]);
xlim([0 x(end)+1]);
lt_plot_zeroline;
set(gca, 'XTick', 1:length(CtxtsToPlots), 'XTickLabel', CtxtsToPlots);
rotateXLabels(gca, 45);


%% ==================== combine all contexts, signifncat?

subtractshuffmean = 1; % for each case, subtract it's own neg control
UseDecodeShuff = 0; % then uses shuf instaed of actual decode...
combineCtxts = 1;
lme = lt_neural_v2_CTXT_BRANCH_DecodePlot_Regr(DatAll, subtractshuffmean, ...
    UseDecodeShuff,LMANonly, combineCtxts);
disp(lme);

%% REGRESSION STEPS:
% 1) get model coefficients

%%

lme = lme2;
%% ================= [REGRESSION] PLOT COEFFICIENTS

fixednames = lme.CoefficientNames;
fixedeff = lme.fixedEffects;

coeff_intercept = fixedeff(strcmp(fixednames, '(Intercept)'));

fixedeff(2:end) = fixedeff(2:end) + coeff_intercept;

fixedeffCI = lme.coefCI;
fixedeffCI(2:end, :) = fixedeffCI(2:end,:) + coeff_intercept;

lt_figure; hold on;
x = 1:length(fixednames);
errorbar(x, fixedeff,  fixedeff-fixedeffCI(:,1), -fixedeff+fixedeffCI(:,2), 'LineStyle', 'none');
xlim([0 x(end)+1]);
lt_plot_zeroline;
set(gca, 'XTick', x, 'XTickLabel', fixednames);
rotateXLabels(gca, 45);

% --------- plot p values
YLIM = ylim;
pvals = lme.Coefficients.pValue;
for x=1:length(pvals)
    if pvals(x)<0.1
        lt_plot_text(x-0.3, YLIM(2), ['p=' num2str(pvals(x))], 'r', 10)
    end
end

ylim([-1 1]);


%% ================== [PLOT] effect sizes (i.e. decode F1) for all

% ------------------
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];



% =============== PLOT
for bregion = BrainRegions'
    plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
    for i=1:Numbirds
        
        Yvals = {};
        Xbranch = [];
        Nnums = {}; % to collect neurons
        birdname = SummaryStruct.birds(i).birdname;
        BranchnameAll = {};
        %         IsSU = {};
        NumCtxts = {}; % i.e. how many to decode
        
        for ii=1:Numbranch
            
            inds = DatAll.AllBirdNum==i & DatAll.AllBranchNum==ii & ...
                strcmp([DatAll.AllBrainRegion], bregion);
            
            if ~any(inds)
                continue
            end
            
            
            % ==== extract all p-vals for this bird/branch (i.e. across
            % branches)
            Yvals = [Yvals DatAll.AllDecode(inds)];
            Xbranch = [Xbranch ii];
            
            % --- what is the name of this branch?
            tmp = find(inds);
            %         branchname = CLASSES.birds(i).neurons(AllNeurNum(tmp(1))).branchnum(ii).regexprstr;
            branchname = DATSTRUCT_BYBRANCH.bird(i).branchID(ii).regexpstr;
            BranchnameAll = [BranchnameAll branchname];
            
            % ------------ plot text of neurons for each datapoint
            Nnums = [Nnums DatAll.AllNeurNum(inds)];
            
            % ============ for ehac neuron, determine how many contexts
            % there were
            NumCtxts = [NumCtxts DatAll.AllNumCtxts(inds)];
            
        end
        
        if isempty(Yvals)
            continue
        end
        
        
        % =========== plot for this bird
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title(birdname);
        xlabel('branch ID');
        ylabel('log10(prob)');
        
        
        % ------------- plot text of neuron num next to point
        if plottext==1
            for j=1:length(Nnums)
                for jj = 1:length(Nnums{j})
                    %                     suthis = IsSU{j}(jj);
                    neur = Nnums{j}(jj);
                    x = Xbranch(j);
                    y = Yvals{j}(jj);
                    %                     if suthis==1
                    %                         lt_plot_text(x, y, [num2str(neur)], 'b')
                    %                     else
                    lt_plot_text(x, y, [num2str(neur)], [0.6 0.6 0.9])
                    %                     end
                end
            end
        end
        
        % -------- plot data
        lt_plot_MultDist(Yvals, Xbranch, 1, plotcol, 1, 0);
        line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
        ylim([0 1]);
        
        set(gca, 'XTickLabel', BranchnameAll);
        rotateXLabels(gca, 45);
        
        
        % -------- proportion of cases overall significant
        %         yvalsall = cell2mat(Yvals);
        %
        %         numSig = sum(yvalsall(:)<log10(0.05));
        %         numTot = length(yvalsall(:));
        %
        %         lt_plot_annotation(1, [num2str(numSig) '/' num2str(numTot) ' (' num2str(numSig/numTot) ') sig (p<0.05)'], 'r')
        %
        
    end
end

linkaxes(hsplots, 'xy');

end

function [X, FracSig, Nsamp] = fn_Nsigdecode_dsamp(Nall, Nsig)
% Nall = by branch, number of total neurons
% Nsig = by branch, number of neruons with sig decoding

% X = number of neurons retained in dataset
% FracSig = same length as X. number of sylables with at least one neuron
% sig.
% Nsamp = number of syllables used in analysis. same length as X



% ========
        X = 1:max(Nall); % number of neurons in dataset
        FracSig = [];
        Nsamp = []; % number of syls (including those that are downsampled)
        for x = X
            
            % --- find all cases with this many or more number of neurons
            indsyls = find(Nall>=x);
            
            
            % --- for each syl, subsample neurons to x. then ask whether
            % any one of those neurons is significant
            Yissig = [];
            for indtmp = indsyls'
                
                numsig = Nsig(indtmp);
                numtot = Nall(indtmp);
                assert(numtot>=x,' sadfasdf'); 
               
                % ---------- create vector of 1 and 0, where length of
                % vector is number of neurons, and 1 means this neuron
                % signfiicnat
                vectmp = zeros(1, numtot);
                vectmp(1:numsig) = 1;
                
                % --------- subsample that vector to the desired sample
                % size
                vectmp_dsamp = vectmp(randperm(numtot, x));
                
                % -------- if any of the sampled neurons are significant,
                % then mark this branch as significant
                if any(vectmp_dsamp)==1
                    issig=1;
                else
                    issig=0;
                end
                    Yissig = [Yissig; issig];
            end
            
            FracSig = [FracSig; sum(Yissig)/length(Yissig)]; % collect number of 
            % sylalbles showing at least one sig neuron.
            Nsamp = [Nsamp; length(Yissig)];
        end
        
end