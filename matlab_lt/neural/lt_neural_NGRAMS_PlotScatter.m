function [AllPairs_Means, AllPairs_Birdnum, AllPairs_Bregions] = ...
    lt_neural_NGRAMS_PlotScatter(OUTSTRUCT, SummaryStruct, plottype, plotON, ...
    PairtypesToplot, dosubtractcontrol, sanitycheckuseneg, plotRawGood, usemedian, ...
    removeBadSyls, minPairs, zscoreMin, labelneur)
%%
assert(sum([plotRawGood plotON])<2, 'cannot do both types of raw plots ...');

%% for globalz
useGlobalNeg = 0; % if 0, then just neg for the desired pairs. default: 1

%% LT 4/24/18 - plots scatter (i.e. mean for pairtype 2 vs. 1)


%% INPUTS

% plottype = 'absfrdiff';
% oneminusrho - 1 minus correlation coefficient of smoothed FR
% absfrdiff - diff in FR, each pair/unit normalized (z) versus its own
% shuffle control
% absfrdiff_globZ - mean abs diff in FR, normalized versus global
% distribtuion of neg control shuffles

% plotON=0; % only works for absfrdiff
% PairtypesToplot = {...
%     '1  1  1', ... % xaxis
%     '1  0  0'}; % yaxis
%
% % ----------- params for one minus rho, specifically
% dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
% sanitycheckuseneg = 0; % uses negative control data instead of data


%% pull out data
% ========================== FOR COMPATIBILITY WITH OLD CODE, EXTRACT ALL
% FIELDS
fnamesthis = fieldnames(OUTSTRUCT);
for j=1:length(fnamesthis)
    eval([fnamesthis{j} ' = OUTSTRUCT.' fnamesthis{j} ';']);
end

%%
numpairtypes = length(PairtypesToplot);
assert(numpairtypes == 2, 'assumes 2, otherwise code below doesnt work (e.g. assumes that YplotAll{3} is negative control');
maxbirds = max(All_birdnum);
maxneur = max(All_neurnum);

pcolors = lt_make_plot_colors(numpairtypes, 0,0);


%%
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% =====================
AllPairs_Means = [];
AllPairs_Bregions = {};
AllPairs_Birdnum = [];
AllPairs_Neur = [];
for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    for ii=1:maxneur
        
        inds = All_birdnum==i & All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        allmeans = nan(1,numpairtypes);
        Neach = nan(1,numpairtypes);
        if strcmp(plottype, 'oneminusrho')
            disp('NOT IMPLEMENTTING NMIN YET (for oneminusrho) !! easy to modify code./..');
            % ========= plot
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-n' num2str(ii) '[' bregion ']']);
            
            AllRhoNeg = [];
            for k =1:numpairtypes
                
                pairtypethis = PairtypesToplot{k};
                pairtypethis = find(strcmp(PairTypesInOrder, pairtypethis));
                
                inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
                
                % ========= collect data
                pairtype = All_diffsyl_PairType(inds);
                oneminusrho_neg = All_OneMinusRho_NEG(inds);
                oneminusrho = All_OneMinusRho(inds);
                
                if sanitycheckuseneg==1
                    oneminusrho = oneminusrho_neg;
                end
                
                % ----- subtract negative control?
                if dosubtractcontrol==1
                    oneminusrho = oneminusrho - oneminusrho_neg;
                end
                
                % =============== plot histograms
                xcenters = 0.05:0.1:1.95;
                if dosubtractcontrol==1
                    xcenters = -1.15:0.1:1.95;
                end
                lt_plot_histogram(oneminusrho, xcenters, 1, 1, '', 1, pcolors{k});
                
                
                % ================ collect mean to then do scatterplot
                ymean = mean(oneminusrho);
                allmeans(k) = ymean;
                
                % ======== collect negative distribution
                AllRhoNeg = [AllRhoNeg; oneminusrho_neg];
            end
            
            % ======== overlay negative distribution
            if dosubtractcontrol==0
                lt_plot_histogram(AllRhoNeg, xcenters, 1, 1, '', 1, [0.7 0.7 0.7]);
            end
            
            
        else
            %%
            if plotON==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([birdname '-n' num2str(ii) '[' bregion ']']);
            end
            
            YplotAll = cell(1,3); % dat1, dat2, neg control
            MotifPairString = cell(1,3);
            for k =1:numpairtypes
                
                pairtypethis = PairtypesToplot{k};
                pairtypethis = find(strcmp(PairTypesInOrder, pairtypethis));
                
                
                if removeBadSyls==1
                    inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis & ...
                        OUTSTRUCT.All_BadSyls == 0;
                else
                    inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
                end
                
                % ########################################## COLLECT DATA
                y = [];
                if strcmp(plottype, 'absfrdiff')
                    % then mean abs diff in FR, normalized within each
                    % pairtype
                    y = All_AbsFRdiff_Zrelshuff(inds);
                    
                elseif strcmp(plottype, 'absfrdiff_globZ')
                    % =========== version with both dat and pos compared to
                    % distribution of negative controls
                    
                    % ----- V1 - NEG = combined from the pairs being analyzed
                    if useGlobalNeg==0 % old version, limited to just the desired pairs. but I think
                        % is better to use all pairs ...
                        [~, pairtypes_neg] = intersect(PairTypesInOrder, PairtypesToplot);
                        indsneg = All_birdnum==i & All_neurnum==ii & ...
                            ismember(All_diffsyl_PairType, pairtypes_neg);
                    elseif useGlobalNeg==1
                        indsneg = All_birdnum==i & All_neurnum==ii;
                    end
                    
                    
                    negmean = mean(All_AbsFRdiff_NEG(indsneg));
                    negstd = std(All_AbsFRdiff_NEG(indsneg));
                    
                    y = (All_AbsFRdiff(inds) - negmean)./negstd;
                    
                    % ================= OUTPUT FOR PLOTTING
                    YplotAll{k} = All_AbsFRdiff(inds);
                    YplotAll{3} = [YplotAll{3}; All_AbsFRdiff_NEG(inds)]; % collect negative controls
                elseif strcmp(plottype, 'absfrdiff_typediff')
                    
                    % ----- V2 - each type compaired to mean of its own neg
                    % (not zscored)
                    y = All_AbsFRdiff(inds) - mean(All_AbsFRdiff_NEG(inds));
                    
                end
                
                % ==================== COLLECT MOTIF NAMES
                motifpairlist = OUTSTRUCT.All_ngramstring_inorder(inds,:);
                MotifPairString{k} = motifpairlist;
                MotifPairString{3} = [MotifPairString{3}; motifpairlist];
                
                % ###################################### PLOT HISTOGRAMS?
                if plotON==1
                    % =============== plot histograms
                    lt_plot_histogram(y, '', 1, 1, '', 1, pcolors{k});
                end
                
                
                % #################### collect mean to then do scatterplot
                if usemedian==1
                    ymean = median(y);
                else
                    ymean = mean(y);
                end
                allmeans(k) = ymean;
                Neach(k) = length(y);
            end
            
            if any(Neach<minPairs)
                disp('SKIPPING - too few pairs...');
                continue
            end
            
            if allmeans(zscoreMin(1))<zscoreMin(2)
                disp('SKIP, signal lower than zscoreMin');
                continue
            end
            
            
        end
            % ======== collect each neuron
            AllPairs_Means = [AllPairs_Means; allmeans];
            AllPairs_Bregions = [AllPairs_Bregions; bregion];
            AllPairs_Birdnum = [AllPairs_Birdnum; i];
            AllPairs_Neur = [AllPairs_Neur; ii];
        
        % =============================== PLOT,
        if plotRawGood==1 & strcmp(plottype, 'absfrdiff_globZ')
            
            % rthen can plot!!
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname '-n' num2str(ii) '[' bregion ']']);
            
            lt_plot_MultDist(YplotAll, [1 2 3], 0, 'k', 0)
            lt_plot_zeroline;
            
            if (0) % old version, plotting histograms , but is too crowded...
                
                xcenters = min(cell2mat(YplotAll)):0.1:max(cell2mat(YplotAll));
                
                % ----- NEG
                indtmp = 3;
                pcol = [0.7 0.6 0.7];
                lt_plot_histogram(YplotAll{indtmp}, xcenters, 1, 1, '', 1, pcol);
                % ----- POS
                indtmp = 1;
                pcol = [0.3 0.3 0.9];
                lt_plot_histogram(YplotAll{indtmp}, xcenters, 1, 1, '', 1, pcol);
                % ----- DAT
                indtmp = 2;
                pcol = 'k';
                lt_plot_histogram(YplotAll{indtmp}, xcenters, 1, 1, 0.2, 0, pcol);
            end
        end
        
    end
end


%% =========== PLOT SCATTER COMPARING TWO CLASSES ACROSS ALL NEURONS
if strcmp(plottype, 'oneminusrho')
    lt_figure; hold on;
    % ----- LMAN
    lt_subplot(2,2,1); hold on;
    title('LMAN');
    indstmp = strcmp(AllPairs_Bregions, 'LMAN');
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ok');
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-1 1], [-1 1]);
    
    % ----- RA
    lt_subplot(2,2,2); hold on;
    title('RA');
    indstmp = strcmp(AllPairs_Bregions, 'RA');
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ok');
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-1 1], [-1 1]);
    
    % --- COMBINED
    lt_subplot(2,2,3); hold on;
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ob');
    
    % RA
    indstmp = strcmp(AllPairs_Bregions, 'RA');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'or');
    
    % FORMATING
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    xlim([-1 1]);
    ylim([-1 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-1 1], [-1 1]);
    
    
    % ============ SAPRATE BY BIRD
    maxbirds = max(AllPairs_Birdnum);
    for j=1:maxbirds
        lt_figure; hold on;
        
        % LMAN
        indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'ob');
        
        % RA
        indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'or');
        
        % FORMATING
        xlabel(PairtypesToplot{1});
        ylabel(PairtypesToplot{2});
        
        xlim([-1 1]);
        ylim([-1 1]);
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
        line([-1 1], [-1 1]);
        
        
        
    end
else
    
    % ================== ONE, EQUAL AXES
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    maxbirds = max(AllPairs_Birdnum);
    for j=1:maxbirds
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([SummaryStruct.birds(j).birdname]);
        hsplots = [hsplots hsplot];
        % LMAN
        indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'ob');
        
        % RA
        indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'or');
        
        % FORMATING
        xlabel(PairtypesToplot{1});
        ylabel(PairtypesToplot{2});
        
        %         xlim([-1 1]);
        %         ylim([-1 1]);
        
        lt_plot_makesquare_plot45line(gca, 'k', -2);
    end
    
    
    % ============ ONE PLOT ALL BIRDS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['ALL BIRDS']);
    hsplots = [hsplots hsplot];
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ob');
    
    % RA
    indstmp = strcmp(AllPairs_Bregions, 'RA');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'or');
    
    % FORMATING
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    %         xlim([-1 1]);
    %         ylim([-1 1]);
    
    lt_plot_makesquare_plot45line(gca, 'k', -2);
    linkaxes(hsplots, 'xy');
    
    
    % ############################## TWO, SEPARATE AXES
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    maxbirds = max(AllPairs_Birdnum);
    for j=1:maxbirds
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([SummaryStruct.birds(j).birdname]);
        % LMAN
        indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'ob');
        
        % RA
        indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
        y = AllPairs_Means(indstmp, :);
        plot(y(:,1), y(:,2), 'or');
        
        % FORMATING
        xlabel(PairtypesToplot{1});
        ylabel(PairtypesToplot{2});
        
        %         xlim([-1 1]);
        %         ylim([-1 1]);
        
        lt_plot_makesquare_plot45line(gca, 'k', -2);
    end
    
    
    % ============ ONE PLOT ALL BIRDS
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['ALL BIRDS']);
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'ob');
    
    % RA
    indstmp = strcmp(AllPairs_Bregions, 'RA');
    y = AllPairs_Means(indstmp, :);
    plot(y(:,1), y(:,2), 'or');
    
    % FORMATING
    xlabel(PairtypesToplot{1});
    ylabel(PairtypesToplot{2});
    
    %         xlim([-1 1]);
    %         ylim([-1 1]);
    
    lt_plot_makesquare_plot45line(gca, 'k', -2);
    
end


%% ========== ancova,
if (0)
    indstmp = AllPairs_Birdnum==6;
    y1 = AllPairs_Means(indstmp, 1);
    y2 = AllPairs_Means(indstmp, 2);
    group = AllPairs_Bregions(indstmp);
    
    aoctool(y1, y2, group)
end

%% ====== OVERLAY DISTRIBUTIONS FOR ALL BIRDS
lt_figure; hold on;
title('each bird one fill col (circle(RA), sqiare(L)');
maxbirds = max(AllPairs_Birdnum);
plotcols = lt_make_plot_colors(maxbirds, 0,0);
xlabel(PairtypesToplot{1});
ylabel(PairtypesToplot{2});
for j=1:maxbirds
    
    % ======= fill with color?
    if length(unique(AllPairs_Bregions(AllPairs_Birdnum==j)))>1
        % then color, since multiple brain regions
        pcol = plotcols{j};
    else
        pcol = [1 1 1];
    end
    
    % ============ LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
    y = AllPairs_Means(indstmp, :);
    % get mean and sem in 2 dimensions
    ymean = mean(y,1);
    ysem_x = lt_sem(y(:,1));
    ysem_y = lt_sem(y(:,2));
    
    plot(ymean(1), ymean(2), 's', 'MarkerFaceColor', pcol, ...
        'Color', 'b', 'MarkerSize', 8);
    
    % =============== RA
    indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
    y = AllPairs_Means(indstmp, :);
    % get mean and sem in 2 dimensions
    ymean = mean(y,1);
    ysem_x = lt_sem(y(:,1));
    ysem_y = lt_sem(y(:,2));
    
    plot(ymean(1), ymean(2), 'o','MarkerFaceColor', pcol, ...
        'Color', 'r', 'MarkerSize', 8);
    
    lt_plot_makesquare_plot45line(gca, 'k', -2);
end

%% ================= GET RATIO OF MEANS FOR EACH BIRD,
maxbirds = max(AllPairs_Birdnum);

Y = nan(maxbirds,2);

for i=1:maxbirds
    
    % ========= LMAN
    x = 1;
    inds = AllPairs_Birdnum==i & strcmp(AllPairs_Bregions, 'LMAN');
    
    % ---- get ratio of means
    y = mean(AllPairs_Means(inds,:), 1);
    Y(i, x) = y(2)/y(1);
    
    
    % ========= RA
    x = 2;
    inds = AllPairs_Birdnum==i & strcmp(AllPairs_Bregions, 'RA');
    
    % ---- get ratio of means
    y = mean(AllPairs_Means(inds,:), 1);
    Y(i, x) = y(2)/y(1);
    
end

lt_figure; hold on;
title('ratio of means, by bird');
xlabel('LMAN -- RA');
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}]);

lt_plot([1 2], Y);

% ---- plot means


% --- plot with lines
plot([1 2], Y(all(~isnan(Y),2),:)', '-k');

xlim([0 3]);

p = ranksum(Y(:,1), Y(:,2));
lt_plot_pvalue(p, 'ranksum', 1);
lt_plot_zeroline;



%% ================= FOR EACH NEURON, DIVIDE DATA INTO POS CONTROL

AllPairs_Ratios = AllPairs_Means(:,2)./AllPairs_Means(:,1);

% ====================== PLOT, each bird compare LMAN and RA
figcount=1;
subplotrows=5;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


maxbirds = max(AllPairs_Birdnum);
for j=1:maxbirds
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([SummaryStruct.birds(j).birdname]);
    hsplots = [hsplots hsplot];
    
    Y = {};
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
    Y{1} = AllPairs_Ratios(indstmp, :);
    
    % LMAN
    indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
    Y{2} = AllPairs_Ratios(indstmp, :);
    
    % === plot
    lt_plot_MultDist(Y, [1 2], 1);
    
    % FORMATING
    ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}])
    xlabel('LMAN -- RA');
    xlim([0 3]);
    
end

% ============ grand mean
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('ALL BIRDS');
hsplots = [hsplots hsplot];

Y = {};
% LMAN
indstmp = strcmp(AllPairs_Bregions, 'LMAN');
Y{1} = AllPairs_Ratios(indstmp, :);

% LMAN
indstmp = strcmp(AllPairs_Bregions, 'RA');
Y{2} = AllPairs_Ratios(indstmp, :);

% === plot
lt_plot_MultDist(Y, [1 2], 1);

% FORMATING
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}])
xlabel('LMAN -- RA');
xlim([0 3]);


linkaxes(hsplots, 'y');


%% ###########################33 COMBINE ALL IN ONE PLOT
lt_figure; hold on;
title('[dat - neg]/[pos - neg]');
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}])
xlabel('Birds (LMAN -- RA)');

maxbirds = max(AllPairs_Birdnum);
Yall = cell(maxbirds, 2);
for j=1:maxbirds
    
    % LMAN
    pcol = 'b';
    x = j*3-2;
    indstmp = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==j;
    if any(indstmp)
        lt_plot_MultDist({AllPairs_Ratios(indstmp, :)}, x, 0, pcol);
        if labelneur==1
           z = AllPairs_Neur(indstmp);
           y = AllPairs_Ratios(indstmp, :);
           for k=1:length(y)
              lt_plot_text(x+0.3, y(k), ['n=' num2str(z(k))], 'm', 8);
           end
        end
    end
    Yall{j, 1} = AllPairs_Ratios(indstmp, :);
    
    % RA
    pcol = 'r';
    x = j*3-1;
    indstmp = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==j;
    if any(indstmp)
        lt_plot_MultDist({AllPairs_Ratios(indstmp, :)}, x, 0, pcol);
            if labelneur==1
           z = AllPairs_Neur(indstmp);
           y = AllPairs_Ratios(indstmp, :);
           for k=1:length(y)
              lt_plot_text(x+0.3, y(k), ['n=' num2str(z(k))], 'm', 8);
           end
        end
end
    Yall{j, 2} = AllPairs_Ratios(indstmp, :);
    
    line(3*[j j], ylim, 'Color', 'k' ,'LineStyle', '--');
end
lt_plot_zeroline;
set(gca, 'XTick', 1:3:3*maxbirds, 'XTickLabel', {SummaryStruct.birds.birdname});
% xticks(1:3:3*maxbirds);
% xticklabels({SummaryStruct.birds.birdname});
rotateXLabels(gca, 90);


%% ########################### PLOT ALL ON ONE PLOT
lt_figure; hold on;
% === bird by bird
lt_subplot(2,2,1); hold on;
title('by bird (mean, std)');
xlabel('LMAN -- RA');
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}]);
for j=1:size(Yall,1)
    
    xplot = [1 3] + 0.8*rand-0.4;
    
    % --lman
    x = 1;
    ymean = mean(Yall{j,x});
    ystd = std(Yall{j,x});
    lt_plot(xplot(x), ymean, {'Errors', ystd});
    
    % -- ra
    x = 2;
    ymean = mean(Yall{j,x});
    ystd = std(Yall{j,x});
    lt_plot(xplot(x), ymean, {'Errors', ystd});
    
    if all(~cellfun(@isempty, Yall(j,:)))
        
        % -- plot line between
        ymeans = cellfun(@mean, Yall(j,:));
        plot(xplot, ymeans, '-or');
    end
end
xlim([0 4]);
lt_plot_zeroline;
y1 = cellfun(@mean, Yall(:,1));
y2 = cellfun(@mean, Yall(:,2));

y1 = y1(~isnan(y1));
y2 = y2(~isnan(y2));

p = ranksum(y1, y2);
lt_plot_pvalue(p, 'rank sum', 1);



% === all datapoints
lt_subplot(2,2,2); hold on;
title('by neuron');
xlabel('LMAN -- RA');
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}]);

Ypoints  = {};
Ypoints{1} = cell2mat(Yall(:,1));
Ypoints{2} = cell2mat(Yall(:,2));

lt_plot_MultDist(Ypoints, [1 2], 0, 'k');
lt_plot_zeroline;

% --- pvalue
p = ranksum(Ypoints{1}, Ypoints{2});
lt_plot_pvalue(p, 'rank sum', 1);


%% ============= PLOT ALL DATA AND OVERLAY BIRDS

% ########################### PLOT ALL ON ONE PLOT
lt_figure; hold on;


lt_subplot(3,2,1); hold on;
% === all datapoints
xlabel('LMAN -- RA');
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}]);

Ypoints  = {};
Ypoints{1} = cell2mat(Yall(:,1));
Ypoints{2} = cell2mat(Yall(:,2));

% lt_plot_MultDist(Ypoints, [1 3], 0, 'k');
% lt_plot_MultDist(Ypoints, [1 3], 0, 'k', 0, 1);
lt_plot_MultDist(Ypoints, [1 3], 1, 'b', 0, 0, 1);
% lt_plot_MultDist(Ypoints, [1 3], 1, 'b', 1, 0)
lt_plot_zeroline;

% --- pvalue
p = ranksum(Ypoints{1}, Ypoints{2});
lt_plot_pvalue(p, 'rank sum', 1);

% === bird by bird
for j=1:size(Yall,1)
    
    xplot = [1 3] + 0.5*rand-0.25;
    
    % --lman
    x = 1;
    ymean = mean(Yall{j,x});
    ystd = std(Yall{j,x});
    lt_plot(xplot(x), ymean, {'Errors', ystd});
    
    % -- ra
    x = 2;
    ymean = mean(Yall{j,x});
    ystd = std(Yall{j,x});
    lt_plot(xplot(x), ymean, {'Errors', ystd});
    
    if all(~cellfun(@isempty, Yall(j,:)))
        
        % -- plot line between
        ymeans = cellfun(@mean, Yall(j,:));
        plot(xplot, ymeans, '-or');
    end
end
xlim([0 4]);
lt_plot_zeroline;
y1 = cellfun(@mean, Yall(:,1));
y2 = cellfun(@mean, Yall(:,2));

y1 = y1(~isnan(y1));
y2 = y2(~isnan(y2));

p = ranksum(y1, y2);
lt_plot_pvalue(p, 'rank sum', 1);



lt_subplot(3,2,2); hold on;
% === all datapoints
xlabel('LMAN -- RA');
ylabel([PairtypesToplot{2} '/' PairtypesToplot{1}]);

Ypoints  = {};
Ypoints{1} = cell2mat(Yall(:,1));
Ypoints{2} = cell2mat(Yall(:,2));

% lt_plot_MultDist(Ypoints, [1 3], 0, 'k');
% lt_plot_MultDist(Ypoints, [1 3], 0, 'k', 0, 1);
% lt_plot_MultDist(Ypoints, [1 3], 1, 'b', 0, 0, 1);
% lt_plot_MultDist(Ypoints, [1 3], 1, 'b', 1, 0)
lt_plot_bar([1 3], cellfun(@mean, Ypoints), {'Errors', cellfun(@lt_sem, Ypoints)});
lt_plot_zeroline;

% --- pvalue
p = ranksum(Ypoints{1}, Ypoints{2});
lt_plot_pvalue(p, 'rank sum', 1);

% === bird by bird
for j=1:size(Yall,1)
    
    xplot = [1 3] + 0.5*rand-0.25;
    
    % --lman
    x = 1;
    ymean = mean(Yall{j,x});
    ystd = std(Yall{j,x});
    lt_plot(xplot(x), ymean, {'Errors', ystd});
    
    % -- ra
    x = 2;
    ymean = mean(Yall{j,x});
    ystd = std(Yall{j,x});
    lt_plot(xplot(x), ymean, {'Errors', ystd});
    
    if all(~cellfun(@isempty, Yall(j,:)))
        
        % -- plot line between
        ymeans = cellfun(@mean, Yall(j,:));
        plot(xplot, ymeans, '-or');
    end
end
xlim([0 4]);
lt_plot_zeroline;
y1 = cellfun(@mean, Yall(:,1));
y2 = cellfun(@mean, Yall(:,2));

y1 = y1(~isnan(y1));
y2 = y2(~isnan(y2));

p = ranksum(y1, y2);
lt_plot_pvalue(p, 'rank sum', 1);



%% ======================= mixed effects model

Yall; % row is bird, column is brain region

% -- expand
bnum = [];
bregion = [];
cmi = [];
for i=1:size(Yall,1)
    
    for ii=1:size(Yall,2)
        
        ythis = Yall{i, ii};
        
        cmi = [cmi; ythis];
        bnum = [bnum; i*ones(size(ythis))];
        bregion = [bregion; ii*ones(size(ythis))];
    end
    
end

% ---- make categorical
bnum = categorical(bnum);
% bnum = flipud(bnum);
bregion = categorical(bregion);

lt_figure; hold on;
plot(bregion, cmi, 'ok');

dat = table(cmi, bnum, bregion);
% formula = 'cmi ~ bregion';
formula = 'cmi ~ bregion + (bregion|bnum)';
% formula = 'cmi ~ bnum';
lme = fitlme(dat, formula)

%% ======================== PLOT EACH CLASS SEPARATELY ACROSS BIRDS

Ymeans = nan(maxbirds,4); % bird x [4 classes below]
Ystd = nan(maxbirds, 4);

for i=1:maxbirds
    % ------ LMAN, TYPE 1
    inds = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==i;
    indtype = 1;
    % get
    Ymeans(i, 1) = mean(AllPairs_Means(inds, indtype));
    Ystd(i, 1) = std(AllPairs_Means(inds, indtype));
    
    % ------ LMAN TYPE 2
    inds = strcmp(AllPairs_Bregions, 'LMAN') & AllPairs_Birdnum==i;
    indtype = 2;
    % get
    Ymeans(i, 2) =mean(AllPairs_Means(inds, indtype));
    Ystd(i, 2) = std(AllPairs_Means(inds, indtype));
    
    % ----- RA, TYPE 1
    inds = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==i;
    indtype = 1;
    % get
    Ymeans(i, 3) = mean(AllPairs_Means(inds, indtype));
    Ystd(i, 3) = std(AllPairs_Means(inds, indtype));
    
    % ----- RA, TYPE 2
    inds = strcmp(AllPairs_Bregions, 'RA') & AllPairs_Birdnum==i;
    indtype = 2;
    % get
    Ymeans(i, 4) = mean(AllPairs_Means(inds, indtype));
    Ystd(i, 4) = std(AllPairs_Means(inds, indtype));
end

% ======== plot
lt_figure; hold on;
title('each bird mean,std');
for i = 1:maxbirds
    
    x = [1:2:7] + 0.8*rand-0.4;
    y = Ymeans(i,:);
    ystd = Ystd(i,:);
    
    lt_plot(x, y, {'Errors', ystd});
    if all(~isnan(Ymeans(i,:)))
        plot(x,y, '-', 'Color', [0.2 0.2 0.8]);
    else
        plot(x,y, '-', 'Color', [0.7 0.7 0.7]);
    end
end
lt_plot([1:2:7]+0.5, nanmean(Ymeans,1), {'Errors', lt_sem(Ymeans), 'Color', 'r'});


lt_plot_zeroline;
xlabel(['LMAN(' PairtypesToplot{1} ') LMAN(' PairtypesToplot{2} ') RA(type1) RA(type2)']);
ylabel('frdiff, z rel shuffle');
