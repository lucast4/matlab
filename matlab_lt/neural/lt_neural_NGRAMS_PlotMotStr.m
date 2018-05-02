function lt_neural_NGRAMS_PlotMotStr(OUTSTRUCT, SummaryStruct, plottype, ...
    PairtypesToplot, birdtoplot, neurtoplot)


%% things that shoudl not be changed
% carryover from other code that this borrowed from

usemedian = 0; % only workds for absfrdiff, absfrdiff_globZ or absfrdiff_typediff
plotON=0; % only works for absfrdiff
Nmin = 2; % number of pairs in class. [DOESNT WORK YET]
plotRawGood = 1; % histogram for pos, negative, dat (only works for plottype = absfrdiff_globZ);
dosubtractcontrol = 1; % then subtracts negative control before plotting [if 0, then overlays neg]
sanitycheckuseneg = 0; % uses negative control data instead of data


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
for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    if ~strcmp(birdname, birdtoplot)
        continue
    end
    
    for ii=1:maxneur
        if ii~=neurtoplot
            continue
        end
        inds = All_birdnum==i & All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        allmeans = nan(1,numpairtypes);
        
        if strcmp(plottype, 'oneminusrho')
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
            
            % ======== collect each neuron
            AllPairs_Means = [AllPairs_Means; allmeans];
            AllPairs_Bregions = [AllPairs_Bregions; bregion];
            AllPairs_Birdnum = [AllPairs_Birdnum; i];
            
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
                
                inds = All_birdnum==i & All_neurnum==ii & All_diffsyl_PairType==pairtypethis;
                
                
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
                    YplotAll{3} = [YplotAll{3} All_AbsFRdiff_NEG(inds)]; % collect negative controls
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
            end
            
            % ======== collect each neuron
            AllPairs_Means = [AllPairs_Means; allmeans];
            AllPairs_Bregions = [AllPairs_Bregions; bregion];
            AllPairs_Birdnum = [AllPairs_Birdnum; i];
            
        end
        
        % =============================== PLOT,
        if plotRawGood==1 & strcmp(plottype, 'absfrdiff_globZ')
            
            % rthen can plot!!
            lt_figure; hold on;
            title([birdname '-n' num2str(ii)]);
            
            lt_plot_MultDist(YplotAll, [1 2 3], 0, 'k', 0, '', 1)
            lt_plot_zeroline;
            
            % ----- plot strings
            for k=1:length(YplotAll)
                ntrials = length(YplotAll{k});
                for kk=1:ntrials
                   y = YplotAll{k}(kk);
                   motifpairstr = [MotifPairString{k}{kk,1} '-' MotifPairString{k}{kk,2}];
                   lt_plot_text(k+0.1, y, motifpairstr, 'r', 7);
                   plot(k, y, 'or');
                end
            end
            
        end
        
    end
end