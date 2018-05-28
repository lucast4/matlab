function lt_neural_NGRAMS_Timecourse_Calc(OUTSTRUCT, SummaryStruct, Params, ...
    PairTypesToCompare, tcoursestyle, useGlobalNeg, numrawplotsperbird)

%%


numpairtypes = length(PairTypesToCompare);
maxbirds = max(OUTSTRUCT.All_birdnum);
maxneur = max(OUTSTRUCT.All_neurnum);

pcolors = lt_make_plot_colors(numpairtypes, 0,0);


%%
% figcount=1;
% subplotrows=5;
% subplotcols=3;
% fignums_alreadyused=[];
% hfigs=[];
% hsplots = [];
%
% =====================
% AllPairs_Means = [];

AllPairs_Means = cell(2,1);
AllPairs_Bregions = {};
AllPairs_Birdnum = [];


for i=1:maxbirds
    birdname = SummaryStruct.birds(i).birdname;
    plotdone = 0;
    for ii=1:maxneur
        
        inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
        
        if ~any(inds)
            continue
        end
        
        % -- brainregion
        bregion = SummaryStruct.birds(i).neurons(ii).NOTE_Location;
        allmeans = [];
        %
        %
        YplotAll = cell(1,3); % dat1, dat2, neg control
        %         MotifPairString = cell(1,3);
        
        skipneuron = 0;
        for k =1:numpairtypes
            
            pairtypethis = PairTypesToCompare{k};
            pairtypethis = find(strcmp(OUTSTRUCT.PairTypesInOrder, pairtypethis));
            
            
            % -------- inds for this pairtype
            inds = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii & ...
                OUTSTRUCT.All_diffsyl_PairType==pairtypethis;
            
            if ~any(inds)
                % then skip this neuron, since it is missing at least one
                % of the pairtypes
                skipneuron =1;
                continue
            end
            
            %% ########################################## COLLECT DATA
            % =================== DATA
            % --- extract all data
            tmp = cell2mat(OUTSTRUCT.AllTcourse_FRdiffDAT(inds));
            % -- how many time bins?
            indtmp = find(inds);
            nbins = size(OUTSTRUCT.AllTcourse_FRdiffDAT{indtmp(1)},1); % number of time bins
            % --- reshape to be tbins x ngram
            datmat = reshape(tmp(:,1), nbins, size(tmp,1)/nbins);
            
            
            
            % =============== FIND NEGATIVE CONTROLS FOR THIS NEURON
            if useGlobalNeg==0 % THIS GOOD version, limited to just the desired pairs.
                [~, pairtypes_neg] = intersect(OUTSTRUCT.PairTypesInOrder, PairTypesToCompare);
                indsneg = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii & ...
                    ismember(OUTSTRUCT.All_diffsyl_PairType, pairtypes_neg);
            elseif useGlobalNeg==1
                % using all pairs
                indsneg = OUTSTRUCT.All_birdnum==i & OUTSTRUCT.All_neurnum==ii;
            end
            
            assert(all(~isempty(OUTSTRUCT.AllTcourse_FRdiffShuff(indsneg))), 'why empty?');
            
            
            
            % ==================== NEG CONTROL
            % --- extract all data
            tmp = cell2mat(OUTSTRUCT.AllTcourse_FRdiffShuff(indsneg));
            % -- how many time bins?
            indtmp = find(indsneg);
            nbins = size(OUTSTRUCT.AllTcourse_FRdiffShuff{indtmp(1)},1); % number of time bins
            % --- reshape to be tbins x ngram
            negmat = reshape(tmp(:,1), nbins, size(tmp,1)/nbins);
            
            
            % ###################### CONVERT DATA ALL TRIALS TO ZSCORE
            % =================== take mean and std of negative control
            negmean = mean(negmat, 2);
            negstd = std(negmat,0,2);
            
            % --------------- for data take running zscore
            y = (datmat - negmean)./negstd;
            % ------------- take mean for this neuron
            ymean = mean(y,2);
            
            
            % ============= for output
            allmeans = [allmeans ymean];
            
            
            % ================= OUTPUT FOR PLOTTING
            YplotAll{k} = datmat;
            YplotAll{3} = negmat; % collect negative controls
            
            
            
            %% ==================== COLLECT MOTIF NAMES
            %             motifpairlist = OUTSTRUCT.All_ngramstring_inorder(inds,:);
            %             MotifPairString{k} = motifpairlist;
            %             MotifPairString{3} = [MotifPairString{3}; motifpairlist];
            
            % ###################################### PLOT HISTOGRAMS?
            %             if plotON==1
            %                 % =============== plot histograms
            %                 lt_plot_histogram(y, '', 1, 1, '', 1, pcolors{k});
            %             end
            
            
            % #################### collect mean to then do scatterplot
            %             if usemedian==1
            %                 ymean = median(y);
            %             else
            %                 ymean = mean(y);
            %             end
            %             allmeans(k) = ymean;
        end
        
        
        if skipneuron ==1
            % skip this neuron since missing at least one pairtype
            continue
        end
           
        %% ==================== PLOT
        pcolors = lt_make_plot_colors(3, 1, [0.8 0.2 0.2]);
        if plotdone<numrawplotsperbird & rand<0.5
            figcount=1;
            subplotrows=4;
            subplotcols=2;
            fignums_alreadyused=[];
            hfigs=[];
            hsplots = [];
            
            titlenames = [PairTypesToCompare 'shuff'];
            
            % ########################## PLOT RAW FR DIFFS
            for k=1:3
                
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title(titlenames{k});
                ylabel('fr diff');
                
                % ============ plot individual
                if size(YplotAll{k},2)>200
                    indsplot = randperm(size(YplotAll{k},2), 200);
                else
                    indsplot = 1:size(YplotAll{k},2);
                end
                
                plot(YplotAll{k}(:,indsplot), '-', 'Color', [0.6 0.6 0.6]);
                
                % ========= plot mean
                ymean = mean(YplotAll{k},2);
                ysem = lt_sem(YplotAll{k}');
                shadedErrorBar(1:length(ymean), ymean, ysem, {'Color', 'k'},1)
                
            end
            
            % ########################## overlay all means
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('all means');
            ylabel('mean FR diff');
            for k=1:3
                
                % ========= plot mean
                ymean = mean(YplotAll{k},2);
                ysem = lt_sem(YplotAll{k}');
                shadedErrorBar(1:length(ymean), ymean, ysem, {'Color', pcolors{k}},1);
            end
            line(1000*[Params.regexpr.motifpredur Params.regexpr.motifpredur], ylim);
            
            
            % ============= plot the zscored means
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('all means');
            ylabel('zscore rel shuffle');
            
            for k=1:2
                plot(allmeans(:,k), 'Color', pcolors{k});
            end
            line(1000*[Params.regexpr.motifpredur Params.regexpr.motifpredur], ylim);
            
            
            % ============= plot contextual modualtion index
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title('contxtual modulation ind');
            xlabel([birdname '-n' num2str(ii) bregion]);
            assert(strcmp(PairTypesToCompare{2}, '1  0  0'), 'asfasd');
            assert(strcmp(PairTypesToCompare{1}, '1  1  1'), 'asfasd');
            plot(allmeans(:,2)./allmeans(:,1), '-k');
            lt_plot_zeroline;
            line(xlim, [1 1], 'Color', 'k');
            line(1000*[Params.regexpr.motifpredur Params.regexpr.motifpredur], ylim);
            
            
            % ==============
            plotdone = plotdone+1;
        end
        
        %% ======== collect each neuron
        AllPairs_Means{1} = [AllPairs_Means{1} allmeans(:,1)];
        AllPairs_Means{2} = [AllPairs_Means{2} allmeans(:,2)];
        %         AllPairs_Means = [AllPairs_Means; allmeans];
        AllPairs_Bregions = [AllPairs_Bregions; bregion];
        AllPairs_Birdnum = [AllPairs_Birdnum; i];
    end
end



%% ========== first, use median filter to smooth data




%% ================
% TRY TWO METHODS - 1) GET CMI FOR EACH UNIT, THEN AVERAGE; 2) AVERAGE
% Z-SCORE ACROSS ALL UNITS, THEN CALC MEAN CMI.
% ALSO TRY REGRESSIONS at each time bin

% ======================= MEAN CMI
% CALCULATE CMI FOR EACH UNIT, THEN AVERAGE
maxbirds = max(AllPairs_Birdnum);
bregionstoplot = {'LMAN', 'RA'};

figcount=1;
subplotrows=5;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];


for bregion = bregionstoplot
    for i=1:maxbirds
        
        inds = strcmp(AllPairs_Bregions, bregion) & AllPairs_Birdnum==i;
        if ~any(inds)
            continue
        end
        
        % ================ 1) PLOT ZSCORES
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['b' num2str(i) ',' bregion{1}]);
        ylabel('zscores');
        
        % ------- ind1
        indtmp = 1;
        y = AllPairs_Means{indtmp}(:, inds);
        ymean = mean(y,2);
        ystd = std(y, 0, 2);
        shadedErrorBar(1:length(ymean), ymean, ystd, {'Color', pcolors{indtmp}},1);
        
        
        % ------- ind2
        indtmp = 2;
        y = AllPairs_Means{indtmp}(:, inds);
        ymean = mean(y,2);
        ystd = std(y, 0, 2);
        shadedErrorBar(1:length(ymean), ymean, ystd, {'Color', pcolors{indtmp}},1);
        
        lt_plot_zeroline;
        
        
        % =============== 2) PLOT CMI
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['b' num2str(i) ',' bregion{1}]);
        ylabel('CMI (ratio of means)');
        
        % ------- ind1
        
        cmi = mean(AllPairs_Means{2}(:, inds),2)./mean(AllPairs_Means{1}(:, inds),2);
        
        plot(cmi, '-k');
        
        lt_plot_zeroline;
        plot(xlim, [1 1], 'Color', 'k');
        
        % =============== 2) PLOT CMI
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['b' num2str(i) ',' bregion{1}]);
        ylabel('CMI (mean of ratios)');
        
        % ------- ind1
        cmi = mean(AllPairs_Means{2}(:, inds)./AllPairs_Means{1}(:, inds), 2);
        plot(cmi, '-k');
        
        lt_plot_zeroline;
        plot(xlim, [1 1], 'Color', 'k');
    end
    
    % %%%%%%%%%%%%%%%%%%%%% ACROSS BIRD MEAN (TAKE MEAN OVER ALL UNITS)
    % ============================ 1) MEAN ZSCORE
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all units, ' bregion]);
    ylabel('zscores');
    
    inds = strcmp(AllPairs_Bregions, bregion);
    
    % ------- ind1
    indtmp = 1;
    y = AllPairs_Means{indtmp}(:, inds);
    ymean = mean(y,2);
    ystd = std(y, 0, 2);
    shadedErrorBar(1:length(ymean), ymean, ystd, {'Color', pcolors{indtmp}},1);
    
    
    % ------- ind2
    indtmp = 2;
    y = AllPairs_Means{indtmp}(:, inds);
    ymean = mean(y,2);
    ystd = std(y, 0, 2);
    shadedErrorBar(1:length(ymean), ymean, ystd, {'Color', pcolors{indtmp}},1);
    
    lt_plot_zeroline;
    
    % =============== 2) PLOT CMI
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all units, ' bregion]);
    ylabel('CMI (ratio of means)');
    
    % ------- ind1
    cmi = mean(AllPairs_Means{2}(:, inds),2)./mean(AllPairs_Means{1}(:, inds),2);
    
    plot(cmi, '-k');
    
    lt_plot_zeroline;
    plot(xlim, [1 1], 'Color', 'k');
    
    
    % =============== 2) PLOT CMI
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['all units, ' bregion]);
    ylabel('CMI (mean of ratios)');
    
    % ------- ind1
    cmi = mean(AllPairs_Means{2}(:, inds)./AllPairs_Means{1}(:, inds), 2);
    plot(cmi, '-k');
    
    lt_plot_zeroline;
    plot(xlim, [1 1], 'Color', 'k');
    
end



%% ====================== FOCUS ON ONE METHOD -
% EACH BIRD GET CMI = RATIO OF MEAN(ZSCORES)
% AVERAGE ACROSS ALL BIRDS TO GET ERROR BARS

% ########################## 1) COLLECT ALL CMI
CMIall =[];
Bregionall = {};
Birdnumall = [];

for bregion = bregionstoplot
    for i=1:maxbirds
        
        inds = strcmp(AllPairs_Bregions, bregion) & AllPairs_Birdnum==i;
        if ~any(inds)
            continue
        end
        
        % =========== EXTRACT CMI (mean of zscores)
        cmi = mean(AllPairs_Means{2}(:, inds),2)./mean(AllPairs_Means{1}(:, inds),2);
        
        % ================= SAVE
        CMIall = [CMIall cmi];
        Bregionall = [Bregionall bregion];
        Birdnumall = [Birdnumall i];
    end
    
    
end


%  ############################# 2) PLOT
lt_figure; hold on;
% =============== each bregion
for j=1:length(bregionstoplot)
    bregionthis = bregionstoplot{j};
    
    lt_subplot(4,2,j); hold on;
    title([bregionthis]);
    ylabel('CMI (ratio of mean z)');
    
    inds = strcmp(Bregionall, bregionthis);
    Y = CMIall(:, inds);
    plot(1:size(Y,1), Y);
    
    % ---------
    lt_plot_zeroline;
    ylim([-0.4 1.6]);
    line(xlim, [1 1]);
    line(1000*[Params.regexpr.motifpredur Params.regexpr.motifpredur], ylim);
end

pcolors = lt_make_plot_colors(2, 0, 0);
% ============== 2 bregions overlaid
lt_subplot(4,2,length(bregionstoplot)+1); hold on;
for j=1:length(bregionstoplot)
    bregionthis = bregionstoplot{j};
    inds = strcmp(Bregionall, bregionthis);
    Y = CMIall(:, inds);
    ymean = mean(Y,2);
    ysem = lt_sem(Y');
    shadedErrorBar(1:length(ymean), ymean, ysem, {'Color', pcolors{j}},1);
    
end

% ---------
lt_plot_zeroline;
ylim([-0.4 1.4]);
line(xlim, [1 1]);
line(1000*[Params.regexpr.motifpredur Params.regexpr.motifpredur], ylim);



%% ============== REGRESSION, DO SEPARATELY IN EACH TIME WINDOW
% (NEURON AS UNIT OF DATA)

windsize = 0.015; % sec
windshift = 0.005; %

% ======================= GET LIST OF TEMPORAL WINDOWS
dur_total = Params.regexpr.motifpostdur + Params.regexpr.motifpredur;

onsets = 0:windshift:dur_total-windsize-windshift; % subtract windshift because sometime lose data at the end
offsets = windsize:windshift:dur_total-windshift;
Windowlist = [onsets' offsets'];

xbins = lt_neural_QUICK_XfromFRmat(AllPairs_Means{1}(:,1:2));


%% [REGRESSION] METHOD 1 - SEPARATE MODELS PER BRAIN REGION
usesimpleregression =0; % if 0, then does mixed effects.
formula = 'ydat ~ -1 + xdat + (xdat-1|birdnum)';

% =============== COLLECT
CMIall =[];
CMISEall = [];
Bregionall = {};
% Birdnumall = [];

for bregion = bregionstoplot
    % ==============
    ind_neur = strcmp(AllPairs_Bregions, bregion);
    cmi = [];
    cmi_SE = [];
    
    for j=1:size(Windowlist)
        windthis = Windowlist(j,:);
        
        ind_x = xbins>=windthis(1) & xbins<windthis(2);
        
        xdat = double(mean(AllPairs_Means{1}(ind_x, ind_neur),1)');
        ydat = double(mean(AllPairs_Means{2}(ind_x, ind_neur),1)');
        birdnum = AllPairs_Birdnum(ind_neur);
        
        % ============= quick - simple linaer regression
        if usesimpleregression==1
            [b,bint] = lt_regress(ydat, xdat, 0);
            slope = b(2);
        else
            tbl = table(xdat, ydat, birdnum);
            lme = fitlme(tbl, formula);
            indcoeff = strcmp(lme.CoefficientNames, 'xdat');
            slope = lme.Coefficients.Estimate(indcoeff);
            slope_SE = lme.Coefficients.SE(indcoeff);
            
        end
        % ============ collect
        cmi = [cmi slope];
        cmi_SE = [cmi_SE slope_SE];
    end
    CMIall = [CMIall; cmi];
    CMISEall = [CMISEall; cmi_SE];
    Bregionall = [Bregionall; bregion{1}];
end


% ---------------------
lt_figure; hold on;
if usesimpleregression==1
ylabel('cmi, from simple regression, pooled over units');
else
    ylabel('cmi, from fitting separate models for each bregion');
end
for j=1:size(CMIall,1)
    
    x = mean(Windowlist,2);
    y = CMIall(j,:);
    ySE = CMISEall(j,:);
    bregion = Bregionall{j};
    
%     plot(x,y);
    shadedErrorBar(x, y, ySE, {'Color', 'k'}, 1);
end

%%
% AllPairs_Bregions_backup = AllPairs_Bregions;
% 
% indtmp = strcmp(AllPairs_Bregions, 'RA');
% 
% AllPairs_Bregions(indtmp) = 'LMAN';
% AllPairs_Bregions(~indtmp) = 'RA';
%% [REGRESSION] METHOD 2- one model for all datapoints 
% (fit effect of bregion on slope)
% still fits separate models for each time point.
formula = 'ydat ~ -1 + xdat + xdat:bregions + (xdat-1|birdnum)';

allslopes = [];
allslopes_SE = [];
allslopediff = [];
allslopediff_SE = [];
for j=1:size(Windowlist)
    windthis = Windowlist(j,:);
    
    ind_x = xbins>=windthis(1) & xbins<windthis(2);
    
    xdat = double(mean(AllPairs_Means{1}(ind_x, :),1)');
    ydat = double(mean(AllPairs_Means{2}(ind_x, :),1)');
    birdnum = AllPairs_Birdnum;
    bregions = AllPairs_Bregions;
    
    % ============== RUN MODEL
    tbl = table(xdat, ydat, birdnum, bregions);
    lme = fitlme(tbl, formula);
    
    % --------- slope
    slope = lme.Coefficients.Estimate(1);
    slope_SE = lme.Coefficients.SE(1);
    
    % --------- bregion:slope
    slopediff = lme.Coefficients.Estimate(2);
    slopediff_SE = lme.Coefficients.SE(2);
    
    % ============ collect
    allslopes = [allslopes slope];
    allslopes_SE = [allslopes_SE slope_SE];
    allslopediff = [allslopediff slopediff];
    allslopediff_SE = [allslopediff_SE slopediff_SE];
    
end



% ------------------------ PLOT
lt_figure; hold on;
    x = mean(Windowlist,2);
% 
% plot(x, allslopes, '-k');
% plot(x, allslopes+allslopediff, '-r');

shadedErrorBar(x, allslopes, allslopes_SE, {'Color', 'k'},1);
shadedErrorBar(x, allslopes+allslopediff, allslopediff_SE, {'Color', 'r'},1);

lt_plot_zeroline;
ylim([-0.4 1.4]);
line(xlim, [1 1]);
line([Params.regexpr.motifpredur Params.regexpr.motifpredur], ylim);


















