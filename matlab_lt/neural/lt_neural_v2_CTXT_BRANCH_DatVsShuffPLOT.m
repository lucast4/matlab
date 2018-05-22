function [DATSTRUCT, DATSTRUCT_BYBRANCH] = lt_neural_v2_CTXT_BRANCH_DatVsShuffPLOT(analyfname, ...
    plotON)
if ~exist('plotON', 'var')
    plotON =1;
end

%% lt 11/29/17 - modified to not use CLASSES (Which takes up al lot of memory ...)


%% lt 10/27/17 - plots CLASSES, after runninglt_neural_v2_CTXT_BRANCH_DatVsShuff
% - plots decoding of dat vs. shuffled
% - currently only works with one time bin (the first) - can easily modify
% to plot multiple time bins.

%% params

plottext = 0; % neuron number

%%

decodestat = 'F1';
plotdecodenegdistr = 0; % will plot each neuron/branch neg and actual dat (lots of plots!!!)

%% load branch
savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

% load([savedir '/CLASSESv2_' analyfname '.mat']);
load([savedir '/ALLBRANCHv2_' analyfname '.mat']);
load([savedir '/SUMMARYv2_' analyfname '.mat']);
try
    load([savedir '/PARAMSv2_' analyfname '.mat']);
catch err
end

%% ======= Filter data (e.g. remove noise, poor labels, etc)

Params.LocationsToKeep = {};
Params.birdstoexclude = {};
Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
Params.durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
Params.durThreshOmega.gappre= [];
Params.durThreshOmega.gappost= [];
Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
% class within branch) is longer than this (sec)
Params.RemoveHandCoded =1 ; % see below

ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params);

%% use first time bin only (can modify)

tt = 1;
apos = 1; % align pos = 1

%% ========= COMPUILE BY ACTUAL BRANCH IDS (based on regexp str)
DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, ...
    analyfname, tt);

disp('MANY EXTRACTION FAILURES BECUASE BRANCH POINTS DOESNT EXIST!! - IS OK')

%% 

numshuffs = length(DATSTRUCT_BYBRANCH.bird(1).branchID(1).DAT.PREMOTORDECODE_struct(1).ConfMatAll_NEG);

%% COLLECT


% ---- to visualize decodeneg distributions
figcount=1;
subplotrows=8;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];


AllPdat = [];
AllBirdNum = [];
AllNeurNum = [];
AllBranchNum = [];
AllDecode = [];
AllDecode_z = [];

AllWindows_RelOnsActual = [];
AllWindows_RelOnsOffDesired = [];

AllBrainRegion = {};
DATSTRUCT.AllSingleUnit = [];

DATSTRUCT.AllNumCtxts = [];

DATSTRUCT.AllDecode_neg_mean = [];

Numbirds = length(DATSTRUCT_BYBRANCH.bird);

% ------------ load decode param
timewindows = load([savedir '/' analyfname '/SHUFFDECODE/Params.mat']);

for i=1:Numbirds
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb = 1:numbranches
        
        numneurons = size(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode,1);
        
        for nn=1:numneurons
            % ------------
            disp(['brd' num2str(i) '-br' num2str(bb) '-neur' num2str(nn)]);
            
            % ------------ find shuffle decode data
            decodestruct = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_struct(nn);
            p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(nn);
            brainregion = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion{nn};
            
            neurnum = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(nn);
            
            % === other features
            AllBirdNum = [AllBirdNum; i];
            AllNeurNum = [AllNeurNum; neurnum];
            AllBranchNum = [AllBranchNum; bb];
            AllBrainRegion = [AllBrainRegion; brainregion];
                        
            % === Pdats
            AllPdat = [AllPdat; p];
            
            % === num contexts
            nctxts = size(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_struct(nn).ConfMatAll_DAT{1}, 1);
            DATSTRUCT.AllNumCtxts = [DATSTRUCT.AllNumCtxts; nctxts];
            
            % === convert decode into z-score rel to null distribution
            cmat = decodestruct.ConfMatAll_NEG;
            n = min([100 length(cmat)]); % cap num trial for z-score at 50.
            decodeneg = [];
            for j=1:n
                sts = lt_neural_ConfMatStats(cmat{j});
                decodeneg = [decodeneg sts.(decodestat)];
            end
            
            % --- same, data
            cmat = decodestruct.ConfMatAll_DAT;
            n = min([100 length(cmat)]); % cap num trial for z-score at 50.
            decodedat = [];
            for j=1:n
                sts = lt_neural_ConfMatStats(cmat{j});
                decodedat = [decodedat sts.(decodestat)];
            end
            decode = mean(decodedat);
            
            
            % ============= is this SU?
            singleunitLogical = ...
                strcmp({SummaryStruct.birds(i).neurons(neurnum).NOTE_is_single_unit}, 'yes');
            DATSTRUCT.AllSingleUnit = [DATSTRUCT.AllSingleUnit; ...
                singleunitLogical];
            
            
            % ################################# OUTPUT
            AllDecode = [AllDecode decode];
            
            nullmean = mean(decodeneg);
            nullstd = std(decodeneg);
            AllDecode_z = [AllDecode_z (decode - nullmean)/nullstd];
            
            
            % --------------- MEAN DECODE NEG, cap N at number of
            % iterations for dat (to equalize variance)
            decode_neg_mean = mean(decodeneg(1:length(decodestruct.ConfMatAll_DAT)));
            DATSTRUCT.AllDecode_neg_mean = [DATSTRUCT.AllDecode_neg_mean;
                decode_neg_mean];
            
            % --------------------------- check neg decode distrubtions
            if plotdecodenegdistr ==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                lt_plot_histogram(decodeneg);
                line([decode decode], ylim, 'Color','r')
            end
            
            
            % ----------------- other params
            windowrelonset = decodestruct.window_relonset;
            AllWindows_RelOnsActual = [AllWindows_RelOnsActual; windowrelonset];
            
            windowrelOnsOff = timewindows.TimeWindows;
            AllWindows_RelOnsOffDesired = [AllWindows_RelOnsOffDesired; windowrelOnsOff];
        end
    end
end

AllSingleUnit = DATSTRUCT.AllSingleUnit;

%%
if plotON==1
    %% plot grand histogram
    for i=1:size(AllPdat,2)
        lt_figure; hold on;
        xlabel('log10(pval)');
        title(['prob of dat rel shuff distribution, ' analyfname]);
        lt_plot_histogram(log10(AllPdat(:,i) + 1/numshuffs));
        line([log10(0.05) log10(0.05)], ylim);
        
        numSig = sum(AllPdat(:)<0.05);
        numTot = length(AllPdat(:));
        
        lt_plot_annotation(1, [num2str(numSig) '/' num2str(numTot) ' (' num2str(numSig/numTot) ') sig (p<0.05)'])
    end
    
    %% pval correlate with zscore
    lt_figure; hold on;
    title([analyfname]);
    xlabel('zscore vs. null');
    ylabel('pval');
    plot(AllDecode_z, AllPdat, '.k');
    
    %%
    Numbirds = max(AllBirdNum);
    Numneur = max(AllNeurNum);
    Numbranch = max(AllBranchNum);
    BrainRegions = unique(AllBrainRegion)';
    
    %% === plot colrs for brain region
    plotcols_brainreg = lt_make_plot_colors(length(BrainRegions), 0, 0);
    
    %% ======== for each brain region/bird, plot distribtion of p-values
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ################### PVAL
    for bregion = BrainRegions
        plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        
        for i=1:Numbirds
            
            % ========================= COLLECT ALL NEURONS FOR THIS BIRD
            Yvals = {};
            XNeurs = [];
            birdname = SummaryStruct.birds(i).birdname;
            for ii=1:Numneur
                
                inds = AllBirdNum==i & AllNeurNum==ii & strcmp(AllBrainRegion, bregion);
                
                if ~any(inds)
                    continue
                end
                
                % ==== extract all p-vals for this bird/neuron (i.e. across
                % branches)
                Yvals = [Yvals log10(AllPdat(inds) + 1/numshuffs)];
                XNeurs = [XNeurs ii];
            end
            
            if isempty(Yvals)
                continue
            end
            
            % =========== plot for this bird
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(birdname);
            
            yvalsall = cell2mat(Yvals');
            lt_plot_histogram(yvalsall, '', 1, 0, [], 1, plotcol);
            line([log10(0.05) log10(0.05)], ylim, 'Color', 'k');
            
            % -------- proportion of cases overall significant
            numSig = sum(yvalsall(:)<log10(0.05));
            numTot = length(yvalsall(:));
            
            lt_plot_annotation(1, [num2str(numSig) '/' num2str(numTot) ' (' num2str(numSig/numTot) ') sig (p<0.05)'], 'r')
            
            % ========== binomial test for significance of fraction
            X = 0:numTot;
            P = 0.05;
            y = binopdf(X, numTot, P);
            
            pval = sum(y(X>=numSig));
            lt_plot_pvalue(pval, 'binomial test',1);
            
            if (0)
                lt_figure; hold on;
                plot(X, y, 'ok-');
                
                lt_plot_pvalue(pval, '',1);
                line([actualN actualN], ylim);
            end
        end
    end
    linkaxes(hsplots, 'xy');
    
    %% ======== [DISTRIBUTIONPLOT] for each brain region/bird, plot distribtion of p-values
    
    % ################### PVAL
    for bregion = BrainRegions
        plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        lt_figure; hold on;
        title(bregion);
        
        Ybybird = {};
        Birdnames = {};
        for i=1:Numbirds
            
            % ========================= COLLECT ALL NEURONS FOR THIS BIRD
            inds = AllBirdNum==i & strcmp(AllBrainRegion, bregion);
            
            yvalsall = AllPdat(inds);
            if isempty(yvalsall)
                continue
            end
            yvalsall = log10(AllPdat(inds) + 1/numshuffs);
            birdname = SummaryStruct.birds(i).birdname;
            
            %
            %
            %             Yvals = {};
            %             birdname = SummaryStruct.birds(i).birdname;
            %             for ii=1:Numneur
            %
            %                 inds = AllBirdNum==i & AllNeurNum==ii & strcmp(AllBrainRegion, bregion);
            %
            %                 if ~any(inds)
            %                     continue
            %                 end
            %
            %                 % ==== extract all p-vals for this bird/neuron (i.e. across
            %                 % branches)
            %                 Yvals = [Yvals log10(AllPdat(inds)+0.0001)];
            %             end
            %
            %             if isempty(Yvals)
            %                 continue
            %             end
            %
            % -------- proportion of cases overall significant
            
            % ============== OUTPUT
            Birdnames = [Birdnames birdname];
            Ybybird = [Ybybird yvalsall];
        end
        
        % ------------------ PLOT PVALS
        lt_subplot(2,1,1); hold on;
        title('pvals');
        ylabel('log10(p)');
        
        lt_plot_MultDist(Ybybird, 1:length(Ybybird), 0, 'k');
        plot(xlim, [log10(0.05) log10(0.05)], 'Color', 'm', 'LineStyle', '--');
        set(gca, 'XTick', 1:length(Ybybird), 'XTickLabel', Birdnames);
        rotateXLabels(gca, 90);
       xlim([0.5 length(Ybybird)+0.5]);
        
        
        
        % --------------- PLOT FARCTIONS note if is significant
        lt_subplot(2,1,2); hold on;
        title('fraction sign');
        ylabel('fraction p<0.05');
        
        Yfrac = [];
        Ypval = [];
        for j = 1:length(Ybybird)
            
            yvalsall = Ybybird{j};
            numSig = sum(yvalsall(:)<log10(0.05));
            numTot = length(yvalsall(:));
            
            lt_plot_text(j-0.2, 0.8, ...
                [num2str(numSig) '/' num2str(numTot)], 'r', 8)
            
            % ========== binomial test for significance of fraction
            X = 0:numTot;
            P = 0.05;
            y = binopdf(X, numTot, P);
            
            pval = sum(y(X>=numSig));
            if pval<0.1
            lt_plot_text(j-0.2, 0.7, '*' , 'r', 10)
            end
            %             lt_plot_pvalue(pval, 'binomial test',1);
            
            if (0)
                lt_figure; hold on;
                plot(X, y, 'ok-');
                
                lt_plot_pvalue(pval, '',1);
                line([actualN actualN], ylim);
            end
            
            % ==================
            Yfrac = [Yfrac numSig/numTot];
            Ypval = [Ypval pval];
            
        end
        
        lt_plot_bar(1:length(Ybybird), Yfrac, {'Color', 'k'});
               set(gca, 'XTick', 1:length(Ybybird), 'XTickLabel', Birdnames);
        rotateXLabels(gca, 90);
               xlim([0.5 length(Ybybird)+0.5]);
        line(xlim, [0.05 0.05]);
    end
    
    
    %% ########################  PLOT BY NEURON NUM [PVAL]
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ################### PVAL
    for bregion = BrainRegions
        plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        
        for i=1:Numbirds
            
            % ========================= COLLECT ALL NEURONS FOR THIS BIRD
            Yvals = {};
            XNeurs = [];
            birdname = SummaryStruct.birds(i).birdname;
            for ii=1:Numneur
                
                inds = AllBirdNum==i & AllNeurNum==ii & strcmp(AllBrainRegion, bregion);
                
                if ~any(inds)
                    continue
                end
                
                % ==== extract all p-vals for this bird/neuron (i.e. across
                % branches)
                Yvals = [Yvals log10(AllPdat(inds) + 1/numshuffs)];
                XNeurs = [XNeurs ii];
            end
            
            if isempty(Yvals)
                continue
            end
            
            % =========== plot for this bird
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(birdname);
            xlabel('neur num');
            ylabel('log10(prob)');
            lt_plot_MultDist(Yvals, XNeurs, 1, plotcol, 1, 0);
            line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
            ylim([-5 1]);
            
            % -------- proportion of cases overall significant
            yvalsall = cell2mat(Yvals');
            
            numSig = sum(yvalsall(:)<log10(0.05));
            numTot = length(yvalsall(:));
            
            lt_plot_annotation(1, [num2str(numSig) '/' num2str(numTot) ' (' num2str(numSig/numTot) ') sig (p<0.05)'], 'r')
        end
    end
    linkaxes(hsplots, 'xy');
    
    
    %% ######################## PLOT BY NEURON NUM [ZSCORE]
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    for bregion = BrainRegions
        plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        
        for i=1:Numbirds
            
            % ========================= COLLECT ALL NEURONS FOR THIS BIRD
            Yvals = {};
            XNeurs = [];
            birdname = SummaryStruct.birds(i).birdname;
            for ii=1:Numneur
                
                inds = AllBirdNum==i & AllNeurNum==ii & strcmp(AllBrainRegion, bregion);
                
                if ~any(inds)
                    continue
                end
                
                % ==== extract all p-vals for this bird/neuron (i.e. across
                % branches)
                Yvals = [Yvals AllDecode_z(inds)];
                XNeurs = [XNeurs ii];
            end
            
            if isempty(Yvals)
                continue
            end
            
            % =========== plot for this bird
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(birdname);
            xlabel('neur num');
            ylabel('decode (zscore)');
            lt_plot_MultDist(Yvals, XNeurs, 1, plotcol, 1, 0);
            %     line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
            %     ylim([-5 1]);
            lt_plot_zeroline;
            
        end
        
    end
    linkaxes(hsplots, 'xy');
    
    
    
    %% ##################### separated by branch num [PVAL]
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ############################# PVAL
    
    for bregion = BrainRegions
        plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        for i=1:Numbirds
            
            Yvals = {};
            Xbranch = [];
            Nnums = {}; % to collect neurons
            birdname = SummaryStruct.birds(i).birdname;
            BranchnameAll = {};
            IsSU = {};
            
            for ii=1:Numbranch
                
                inds = AllBirdNum==i & AllBranchNum==ii & strcmp(AllBrainRegion, bregion);
                
                if ~any(inds)
                    continue
                end
                
                % ==== extract all p-vals for this bird/branch (i.e. across
                % branches)
                Yvals = [Yvals log10(AllPdat(inds) + 1/numshuffs)];
                Xbranch = [Xbranch ii];
                
                % --- what is the name of this branch?
                tmp = find(inds);
                %         branchname = CLASSES.birds(i).neurons(AllNeurNum(tmp(1))).branchnum(ii).regexprstr;
                branchname = DATSTRUCT_BYBRANCH.bird(i).branchID(ii).regexpstr;
                BranchnameAll = [BranchnameAll branchname];
                
                % ------------ plot text of neurons for each datapoint
                Nnums = [Nnums AllNeurNum(inds)];
                
                IsSU = [IsSU AllSingleUnit(inds)];
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
                        suthis = IsSU{j}(jj);
                        neur = Nnums{j}(jj);
                        x = Xbranch(j);
                        y = Yvals{j}(jj);
                        if suthis==1
                            lt_plot_text(x, y, [num2str(neur)], 'b')
                        else
                            lt_plot_text(x, y, [num2str(neur)], [0.6 0.6 0.9])
                        end
                    end
                end
            end
            
            % -------- plot data
            lt_plot_MultDist(Yvals, Xbranch, 1, plotcol, 1, 0);
            line(xlim, [log10(0.05) log10(0.05)], 'Color','r');
            ylim([-5 1]);
            
            set(gca, 'XTickLabel', BranchnameAll);
            rotateXLabels(gca, 45);
            
            
            % -------- proportion of cases overall significant
            yvalsall = cell2mat(Yvals');
            
            numSig = sum(yvalsall(:)<log10(0.05));
            numTot = length(yvalsall(:));
            
            lt_plot_annotation(1, [num2str(numSig) '/' num2str(numTot) ' (' num2str(numSig/numTot) ') sig (p<0.05)'], 'r')
            
            
        end
    end
    
    linkaxes(hsplots, 'xy');
    
    
    %% ##################### separated by branch num [ZSCORE]
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    % ############################# PVAL
    for bregion = BrainRegions
        plotcol = plotcols_brainreg{strcmp(BrainRegions, bregion)};
        for i=1:Numbirds
            
            Yvals = {};
            Xbranch = [];
            Nnums = {}; % to collect neurons
            birdname = SummaryStruct.birds(i).birdname;
            BranchnameAll = {};
            IsSU = {};
            
            for ii=1:Numbranch
                
                inds = AllBirdNum==i & AllBranchNum==ii & strcmp(AllBrainRegion, bregion);
                
                if ~any(inds)
                    continue
                end
                
                % ==== extract all p-vals for this bird/branch (i.e. across
                % branches)
                Yvals = [Yvals AllDecode_z(inds)];
                Xbranch = [Xbranch ii];
                
                % --- what is the name of this branch?
                tmp = find(inds);
                %         branchname = CLASSES.birds(i).neurons(AllNeurNum(tmp(1))).branchnum(ii).regexprstr;
                branchname = DATSTRUCT_BYBRANCH.bird(i).branchID(ii).regexpstr;
                BranchnameAll = [BranchnameAll branchname];
                
                % -- is SU?
                IsSU = [IsSU AllSingleUnit(inds)];
                
                % ------------ plot text of neurons for each datapoint
                Nnums = [Nnums AllNeurNum(inds)];
            end
            
            if isempty(Yvals)
                continue
            end
            
            
            % =========== plot for this bird
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title(birdname);
            xlabel('branch ID');
            ylabel('Z-SCORE');
            
            
            % ------------- plot text of neuron num next to point
            if plottext==1
                for j=1:length(Nnums)
                    for jj = 1:length(Nnums{j})
                        suthis = IsSU{j}(jj);
                        
                        neur = Nnums{j}(jj);
                        x = Xbranch(j);
                        y = Yvals{j}(jj);
                        if suthis==1
                            lt_plot_text(x, y, [num2str(neur)], 'b')
                        else
                            lt_plot_text(x, y, [num2str(neur)], [0.6 0.6 0.9])
                        end
                    end
                end
            end
            
            % -------- plot data
            lt_plot_MultDist(Yvals, Xbranch, 1, plotcol, 1, 0);
            lt_plot_zeroline;
            
            set(gca, 'XTickLabel', BranchnameAll);
            rotateXLabels(gca, 45);
            
        end
    end
    
    linkaxes(hsplots, 'xy');
    
    %% -- separate plot for each bird, matrix of neuron x branch
    Numbirds = max(AllBirdNum);
    Numneur = max(AllNeurNum);
    Numbranch = max(AllBranchNum);
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    
    
    for i=1:Numbirds
        
        for ii=1:Numbranch
            
            
        end
        
    end
    
    
end

%% =========== for output

DATSTRUCT.AllPdat = AllPdat;
DATSTRUCT.AllBirdNum = AllBirdNum;
DATSTRUCT.AllNeurNum = AllNeurNum;
DATSTRUCT.AllBranchNum = AllBranchNum;
DATSTRUCT.AllDecode = AllDecode;
DATSTRUCT.AllDecode_z = AllDecode_z;
DATSTRUCT.AllBrainRegion = AllBrainRegion;
DATSTRUCT.AllWindows_RelOnsActual = AllWindows_RelOnsActual;
DATSTRUCT.AllWindows_RelOnsOffDesired = AllWindows_RelOnsOffDesired;







