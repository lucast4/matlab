function lt_neural_v2_CTXT_BRANCH_PlotByBranchID(ALLBRANCH, BrainRegions, ...
    BirdToPlot, useDprime, ExptToPlot, analyfname, sortbypremotorpval, ...
    plotAllZscoreFR)
apos =1;

%% lt 11/29/17 - plots, sepaated by unique branch points (i.e. regexpstr)

if strcmp(ALLBRANCH.alignpos(apos).ParamsFirstIter.Extract.strtype, 'xaa')
    wh44tempfix = 1;
else
    wh44tempfix = 0;
end


%%

if sortbypremotorpval==1
    assert(~isempty(analyfname), 'need analyfname to extract premotor decode pvalue');
end


%% allowed to specify filename and will open automatically + will also
% extract data about decode (using premotor window)

% --- if filename is not empty, then load it.
if ~isempty(analyfname)
    disp('CLEARING ALLBRANCH AND REOPENING (using analyfname)');
    clear ALLBRANCH
    savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';
    
    load([savedir '/ALLBRANCHv2_' analyfname '.mat']);
    load([savedir '/SUMMARYv2_' analyfname '.mat']);
    try
        load([savedir '/PARAMSv2_' analyfname '.mat']);
    catch err
    end
    
end


%%
strtype = ALLBRANCH.alignpos(apos).ParamsFirstIter.Extract.strtype;

if strcmp(ALLBRANCH.alignpos(apos).ParamsFirstIter.Extract.strtype, 'xaa')
    % datwindow_xcov = [-0.1 0.05]; % data to use, rel onset
    datwindow_xcov = [-0.12 0.07]; % data to use, rel onset
elseif strcmp(ALLBRANCH.alignpos(apos).ParamsFirstIter.Extract.strtype, 'aax') % then is divergent
    datwindow_xcov = [-0.05 0.14]; % data to use, rel onset
end



%% ======= Filter data (e.g. remove noise, poor labels, etc)

Params.LocationsToKeep = BrainRegions;
Params.birdstoexclude = {}; % performs 1st
Params.birdstokeep = BirdToPlot; % performs 2nd
Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
Params.durThreshOmega.syl = []; % omega2 (will only keep if lower) [leave empty to ignore]
Params.durThreshOmega.gappre= [];
Params.durThreshOmega.gappost= [];
Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
% class within branch) is longer than this (sec)
Params.RemoveHandCoded =1 ; % see below
Params.expttokeep = ExptToPlot;

if strcmp(ALLBRANCH.alignpos(apos).ParamsFirstIter.Extract.strtype, 'aax') % then is divergent
    Params.GapDurPostMax = 0.1; %
end
ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params);

%%
if useDprime==1
    % DPRIME STUFF
    Niter = 3; % shuffles
    Nmin = 3; % min sample size;
    % Nmin = ALLBRANCH.alignpos(1).ParamsFirstIter.minN; % min sample size;
    
    % DprimeNegVersion = 'Wohl'; % DONT USE THIS - is biased to be large,  because of lower sample size (this splits data in half)
    DprimeNegVersion = 'shuff'; % correct version, shuffles and resplits, maintaining sampel size.
    
    ALLBRANCH = lt_neural_v2_CTXT_AllBranchDprime(ALLBRANCH, Nmin, Niter, ...
        DprimeNegVersion);
end

%% ===== organize data by unique branches
if ~isempty(analyfname)
    DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, analyfname, 1, ...
        useDprime);
    
else
    DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, '', '', ...
        useDprime);
end
%% ===== FOR EACH BIRD PLOT EACH BRANCH, SEPARATING BY BRAIN REGION
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=max([2 length(locationstoplot)]);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        for loc = locationstoplot
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            % -------------- PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
            
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
            
            sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
            sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
            sylcontours_mean = mean(sylcontours_mean,1);
            sylcontours_x = mean(sylcontours_x,1);
            
            % ----- get premotor decode if exists
            predecode_p = [];
            if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
            end
            
            if ~isempty(ydecode)
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % -- lines
                line([0 0], ylim);
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                
                
                plot(xdecode, ydecode_neg, '-', 'Color', [0.6 0.6 0.6]);
                plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                plot(xdecode, ydecode_pos, '-b');
                plot(xdecode, mean(ydecode_pos,1), '-b', 'LineWidth',2);
                
                if sortbypremotorpval ==1
                    % then color differently based on whether premotor
                    % decode was significant
                    
                    % -- significant
                    if any(predecode_p<0.05)
                        plot(xdecode, ydecode(predecode_p<0.05, :), '-', 'Color', 'r');
                        plot(xdecode, mean(ydecode(predecode_p<0.05, :),1), '-r', 'LineWidth', 2)
                    end
                    % -- not significant
                    if any(predecode_p>=0.05)
                        plot(xdecode, ydecode(predecode_p>=0.05, :), '-', 'Color', 'm');
                    end
                    % mean of all
                    
                else
                    plot(xdecode, ydecode, '-', 'Color', 'r');
                    plot(xdecode, mean(ydecode,1), '-r', 'LineWidth', 2)
                end
                
                
                ymax = max(ydecode(:));
                plot(sylcontours_x, 0.25*sylcontours_mean+1, '-m');
            end
            
            % ---- note down neuron/channel number
            neuronOrig = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
            for j=1:length(neuronOrig)
                
                nn = neuronOrig(j);
                
                chan = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).channel;
                lt_plot_text(xdecode(end), ydecode(j,end), ['ch' num2str(chan)], 'b', 7);
                
            end
            
            % ===================== TEST WHETHER SIGNIFICANT DECODING
            % ACCURACY
            if (0)
                for tt = 1:length(xdecode)
                    % in each time bin do ttest comparing data to shuffle
                    
                    [h, p] = ttest(ydecode(:,tt), ydecode_neg(:,tt));
                    if p<0.05
                        plot(xdecode, 1, 'ob');
                    end
                end
            else
                
                % ---------------------- ttest of entire trajectory vs. neg
                % control
                [h, p] = ttest(mean(ydecode,2), mean(ydecode_neg,2));
                lt_plot_pvalue(p, 'mean decode vs. neg', 1);
            end
            
        end
        
    end
    axis tight
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
end

%% ===== FOR EACH BIRD PLOT EACH BRANCH, SEPARATING BY BRAIN REGION
% ======== ALL SUBTRACTING NEGATIVE CONTROL

numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=max([2 length(locationstoplot)]);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        for loc = locationstoplot
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            % -------------- PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
            ylabel('minus neg control');
            
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
            
            sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
            sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
            sylcontours_mean = mean(sylcontours_mean,1);
            sylcontours_x = mean(sylcontours_x,1);
            
            % =============== subtract negative control
            ydecode = ydecode - ydecode_neg;
            ydecode_pos = ydecode_pos - ydecode_neg;
            
            % ----- get premotor decode if exists
            predecode_p = [];
            if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
            end
            
            if ~isempty(ydecode)
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % -- lines
                line([0 0], ylim);
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                
                
%                 plot(xdecode, ydecode_neg, '-', 'Color', [0.6 0.6 0.6]);
%                 plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                plot(xdecode, ydecode_pos, '-b');
                plot(xdecode, mean(ydecode_pos,1), '-b', 'LineWidth',2);
                
                if sortbypremotorpval ==1
                    % then color differently based on whether premotor
                    % decode was significant
                    
                    % -- significant
                    if any(predecode_p<0.05)
                        plot(xdecode, ydecode(predecode_p<0.05, :), '-', 'Color', 'r');
                        plot(xdecode, mean(ydecode(predecode_p<0.05, :),1), '-r', 'LineWidth', 2)
                    end
                    % -- not significant
                    if any(predecode_p>=0.05)
                        plot(xdecode, ydecode(predecode_p>=0.05, :), '-', 'Color', 'm');
                    end
                    % mean of all
                    
                else
                    plot(xdecode, ydecode, '-', 'Color', 'r');
                    plot(xdecode, mean(ydecode,1), '-r', 'LineWidth', 2)
                end
                
                
                ymax = max(ydecode(:));
                plot(sylcontours_x, 0.25*sylcontours_mean+1, '-m');
            end
        end
        
    end
    axis tight
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
end


%% ============ FOR EACH NEURON, COLLECT ACROSS ALL BRANCHES


for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    
    figcount=1;
    subplotrows=4;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    AllXdecode = [];
    AllYdecode = [];
    AllYdecode_neg = [];
    AllYdecode_pos = [];
    AllNeurInds = [];
    AllBregions = {};
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        
        
        % ========== COLLECT
        xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
        ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode;
        ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg;
        ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos;
        neurInds = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron;
        bregions = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion;
        
        sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean;
        sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x;
        sylcontours_mean = mean(sylcontours_mean,1);
        sylcontours_x = mean(sylcontours_x,1);
        
        % ===========
        AllXdecode = [AllXdecode; xdecode];
        AllYdecode = [AllYdecode; ydecode];
        AllYdecode_neg = [AllYdecode_neg; ydecode_neg];
        AllYdecode_pos = [AllYdecode_pos; ydecode_pos];
        AllNeurInds = [AllNeurInds; neurInds];
        AllBregions = [AllBregions; bregions];
    end
    
    Xdecode = AllXdecode(1,:);
    
    % ================= GO THRU EACH NEURON AND PLOT
    for loc = locationstoplot
        
        neuronstoplot = AllNeurInds(strcmp(AllBregions, loc));
        neuronstoplot = unique(neuronstoplot); % list of neurons that have data for this bregion
        
        % ===========
        for n=1:length(neuronstoplot)
            nn = neuronstoplot(n);
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-n' num2str(nn) ' [' loc{1} ']']);
            
            
            indstmp = AllNeurInds == nn;
            
            x = Xdecode;
            ymat = AllYdecode(indstmp, :);
            yneg = AllYdecode_neg(indstmp, :);
            ypos = AllYdecode_pos(indstmp, :);
           
            plot(x, ymat, '-k');
            plot(x, mean(ymat,1), '-k', 'LineWidth', 2);
            
            plot(x, yneg, '-r');
            plot(x, mean(yneg,1), '-r', 'LineWidth', 2);
            
            plot(x, ypos, '-b');
            plot(x, mean(ypos,1), '-b', 'LineWidth', 2);
            
            line([0 0], ylim, 'LineStyle', '--', 'Color', 'y');
        end
    end
end


%% ===== FOR EACH BIRD AND BRANCH, OVERLAY REPRESENTATION OF ALL FR
if plotAllZscoreFR ==1
    numbirds = length(DATSTRUCT_BYBRANCH.bird);
    locationstoplot = BrainRegions;
    motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;
    for i=1:numbirds
        birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        figcount=1;
        subplotrows=2;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots =[];
        
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            % ========= plot for this branch
            hsplots = [];
            for loc = locationstoplot
                
                inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
                if ~any(inds)
                    continue
                end
                
                % -------------- PLOT
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
                
                
                xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
                ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
                %             ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
                %             ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
                
                %             sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
                %             sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
                %             sylcontours_mean = mean(sylcontours_mean,1);
                %             sylcontours_x = mean(sylcontours_x,1);
                
                % ----- get premotor decode if exists
                predecode_p = [];
                if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                    predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
                end
                
                
                % ============== FOR EACH NEURON, PLOT FR
                branchall = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_branch(inds);
                neuronall = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
                sylcontall = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
                sylcontx = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(1,:);
                sylcontx = sylcontx+motifpredur;
                
                if (1)
                    if sortbypremotorpval==1
                        branchall = branchall(predecode_p<0.05);
                        neuronall = neuronall(predecode_p<0.05);
                        sylcontall = sylcontall(predecode_p<0.05,:);
                        ydecode = ydecode(predecode_p<0.05, :);
                    end
                end
                
                % =========== for each neuron, collect FR across classes
                BranchNamesTmp = {};% -- collect branch names
                for n = 1:length(neuronall)
                    nn = neuronall(n);
                    branchthis = branchall(n);
                    
                    dattmp = ALLBRANCH.alignpos(apos).bird(i).branch(branchthis).neuron(nn);
                    BranchNamesTmp = [BranchNamesTmp dattmp.prms_regexpstr];
                    
                    % ==== collect original FR
                    numclasses = length(dattmp.FR.classnum);
                    FRmean_all = [];
                    for cc = 1:numclasses
                        
                        frmat = dattmp.FR.classnum(cc).FRsmooth_rate_CommonTrialDur;
                        frx = dattmp.FR.classnum(cc).FRsmooth_xbin_CommonTrialDur;
                        
                        frmean = mean(frmat,2);
                        
                        if ~isempty(FRmean_all)
                            if size(FRmean_all,2) ~= length(frmean)
                                veclength = min([size(FRmean_all,2), length(frmean)]);
                                FRmean_all = FRmean_all(:,1:veclength);
                                frmean = frmean(1:veclength);
                            end
                        end
                        
                        FRmean_all = [FRmean_all; frmean'];
                    end
                    
                    % ==== plot for this neuron
                    % -- first take global zscore
                    frmean = mean(FRmean_all(:));
                    frstd = std(FRmean_all(:));
                    Y = (FRmean_all - frmean)./frstd;
                    plot(frx, Y+3*(n-1), '-', 'Color', 0.2+0.6*[rand rand rand], 'LineWidth', 2);
                    
                end
                assert(length(unique(BranchNamesTmp))<2, 'diff branches ...');
                % === overlay syl contour
                plot(sylcontx, mean(sylcontall,1)+3*(length(neuronall)-1)+6, '--k')
                line([motifpredur motifpredur], ylim);
                
                % === overlay all decode traces
                if (0)
                    % - first take mean then take z
                    ydecode = mean(ydecode,1);
                    tmpmean = mean(ydecode(:));
                    tmpstd = std(ydecode(:), 0,1);
                    ydecode_z = (ydecode - tmpmean)./tmpstd;
                    
                    plot(xdecode+motifpredur, ydecode_z'+3*(length(neuronall)-1)+3, ':r', 'LineWidth', 1.5);
                    %             line(xlim, [3*(length(neuronall)-1)+3.5 3*(length(neuronall)-1)+3.5], 'LineStyle', ':');
                    %             line(xlim, [3*(length(neuronall)-1)+2.5 3*(length(neuronall)-1)+2.5], 'LineStyle', ':');
                end
                
            end
            linkaxes(hsplots, 'xy');
            if strcmp(strtype, 'xaa')
                xlim([0.03 0.22]);
            else
            end
        end
    end
    
end
%% ===== FOR EACH BIRD PLOT EACH BRANCH, SEPARATING BY BRAIN REGION
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=max([2 length(locationstoplot)]);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        for loc = locationstoplot
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            % -------------- PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
            
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
            
            sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
            sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
            sylcontours_mean = mean(sylcontours_mean,1);
            sylcontours_x = mean(sylcontours_x,1);
            
            % ----- get premotor decode if exists
            predecode_p = [];
            if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
            end
            
            if ~isempty(ydecode)
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % -- lines
                line([0 0], ylim);
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                
                plot(xdecode, ydecode_neg, '-', 'Color', [0.6 0.6 0.6]);
                plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                if sortbypremotorpval ==1
                    % then color differently based on whether premotor
                    % decode was significant
                    
                    % -- significant
                    if any(predecode_p<0.05)
                        plot(xdecode, ydecode(predecode_p<0.05, :), '-', 'Color', 'r');
                        plot(xdecode, mean(ydecode(predecode_p<0.05, :),1), '-r', 'LineWidth', 2)
                    end
                    % -- not significant
                    if any(predecode_p>=0.05)
                        plot(xdecode, ydecode(predecode_p>=0.05, :), '-', 'Color', 'm');
                    end
                    % mean of all
                    
                else
                    plot(xdecode, ydecode, '-', 'Color', 'r');
                    plot(xdecode, mean(ydecode,1), '-r', 'LineWidth', 2)
                end
                
                
                ymax = max(ydecode(:));
                plot(sylcontours_x, 0.25*sylcontours_mean+1, '-m');
            end
            
            % ---- note down neuron/channel number
            neuronOrig = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
            for j=1:length(neuronOrig)
                
                nn = neuronOrig(j);
                
                chan = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).channel;
                lt_plot_text(xdecode(end), ydecode(j,end), ['ch' num2str(chan)], 'b', 7);
                
            end
            
            % ===================== TEST WHETHER SIGNIFICANT DECODING
            % ACCURACY
            if (0)
                for tt = 1:length(xdecode)
                    % in each time bin do ttest comparing data to shuffle
                    
                    [h, p] = ttest(ydecode(:,tt), ydecode_neg(:,tt));
                    if p<0.05
                        plot(xdecode, 1, 'ob');
                    end
                end
            else
                
                % ---------------------- ttest of entire trajectory vs. neg
                % control
                [h, p] = ttest(mean(ydecode,2), mean(ydecode_neg,2));
                lt_plot_pvalue(p, 'mean decode vs. neg', 1);
            end
            
        end
        
    end
    axis tight
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
end







%% ===== FOR EACH BIRD AND BRANCH, OVERLAY REPRESENTATION OF ALL FR

if (0)
    % this plots STD of FR, but really noisy ...
    
    numbirds = length(DATSTRUCT_BYBRANCH.bird);
    locationstoplot = BrainRegions;
    
    for i=1:numbirds
        birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        figcount=1;
        subplotrows=4;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        hsplots =[];
        
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            % ========= plot for this branch
            for loc = locationstoplot
                
                inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
                if ~any(inds)
                    continue
                end
                
                % -------------- PLOT
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
                
                
                %             xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
                %             ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
                %             ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
                %             ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
                %
                %             sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
                %             sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
                %             sylcontours_mean = mean(sylcontours_mean,1);
                %             sylcontours_x = mean(sylcontours_x,1);
                
                % ----- get premotor decode if exists
                predecode_p = [];
                if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                    predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
                end
                
                
                % ============== FOR EACH NEURON, PLOT FR
                branchall = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_branch(inds);
                neuronall = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
                
                if sortbypremotorpval==1
                    branchall = branchall(predecode_p<0.05);
                    neuronall = neuronall(predecode_p<0.05);
                end
                
                % =========== for each neuron, collect FR across classes
                for n = 1:length(neuronall)
                    nn = neuronall(n);
                    branchthis = branchall(n);
                    
                    dattmp = ALLBRANCH.alignpos(apos).bird(i).branch(branchthis).neuron(nn);
                    
                    % ==== collect original FR
                    numclasses = length(dattmp.FR.classnum);
                    FRmean_all = [];
                    for cc = 1:numclasses
                        
                        frmat = dattmp.FR.classnum(cc).FRsmooth_rate_CommonTrialDur;
                        frx = dattmp.FR.classnum(cc).FRsmooth_xbin_CommonTrialDur;
                        
                        frmean = mean(frmat,2);
                        
                        FRmean_all = [FRmean_all; frmean'];
                    end
                    
                    % ================= PLOT FR
                    % --- z-transform (relative to global mean and sd)
                    frmean = mean(FRmean_all(:));
                    frstd = std(FRmean_all(:));
                    Y = (FRmean_all - frmean)./frstd;
                    if (0)
                        plot(frx, Y, '-', 'Color', 0.2+0.6*[rand rand rand])
                        
                    else
                        % --- running std
                        Ystd = std(Y,0,1);
                        plot(frx, Ystd, '-', 'Color', 0.2+0.6*[rand rand rand])
                    end
                end
                
                
            end
        end
    end
end

%% ===== FOR EACH BIRD AND BRANCH, OVERLAY MEAN OF DIFF BRAIN REGIONS

numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;
plotcols = lt_make_plot_colors(length(locationstoplot), 0,0);

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=max([2 length(locationstoplot)]);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title([birdname '-' thisbranch{1}]);
        
        % ========= plot for this branch
        for ll = 1:length(locationstoplot)
            loc = locationstoplot{ll};
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            % -------------- PLOT
            
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
            sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
            sylcontours_mean = mean(sylcontours_mean,1);
            sylcontours_x = mean(sylcontours_x,1);
            
            % ----- get premotor decode if exists
            predecode_p = [];
            if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
            end
            
            if ~isempty(ydecode)
                % -- lines
                line([0 0], ylim);
                
                % -- plot negative control
                if sortbypremotorpval ==1
                    if any(predecode_p<0.05)
                        plot(xdecode, mean(ydecode_neg(predecode_p<0.05,:), 1), '--', 'LineWidth', 1, 'Color', plotcols{ll});
                    end
                else
                    plot(xdecode, mean(ydecode_neg, 1), '--', 'LineWidth', 1, 'Color', plotcols{ll});
                end
                
                % --- -plot data
                if sortbypremotorpval ==1
                    % then color differently based on whether premotor
                    % decode was significant
                    
                    % -- significant
                    if any(predecode_p<0.05)
                        ymean = mean(ydecode(predecode_p<0.05, :),1);
                        ysem = lt_sem(ydecode(predecode_p<0.05, :));
                        
                        if length(ysem)>1
                            shadedErrorBar(xdecode, ymean, ysem, {'Color', plotcols{ll}}, 1);
                        end
                        plot(xdecode, ymean, '-', 'LineWidth', 2, ...
                            'Color', plotcols{ll});
                    end
                    
                else
                    ymean = mean(ydecode,1);
                    ysem = lt_sem(ydecode);
                    
                    if length(ysem)>1
                        shadedErrorBar(xdecode, ymean, ysem, {'Color', plotcols{ll}}, 1);
                    end
                    plot(xdecode, ymean, '-', 'LineWidth', 2, ...
                        'Color', plotcols{ll});
                    
                end
                
                ymax = max(ydecode(:));
                plot(sylcontours_x, 0.25*sylcontours_mean+0.8, '-m');
            end
            
            % ---- note down neuron/channel number
            neuronOrig = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
            for j=1:length(neuronOrig)
                
                nn = neuronOrig(j);
                
                chan = ALLBRANCH.SummaryStruct.birds(i).neurons(nn).channel;
                lt_plot_text(xdecode(end), ydecode(j,end), ['ch' num2str(chan)], 'b', 7);
                
            end
            
            % ===================== TEST WHETHER SIGNIFICANT DECODING
            % ACCURACY
            if (0)
                for tt = 1:length(xdecode)
                    % in each time bin do ttest comparing data to shuffle
                    
                    [h, p] = ttest(ydecode(:,tt), ydecode_neg(:,tt));
                    if p<0.05
                        plot(xdecode, 1, 'ob');
                    end
                end
            else
                
                % ---------------------- ttest of entire trajectory vs. neg
                % control
                [h, p] = ttest(mean(ydecode,2), mean(ydecode_neg,2));
                lt_plot_pvalue(p, 'mean decode vs. neg', 1);
            end
            
        end
        
    end
    axis tight
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
end

%% ===== FOR EACH BIRD AND BRANCH, OVERLAY MEAN OF DIFF BRAIN REGIONS
% PLOT USING ZSCORE

numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
plotcols = lt_make_plot_colors(length(locationstoplot), 0,0);

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=max([2 length(locationstoplot)]);
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        hsplots = [hsplots hsplot];
        title([birdname '-' thisbranch{1}]);
        
        % ========= plot for this branch
        for ll = 1:length(locationstoplot)
            loc = locationstoplot{ll};
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            % -------------- PLOT
            
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            
            indtmp = xdecode>=datwindow_xcov(1) & xdecode<=datwindow_xcov(2);
            xdecode = xdecode(indtmp);
            ydecode = ydecode(:, indtmp);
            
            
            % ----- get premotor decode if exists
            predecode_p = [];
            if isfield(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT, 'PREMOTORDECODE_pval')
                predecode_p = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
            end
            
            if ~isempty(ydecode)
                % -- lines
                line([0 0], ylim);
                
                %                 % -- plot negative control
                %                 if sortbypremotorpval ==1
                %                     if any(predecode_p<0.05)
                %                     plot(xdecode, mean(ydecode_neg(predecode_p<0.05,:), 1), '--', 'LineWidth', 1, 'Color', plotcols{ll});
                %                     end
                %                 else
                %                     plot(xdecode, mean(ydecode_neg, 1), '--', 'LineWidth', 1, 'Color', plotcols{ll});
                %                 end
                %
                % --- -plot data
                if sortbypremotorpval ==1
                    % then color differently based on whether premotor
                    % decode was significant
                    
                    % -- significant
                    if any(predecode_p<0.05)
                        ymean = mean(ydecode(predecode_p<0.05, :),1);
                        ymean = (ymean - mean(ymean))./std(ymean);
                        plot(xdecode, ymean, '-', 'LineWidth', 2, ...
                            'Color', plotcols{ll});
                    end
                else
                    ymean = mean(ydecode,1);
                    ymean = (ymean - mean(ymean))./std(ymean);
                    plot(xdecode, ymean, '-', 'LineWidth', 2, ...
                        'Color', plotcols{ll});
                    
                    
                end
            end
            
            
        end
        
    end
    axis tight
    if ~isempty(hsplots)
        linkaxes(hsplots, 'xy');
    end
    ylabel('zscore');
end


%% ======================= SUMMARY PLOTS FOR EACH REGION
figcount=1;
subplotrows=4;
subplotcols=max([2 length(locationstoplot)]);
fignums_alreadyused=[];
hfigs=[];
hsplots =[];



%% =========================== BY BIRD
numbirds = length(DATSTRUCT_BYBRANCH.bird);

for loc = locationstoplot
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BY BIRD]' loc{1}]);
    
    Yall = [];
    for i=1:numbirds
        birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        Ythisbird = [];
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            %  ================ GET brain region
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            % -------------- GET DAT
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            Ythisbird = [Ythisbird; ydecode];
            
        end
        
        % --- plot for this bird
        Ythisbirdmean = mean(Ythisbird,1);
        if ~isempty(Ythisbirdmean)
            plot(xdecode, Ythisbirdmean, '-k');
            % ==================== COLLECT ACROSS BIRDS
            Yall = [Yall; Ythisbirdmean];
            disp(birdname);
        end
        
    end
    
    % -- plot mean and sem
    Yallmean = mean(Yall,1);
    Yallsem = lt_sem(Yall);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, Yallmean, Yallsem, {'Color', 'r'},1)
    end
end



%% =========================== BY BRANCH
numbirds = length(DATSTRUCT_BYBRANCH.bird);

for loc = locationstoplot
    disp(loc)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BRANCHES]' loc{1}]);
    
    Yall = [];
    Yallneg = [];
    for i=1:numbirds
        %         birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            %  ================ GET brain region
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            % -------------- GET DAT
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            %% take mean across all neurons for this branch
            ydecode = mean(ydecode,1);
            ydecode_neg = mean(ydecode_neg, 1);
            
            %%
            
            if ~isempty(ydecode)
                plot(xdecode, ydecode, '-', 'Color', 'k');
                
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                %                 plot(xdecode, ydecode_neg, '-', 'Color', 'k');
                %                 plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                % -- lines
                line([0 0], ylim);
                
            end
            
            
            % ==================== COLLECT
            Yall = [Yall; ydecode];
            Yallneg = [Yallneg; ydecode_neg];
        end
    end
    % -- plot mean and sem
    Yallmean = mean(Yall,1);
    Yallsem = lt_sem(Yall);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, Yallmean, Yallsem, {'Color', 'r'},1)
    end
    
    % -- neg
    YallmeanNEG = mean(Yallneg,1);
    YallsemNEG = lt_sem(Yallneg);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, YallmeanNEG, YallsemNEG, {'Color', 'k'},1)
    end
    
end


%% =========================== BY BRANCH/NEURON
numbirds = length(DATSTRUCT_BYBRANCH.bird);


Ymean_loc =[];
Ysem_loc = [];
for loc = locationstoplot
    disp(loc)
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    hsplots = [hsplots hsplot];
    title(['[BRANCHES/NEURON]' loc{1}]);
    
    Yall = [];
    Yall_neg = [];
    for i=1:numbirds
        %         birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
            
            %  ================ GET brain region
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            % -------------- GET DAT
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            
            %%
            
            if ~isempty(ydecode)
                plot(xdecode, ydecode, '-', 'Color', [0.6 0.6 0.6]);
                
                %                 % -- plot mean
                %                 if size(ydecode,1)>1
                %                     ymean = mean(ydecode,1);
                %                     ysem = lt_sem(ydecode);
                %                     shadedErrorBar(xdecode, ymean, ysem, {'Color', 'r'}, 1);
                %                 end
                
                % --- plot mean + SD of neg control
                %                 ymean_neg = mean(ydecode_neg, 1);
                %                 ystd_neg = std(ydecode_neg, 0,1);
                %                 shadedErrorBar(xdecode, ymean_neg, ystd_neg, {'Color', 'k'}, 1);
                %                 plot(xdecode, ydecode_neg, '-', 'Color', 'k');
                %                 plot(xdecode, mean(ydecode_neg, 1), '-k', 'LineWidth', 2);
                
                % -- lines
                line([0 0], ylim);
                
            end
            
            
            % ==================== COLLECT
            Yall = [Yall; ydecode];
            Yall_neg = [Yall_neg; ydecode_neg];
        end
    end
    % -- plot mean and sem
    Yallmean = mean(Yall,1);
    Yallsem = lt_sem(Yall);
    if size(Yall,1)>1
        shadedErrorBar(xdecode, Yallmean, Yallsem, {'Color', 'r'},1)
    end
    
    % --- put neg control
    Yallneg_mean = mean(Yall_neg,1);
    Yallneg_std = std(Yall_neg,0,1);
    Yallneg_sem = lt_sem(Yall_neg);
    shadedErrorBar(xdecode, Yallneg_mean, Yallneg_sem, {'Color', 'm'},1)
    
    
    % ================== COLLECT FOR THIS LOC
    Ymean_loc = [Ymean_loc; Yallmean];
    Ysem_loc = [Ysem_loc; Yallsem];
end

% ==================== PLOT ALL LOC COMBINED
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title(['[COMBINED]']);

for i=1:size(Ymean_loc,1)
    shadedErrorBar(xdecode, Ymean_loc(i,:), Ysem_loc(i,:), {'Color', 'r'}, 1);
end


linkaxes(hsplots, 'xy');

%%
if length(BrainRegions)==2
    %% ==== cross correlation between each pair of brain regions
    if (0) % OLD VERSION !!!!!!!!!!!!!!!!!
        %         datwindow = [-0.1 0.05]; % data to use, rel onset
        %         windowmax = 0.05; % in sec
        datwindow = []; % data to use, rel onset
        windowmax = 0.1; % in sec
        
        YallMat = [];
        binsizeholder = [];
        figure; hold on;
        for i=1:numbirds
            numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
            
            for bb=1:numbranches
                
                thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
                Yall = cell(1, length(locationstoplot));
                
                for k = 1:length(locationstoplot)
                    loc = locationstoplot{k};
                    
                    %  ================ GET brain region
                    inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc); % neurons
                    if ~any(inds)
                        continue
                    end
                    % -------------- GET DAT
                    xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
                    binsize = xdecode(2)-xdecode(1);
                    binsizeholder = [binsizeholder binsize];
                    ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
                    ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
                    
                    if ~isempty(datwindow)
                        indstmp = xdecode>=datwindow(1) & xdecode<=datwindow(2);
                        
                        ydecode = ydecode(:,indstmp);
                        ydecode_neg = ydecode_neg(:, indstmp);
                        
                    end
                    % ==================== COLLECT
                    Yall{k} = ydecode;
                    
                end
                
                % ==== CALCULATE ALL CROSS CORRELATIONS (between all pairs of
                % neurons)
                nneur1 = size(Yall{1},1);
                nneur2 = size(Yall{2},1);
                for k = 1:nneur1
                    for kk = 1:nneur2
                        
                        dat1 = Yall{1}(k,:);
                        dat2 = Yall{2}(kk,:);
                        
                        [cc, lags] = xcov(dat1, dat2, ceil(windowmax/binsize), 'coeff');
                        %                         [cc, lags] = xcorr(dat1, dat2, ceil(windowmax/binsize));
                        plot(lags*binsize, cc);
                        YallMat = [YallMat; cc];
                    end
                end
                
                
            end
        end
        
        assert(length(unique(binsizeholder))==1, 'asdfsd');
        
        plot(lags*binsize, mean(YallMat, 1), 'k');
        xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
        
    end
    
    %% ==== cross correlation between each pair of brain regions
    numbirds = length(DATSTRUCT_BYBRANCH.bird);
    % #################################### 1) COLLECT ALL DATA INTO ARRAYS
    YallMat = [];
    XallMat = [];
    BranchID = [];
    LocationAll = [];
    BirdID = [];
    
    binsizeholder = [];
    
    for i=1:numbirds
        numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
        
        for bb=1:numbranches
            
            for k = 1:length(locationstoplot)
                loc = locationstoplot{k};
                
                %  ================ GET neurons
                inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc); % neurons
                if ~any(inds)
                    continue
                end
                
                % --------------- COLLECT DAT
                xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(inds,:);
                ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
                ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
                
                % --------------- IF DESIRED, ONLY KEEP THOSE WITH
                % SIGNIFICANT PREMOTOR DECODE
                if sortbypremotorpval==1
                    decodepval = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.PREMOTORDECODE_pval(inds);
                    xdecode = xdecode(decodepval < 0.05, :);
                    ydecode = ydecode(decodepval < 0.05, :);
                    ydecode_neg = ydecode_neg(decodepval < 0.05, :);
                    
                    % ---- skip if this throws out all data
                    if isempty(ydecode)
                        continue
                    end
                end
                
                % --------- withold to premotor window here
                if ~isempty(datwindow_xcov)
                    indstmp = xdecode(1,:)>=datwindow_xcov(1) & ...
                        xdecode(1,:)<=datwindow_xcov(2);
                    xdecode = xdecode(:,indstmp);
                    ydecode = ydecode(:,indstmp);
                    ydecode_neg = ydecode_neg(:, indstmp);
                end
                
                % ---- collect all binsizes - will make sure all are equal
                binsize = xdecode(1,2)-xdecode(1,1);
                binsizeholder = [binsizeholder binsize];
                
                % ==================== COLLECT
                YallMat = [YallMat; ydecode];
                XallMat = [XallMat; xdecode];
                BranchID = [BranchID; bb*ones(size(ydecode,1),1)];
                LocationAll = [LocationAll; k*ones(size(ydecode,1),1)];
                BirdID = [BirdID; i*ones(size(ydecode,1),1)];
                
            end
            
            %         % ==== CALCULATE ALL CROSS CORRELATIONS (between all pairs of
            %         % neurons)
            %         nneur1 = size(Yall{1},1);
            %         nneur2 = size(Yall{2},1);
            %         for k = 1:nneur1
            %             for kk = 1:nneur2
            %
            %                 dat1 = Yall{1}(k,:);
            %                 dat2 = Yall{2}(kk,:);
            %
            %                 [cc, lags] = xcov(dat1, dat2, ceil(windowmax/binsize), 'coeff');
            %                 plot(lags*binsize, cc);
            %                 YallMat = [YallMat; cc];
            %             end
            %         end
            %
            %
        end
    end
    
    assert(length(unique(binsizeholder))==1, 'asdfsd');
    
    
    
    %% ============ LIMIT TO BRANCHES THAT ARE NOT DEPENDENT MOTIFS
    if wh44tempfix==1
        if (1)
            birdnumthis = find(strcmp({ALLBRANCH.SummaryStruct.birds.birdname}, 'wh44wh39'));
            
            indstoremove = ismember(BranchID, [1 2 3]) & ...
                BirdID == birdnumthis;
            
            BirdID(indstoremove) = [];
            BranchID(indstoremove) = [];
            LocationAll(indstoremove) = [];
            YallMat(indstoremove, :) = [];
            XallMat(indstoremove, :) = [];
            
        else
            % --
            indstokeep = ismember(BranchID, [4 5 6 7 8]);
            
            BirdID = BirdID(indstokeep);
            BranchID = BranchID(indstokeep);
            LocationAll = LocationAll(indstokeep);
            YallMat = YallMat(indstokeep, :);
            XallMat = XallMat(indstokeep, :);
        end
    end
    %% #####################################
    
    %% ==== do separately for each bird
    windowmax = 0.075; % sec
    doCConMeans = 0; % then uses mean for brain region and branch
    for j=1:max(BirdID)
        [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID(BirdID==j), ...
            BranchID(BirdID==j), LocationAll(BirdID==j), YallMat(BirdID==j, :), ...
            XallMat(BirdID==j, :), locationstoplot, 0, doCConMeans);
    end
    
    %% ============== COMBINE ACROSS BIRDS
    windowmax = 0.075; % sec
    doCConMeans = 0; % then uses mean for brain region and branch
    [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID, BranchID, LocationAll, ...
        YallMat, XallMat, locationstoplot, 0, doCConMeans);
    
    
    %%
    
    %% CALCULATE GRAND MEAN TO GET PEAK TIMING, AND BOOTSTRAP TO GET VARIABILITY
    
    % -------------- DAT
    ccmean_dat = mean(CCall);
    
    % -------------- SHUFF
    nshuff = 1000;
    N = size(CCall,1);
    ccmean_shuff = [];
    ccmean_shuff_interp = [];
    for nn=1:nshuff
        
        indtmp = randi(N, 1, N);
        
        ccmean_shuff = [ccmean_shuff; mean(CCall(indtmp,:),1)];
        
        % ---  smooth
        ccmean_shuff_interp = [ccmean_shuff_interp; ...
            smooth(mean(CCall(indtmp,:),1), 9)'];
        
    end
    
    
    % ======================== PLOT (NOT INTERPOLATED)
    lt_figure; hold on;
    lt_subplot(2,2,1); hold on;
    title('dat/shuff, no interp');
    plot(LagsAll(1,:).*binsize, ccmean_shuff, '-', 'Color', [0.7 0.7 0.7]);
    plot(LagsAll(1,:).*binsize, ccmean_dat, '-k', 'LineWidth', 3);
    
    % ---------------------
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    % ============ same, but with interpolation
    lt_subplot(2,2,2); hold on;
    title('dat/shuff, with interp');
    
    % --- shuff
    %
    %    ccmean_shuff_interp = interp1(1:length(ccmean_dat), ccmean_shuff', 1:length(ccmean_dat), ...
    %        'spline');
    plot(LagsAll(1,:).*binsize, ccmean_shuff_interp, '-', 'Color', [0.7 0.7 0.7]);
    
    % --- data
    %    xfit = LagsAll(1,:);
    
    %    xfit = 1:0.25:length(ccmean_dat);
    %    ccmean_dat_interp = interp1(1:length(ccmean_dat), ccmean_dat', xfit);
    ccmean_dat_interp = smooth(ccmean_dat, 9);
    %    x = 1:length(ccmean_dat);
    %    ccmean_dat_interp = fit(x', ccmean_dat', 'smoothingspline');
    plot(LagsAll(1,:).*binsize, ccmean_dat_interp, '-k', 'LineWidth', 3);
    
    % ---------------------
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    
    % ============================ PLOT DISTRIBUTION OF TIMINGS AND PEAKS
    % ----- DATA
    [pks, locs] = findpeaks(ccmean_dat_interp, 'SortStr', 'descend');
    pks_dat = pks(1);
    locs_dat = locs(1);
    
    % ---- SHUFF
    pks_shuff = [];
    locs_shuff = [];
    for j=1:size(ccmean_shuff_interp, 1)
        
        [pks, locs] = findpeaks(ccmean_shuff_interp(j, :), 'SortStr', 'descend');
        
        pks_shuff = [pks_shuff pks(1)];
        locs_shuff = [locs_shuff locs(1)];
        
    end
    
    lt_subplot(2,2,3); hold on;
    title('loc and size of peaks (dat, shuff)');
    ylabel(['nshuff=' num2str(nshuff)]);
    % -- SHUFF
    plot(LagsAll(1, locs_shuff).*binsize, pks_shuff, 'x', 'Color', [0.7 0.7 0.7]);
    % --- DAT
    plot(LagsAll(1,locs_dat).*binsize, pks_dat, 'ok', 'LineWidth', 3);
    
    xlim([LagsAll(1,1).*binsize LagsAll(1,end).*binsize]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ========================== CALCULATE MEAN AND STD OF SHUFF
    %    pkmean_shuff = mean(pks_shuff);
    %    pkstd_shuff = std(pks_shuff);
    %    locmean_shuff = mean(locs_shuff);
    %    locstd_shuff = std(locs_shuff);
    %
    %    lt_subplot(2,2,4); hold on;
    %    xlim([LagsAll(1,1).*binsize LagsAll(1,end).*binsize]);
    %    lt_plot_zeroline;
    %    lt_plot_zeroline_vert;
    %
    %    plot(LagsAll(1,locmean_shuff)
    %
    
    %% =================== SAVE LOCATION OF PEAK IN GRAND MEAN (PER BIRD)
    
    CCmean_dat = nan(size(CCall,2), numbirds);
    TimeOfPeak_dat = nan(1, numbirds);
    
    for i=1:numbirds
        
        inds = PairedBirdID==i;
        
        if ~any(inds)
            continue
        end
        
        % ------ collect mean CC
        ccmean = mean(CCall(inds,:),1);
        CCmean_dat(:, i) = ccmean;
        
        
        % ----------- find main peak of mean xcov
        [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
        timeofpeak = LagsAll(1,loctmp)*binsize;
        TimeOfPeak_dat(i) = timeofpeak;
    end
    
    
    
    %% ########################### NEGATIVE CONTROL - SHUFFLE LOCATION LABEL
    % ----
    nshuff = 1000;
    TimeOfPeak_Shuff = nan(nshuff, numbirds);
    CCmean_Shuff = nan(size(CCall,2), nshuff, numbirds);
    
    for nn=1:nshuff
        disp(['shuff ' num2str(nn)]);
        LocationAll_RAND = nan(size(LocationAll));
        
        % ------------ RANDOMIZE NEURON LOCATION (MAINTAINING BIRD/BRANCH ID)
        for i=1:numbirds
            numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
            for bb = 1:numbranches
                %                 disp([num2str(i) '-' num2str(bb)]);
                inds = find(BirdID==i & BranchID==bb);
                
                %                 if isempty(inds)
                % %                     LocationAll_RAND(inds) = 97; % put some crazy number. don't want to be nan, but don't want it to affect actual data
                %                     continue
                %                 end
                
                indsperm = inds(randperm(length(inds)));
                
                %                 assert(all(isnan(LocationAll_RAND(inds))), 'asdf');
                LocationAll_RAND(inds) = LocationAll(indsperm);
            end
        end
        
        assert(~any(isnan(LocationAll_RAND)), 'asdfds');
        
        
        % =========================== CALC COV
        if nn>0
            suppressplots = 1;
        else
            suppressplots=0;
        end
        [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID, BranchID, LocationAll_RAND, ...
            YallMat, XallMat, locationstoplot, suppressplots, doCConMeans);
        
        % =========================== CALCUALTE TIME OF PEAK (and collect
        % mean timecourse
        for i=1:numbirds
            
            inds = PairedBirdID==i;
            
            if ~any(inds)
                continue
            end
            
            % ------ collect mean CC
            ccmean = mean(CCall(inds,:),1);
            CCmean_Shuff(:, nn, i) = ccmean;
            
            % ----------- find main peak of mean xcov
            [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
            timeofpeak = LagsAll(1,loctmp)*binsize;
            
            TimeOfPeak_Shuff(nn, i) = timeofpeak;
        end
    end
    
    %% ====================== COMPARE SHUFFLES TO DATA
    
    
    % =========================1 ) time of peak
    figcount=1;
    subplotrows=2;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    
    for i=1:numbirds
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
        xcenters = min(TimeOfPeak_Shuff(:,i)-0.005):0.005:(max(TimeOfPeak_Shuff(:,i))+0.005);
        lt_plot_histogram(TimeOfPeak_Shuff(:,i), xcenters, 1, 1, '', 1, 'k');
        line([TimeOfPeak_dat(i) TimeOfPeak_dat(i)], ylim, 'Color', 'r', 'LineWidth', 2);
        
        % --- p value
        p = (sum(abs(TimeOfPeak_Shuff(:,i)) > abs(TimeOfPeak_dat(i))) + 1)./(length(TimeOfPeak_Shuff(:,i)) + 1);
        lt_plot_pvalue(p, '2-tailed',1);
        xlim(binsize*[LagsAll(1,1) LagsAll(1,end)]);
    end
    
    % ======================= 2) overlay data mean CC on shuffles
    figcount=1;
    subplotrows=2;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    for i=1:numbirds
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
        % ---- overlay all shuffle dat
        ccmeanall = CCmean_Shuff(:,:,i);
        plot(LagsAll(1,:)*binsize, ccmeanall', ':', 'Color', [0.7 0.7 0.7]);
        plot(LagsAll(1,:)*binsize, CCmean_dat(:, i), 'Color', 'r');
        xlim(binsize*[LagsAll(1,1) LagsAll(1,end)]);
    end
    
    % ======================== 3) overlay mean on shuff mean + sd
    figcount=1;
    subplotrows=2;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    for i=1:numbirds
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
        % ---- overlay all shuffle dat
        ccmeanall = CCmean_Shuff(:,:,i);
        ccmeanmean = mean(ccmeanall, 2);
        ccmeanstd = std(ccmeanall, 0, 2);
        ccmean95 = prctile(ccmeanall', [2.5 97.5]);
        %         ccmean95 = prctile(ccmeanall', [2.5 97.5]);
        %         ccmean95 =  ccmean95 - repmat(ccmeanmean', 2, 1);
        plot(LagsAll(1,:)*binsize, ccmeanmean, '-k', 'LineWidth', 2);
        plot(LagsAll(1,:)*binsize, ccmean95(1,:)', '-k');
        plot(LagsAll(1,:)*binsize, ccmean95(2,:)', '-k');
        %
        %
        %         shadedErrorBar(LagsAll(1,:)*binsize, ccmeanmean, ccmeanstd, {'Color', 'k'}, 1);
        plot(LagsAll(1,:)*binsize, CCmean_dat(:, i), '-r', 'LineWidth', 3);
        xlim(binsize*[LagsAll(1,1) LagsAll(1,end)]);
    end
    
    % ============================ 4) SUM POSITIVE AND NEGATIVE TIME LAGS
    figcount=1;
    subplotrows=2;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    for i=1:numbirds
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(['bird ' num2str(i)]);
        
        % ---- overlay all shuffle dat
        ccmeanall = CCmean_Shuff(:,:,i);
        
        % --- negative sums
        indtmp = LagsAll(1,:)<0;
        yneg = mean(ccmeanall(indtmp, :), 1);
        
        % --- positive sums
        indtmp = LagsAll(1,:)>0;
        ypos = mean(ccmeanall(indtmp, :), 1);
        
        % --- take difference
        y_shuff = ypos - yneg;
        
        
        % === data
        ypos = mean(CCmean_dat(LagsAll(1,:)>0, i));
        yneg = mean(CCmean_dat(LagsAll(1,:)<0, i));
        y_dat = ypos - yneg;
        
        % ======= plot
        lt_plot_histogram(y_shuff);
        line([y_dat y_dat], ylim);
        
        p = (sum(abs(y_shuff) > abs(y_dat)) + 1)./(length(y_shuff) + 1);
        lt_plot_pvalue(p, '2-tailed',1);
    end
    
    
    
    
end
end




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS


function [CCall, LagsAll, PairedBirdID, PairedBranchID, binsize] = fn_calcxcov(windowmax, BirdID, BranchID, LocationAll, ...
    YallMat, XallMat, locationstoplot, suppressplots, doCConMeans)
if ~exist('suppressplots', 'var')
    suppressplots=0;
end

if ~exist('doCConMeans', 'var')
    doCConMeans=0;
end

%% ############################### CALCULATE CROSS CORRELATIONS
% for a given bird/branchID, all pairwise between neurons
% windowmax = 0.1; % in sec

numbirds = max(BirdID);
numbranches = max(BranchID);
numloc = max(LocationAll);

% --- for output
CCall = [];
LagsAll = [];
PairedBirdID = [];
PairedBranchID = [];


% --- RUN
for i=1:numbirds
    for bb = 1:numbranches
        
        inds = BirdID==i & BranchID==bb;
        if ~any(inds)
            continue
        end
        
        % ======================== EXTRACT DAT FOR THIS BRANCH
        loc = LocationAll(inds);
        ymat = YallMat(inds,:);
        xmat = XallMat(inds, :);
        binsize = xmat(1,2) - xmat(1,1);
        
        if all(isnan(loc))
            continue
        end
        
        % ======================= MAKE SURE DAT IS USABLE.
        if length(unique(loc))==1
            %  then this branch only has data from one brain region ...
            continue
        end
        
        if length(unique(loc)) > 2
            % -- then this branch has data from too manyr egions ...
            disp(unique(loc));
            keyboard
            disp('=== NOTE: skipped a branch because has too many brain regions (>2)');
            continue
        end
        
        % ====================== SEPARATE INTO 2 BRAIN REGIONS
        loclist = unique(loc)';
        
        % ----------- separation matrix for each brain region
        y1 = ymat(loc==loclist(1), :);
        y2 = ymat(loc==loclist(2), :);
        
        if doCConMeans==0
            for j = 1:size(y1,1)
                
                for jj=1:size(y2,1)
                    
                    [cc, lags] = xcov(y1(j,:), y2(jj,:), ceil(windowmax/binsize), 'coeff');
                    %                 [cc, lags] = xcorr(y1(j,:), y2(jj,:), ceil(windowmax/binsize), 'unbiased');
                    
                    
                    % ------------------------- OUTPUT
                    CCall = [CCall; cc];
                    LagsAll = [LagsAll; lags];
                    
                    PairedBirdID = [PairedBirdID; i];
                    PairedBranchID = [PairedBranchID; bb];
                    
                    
                end
                
            end
        else
            y1 = mean(y1, 1);
            y2 = mean(y2, 1);
            
            [cc, lags] = xcov(y1, y2, ceil(windowmax/binsize), 'coeff');
            %                 [cc, lags] = xcorr(y1, y2, ceil(windowmax/binsize), 'unbiased');
            
            % ------------------------- OUTPUT
            CCall = [CCall; cc];
            LagsAll = [LagsAll; lags];
            
            PairedBirdID = [PairedBirdID; i];
            PairedBranchID = [PairedBranchID; bb];
            
        end
    end
end

if suppressplots==0
    %% ======= PLOTS (BY BIRD)
    for i=1:numbirds
        figcount=1;
        subplotrows=3;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        
        
        % ========== all dat
        inds = PairedBirdID==i;
        if ~any(inds)
            continue
        end
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('all neuron pairs');
        plot(lags*binsize, CCall(inds,:)', 'Color', [0.7 0.7 0.7]);
        
        shadedErrorBar(lags*binsize, mean(CCall(inds,:), 1), lt_sem(CCall(inds,:)), {'Color', 'r'},1);
        xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
        lt_plot_zeroline_vert;
        
        title(['bird ' num2str(i)]);
        
        % ----------- find and plot peaks
        [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
        line([lags(loctmp)*binsize lags(loctmp)*binsize], ylim, 'Color', 'b');
        
        
        % ============== DO THE SAME FOR EACH BRANCH
        for j=1:numbranches
            
            
            inds = PairedBirdID==i & PairedBranchID==j;
            
            if all(inds==0)
                continue
            end
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['branch ' num2str(j)]);
            plot(lags*binsize, CCall(inds,:)', 'Color', [0.7 0.7 0.7]);
            
            try
                shadedErrorBar(lags*binsize, mean(CCall(inds,:), 1), lt_sem(CCall(inds,:)), {'Color', 'r'},1);
            catch err
            end
            xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
            lt_plot_zeroline_vert
            
            % ----------- find and plot peaks
            [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
            if ~isempty(loctmp)
                line([lags(loctmp)*binsize lags(loctmp)*binsize], ylim, 'Color', 'b');
            end
            
        end
        
        
        % ================== SAME, BUT JUST PLOT THE MEAN
        figcount=1;
        subplotrows=3;
        subplotcols=2;
        fignums_alreadyused=[];
        hfigs=[];
        
        inds = PairedBirdID==i;
        if ~any(inds)
            continue
        end
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title('all neuron pairs');
        %         plot(lags*binsize, CCall(inds,:)', 'Color', [0.7 0.7 0.7]);
        
        shadedErrorBar(lags*binsize, mean(CCall(inds,:), 1), lt_sem(CCall(inds,:)), {'Color', 'r'},1);
        xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
        lt_plot_zeroline_vert;
        
        title(['bird ' num2str(i)]);
        
        % ----------- find and plot peaks
        [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
        line([lags(loctmp)*binsize lags(loctmp)*binsize], ylim, 'Color', 'b');
        
        
        for j=1:numbranches
            
            inds = PairedBirdID==i & PairedBranchID==j;
            
            if all(inds==0)
                continue
            end
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['branch ' num2str(j)]);
            %             plot(lags*binsize, CCall(inds,:)', 'Color', [0.7 0.7 0.7]);
            
            try
                shadedErrorBar(lags*binsize, mean(CCall(inds,:), 1), lt_sem(CCall(inds,:)), {'Color', 'r'},1);
            catch err
            end
            xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
            lt_plot_zeroline_vert
            
            % ----------- find and plot peaks
            [~, loctmp] = findpeaks(mean(CCall(inds,:), 1), 'NPEAKS', 1, 'SORTSTR', 'descend');
            if ~isempty(loctmp)
                line([lags(loctmp)*binsize lags(loctmp)*binsize], ylim, 'Color', 'b');
            end
            
        end
        
        
    end
    
    
    %% ======= PLOTS (GRAND AVE)
    lt_figure; hold on;
    
    lt_subplot(2,2,1); hold on;
    title('all neuron pairs');
    plot(lags*binsize, CCall', 'b');
    xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
    
    lt_subplot(2,2,2); hold on;
    title('mean(sem)');
    if length(lt_sem(CCall))>1
        shadedErrorBar(lags*binsize, mean(CCall, 1), lt_sem(CCall), {'Color', 'k'},1);
        xlabel([locationstoplot{1} ' leads <---> ' locationstoplot{2} ' leads']);
    end
    
    % =========== GET MEAN, WEIGHTED BY SAMPEL SIZE
    if doCConMeans==0
        
        % === group: each level is single branch in single bird
        tmp = num2str([PairedBirdID PairedBranchID]);
        [GrpInds] = grp2idx(tmp);
        
        % --- count sample size for each group
        tmp = tabulate(GrpInds);
        weights = 1./tmp(:,2); % -- weighting factors are 1/N
        
        % -- normalize weights so that will add to one
        weights = weights./sum(weights(GrpInds));
        
        % -- take dot product of CC with weights
        weightvec = weights(GrpInds);
        
        lt_subplot(2,2,3); hold on;
        title('all branches equal weight');
        if (1)
            ccmean = CCall'*weightvec;
            %                 title('weighted mean (1/N)');
            plot(lags*binsize, ccmean);
            lt_plot_zeroline_vert;
        else
            for j=1:size(CCall,2)
                cc =  CCall(:, j)' * weightvec;
                cc = sum(CCall(:, j).*weightvec);
                plot(j, cc, '-ok');
            end
        end
    end
end


end

