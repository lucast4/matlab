function lt_neural_v2_CTXT_PlotAllRaw(ALLBRANCH, BirdToPlot, ExptToPlot, BrainRegions, ...
    plot_negshuff_z)
useDprime = 0;
analyfname = '';

%% lt 4/16/18 - plots each branch and neuron, decode accuracy, FR, pos and neg controls

apos =1;

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


%% ============ PLOT
if plot_negshuff_z ==1
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;

motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=4;
    subplotcols=3;
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
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_z = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.y_decode_z(inds,:);
            
            ydecode_neg_z_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.y_decode_negshuff_mean(inds,:);
            ydecode_neg_z_std = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.y_decode_negshuff_std(inds,:);
            ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
            
            sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(inds,:);
            sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(inds,:);
            sylcontours_mean = mean(sylcontours_mean,1);
            sylcontours_x = mean(sylcontours_x,1);
            
            % ############### 1) PLOT DECODE
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
            ylabel('decode(r), negmean(k), pos(g), 10*negstd(b)');
            
            if isempty(ydecode)
                continue
            end
            
            % -- lines
            line([0 0], ylim);
            
            % -- neg control
            plot(xdecode, ydecode_neg_z_mean, '-', 'Color', [0.5 0.5 0.5]);
            plot(xdecode, 10*ydecode_neg_z_std, '-', 'Color', 'b');
            
            % -- pos control
            plot(xdecode, ydecode_pos, '-', 'Color', [0.3 0.8 0.2]);
            
            % -- data
            plot(xdecode, ydecode, '-', 'Color', 'r');
            plot(xdecode, mean(ydecode,1), '-r', 'LineWidth', 2)
            
            plot(sylcontours_x, 0.25*sylcontours_mean+1, '-m');
            
            
            
            % ###################### 2) PLOT DECODE USING ZSCORE
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            ylabel('z-scored decode (r) [medinan]');
            
            line([0 0], ylim);
            
            plot(xdecode, ydecode_z, '-', 'Color', 'r');
            plot(xdecode, median(ydecode_z, 1), '-r', 'LineWidth', 2)
            
            
            % ##################### 3) PLOT FIRING RATES
            neurIDs = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(inds);
            branchall = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_branch(inds);
            
            FRmean_allall = [];
            FRstd_allall = [];
            for n = 1:length(neurIDs)
                nn = neurIDs(n);
                branchthis = branchall(n);
                
                dattmp = ALLBRANCH.alignpos(apos).bird(i).branch(branchthis).neuron(nn);
                
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
                frmeanTMP = mean(FRmean_all,1);
                frstdTMP = std(FRmean_all, 0, 1);
                
                FRmean_allall = [FRmean_allall; frmeanTMP];
                FRstd_allall = [FRstd_allall; frstdTMP];
            end
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            
            ylabel('mean/STD FR (across branches and neurons for each timepoint');
            y = mean(FRmean_allall,1);
            if size(FRstd_allall,1)>1
                ystd = mean(FRstd_allall, 1);
            end
            x = frx(1:length(mean(FRmean_allall,1)));
            shadedErrorBar(x-motifpredur, y, ystd, {'Color', 'k'},1)
            
            
        end
        linkaxes(hsplots, 'x');
    end
    axis tight
end

else
%% ============ PLOT
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;
apos = 1;

motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;

for i=1:numbirds
    birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
    figcount=1;
    subplotrows=5;
    subplotcols=4;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots =[];
    
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        for loc = locationstoplot
            
            inds = find(strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc));
            if ~any(inds)
                continue
            end
            
            for nn=inds'
                
                % -------------- PLOT
                xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
                ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(nn,:);
                ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(nn,:);
                ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(nn,:);
                neurind = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron(nn);
                branchid = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_branch(nn);
                
                sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean(nn,:);
                sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x(nn,:);
                
                
                % --------- count number of branches here
                Nclass = size(ALLBRANCH.alignpos(apos).bird(i).branch(branchid).neuron(neurind).ConfMatAll,1);
                
                % ############### 1) PLOT DECODE
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                title([birdname '-' thisbranch{1} '-n' num2str(neurind) '[' loc{1} ']']);
                ylabel('decode(r), neg(k), pos(g)');
                
                if isempty(ydecode)
                    continue
                end
                
                % -- lines
                line([0 0], ylim);
                
                % -- neg control
                plot(xdecode, ydecode_neg, '-', 'Color', [0.5 0.5 0.5]);
                
                % -- pos control
                plot(xdecode, ydecode_pos, '-', 'Color', [0.3 0.8 0.2]);
                
                % -- data
                plot(xdecode, ydecode, '-', 'Color', 'r');
                
                plot(sylcontours_x, 0.25*sylcontours_mean+1, '-m');
                
                lt_plot_text(0, 0.8, ['N=' num2str(Nclass)], 'm');
                
                % ##################### 3) PLOT FIRING RATES
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                hsplots = [hsplots hsplot];
                
                dattmp = ALLBRANCH.alignpos(apos).bird(i).branch(branchid).neuron(neurind);
                
                % ==== collect original FR
                numclasses = length(dattmp.FR.classnum);
                plotcols = lt_make_plot_colors(numclasses,0,0);
                for cc = 1:numclasses
                    
                    frmat = dattmp.FR.classnum(cc).FRsmooth_rate_CommonTrialDur;
                    frx = dattmp.FR.classnum(cc).FRsmooth_xbin_CommonTrialDur;
                    
                    frmean = mean(frmat,2);
                    frstd = std(frmat,0, 2);
                    
                    shadedErrorBar(frx-motifpredur, frmean, frstd, {'Color', plotcols{cc}},1)
                    
                end
                
            end
        end
        linkaxes(hsplots, 'x');
    end
    axis tight
end
end