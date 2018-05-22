function lt_neural_v2_BRANCH_DecodePlot(DecodeStruct)
anum = 1;
plottext =0;

%% ================== [PLOT] effect sizes (i.e. decode F1) for all

% ------------------
DatAll = DecodeStruct.analynum(anum).dat;
DATSTRUCT_BYBRANCH = DecodeStruct.analynum(anum).datbybranch;
SummaryStruct = DecodeStruct.analynum(anum).SummaryStruct;

Numbirds = length(DecodeStruct.analynum(anum).datbybranch.bird);
Numbranch = max(DatAll.AllBranchNum);

BrainRegions = unique(DecodeStruct.analynum(anum).dat.AllBrainRegion);
plotcols_brainreg = lt_make_plot_colors(length(BrainRegions), 0,0);


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
        IsSU = {};
        
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