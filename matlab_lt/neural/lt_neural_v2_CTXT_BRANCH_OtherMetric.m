function lt_neural_v2_CTXT_BRANCH_PosNegCompare(ALLBRANCH, BrainRegions, ...
    BirdToPlot, useDprime, ExptToPlot, analyfname, sortbypremotorpval, ...
    plotAllZscoreFR, averageWithinNeuron, minbranches)
%% lt 4/23/18 - instead of using decode, use other metric, calculated 
% using premotor window

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

%% ===== COLLECT DATA INTO VECTOR
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;

timetoplot = -0; % relative to syl onset

All_Y = [];
All_Bregion = {};
All_branchID = [];
All_birdID = [];


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
        Yall = [];
        Bregionall = {};
        
        for loc = locationstoplot
            
            inds = strcmp(DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion, loc);
            if ~any(inds)
                continue
            end
            
            bregions = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion(inds);
            xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
            ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode(inds,:);
            ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg(inds,:);
            ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos(inds,:);
            % --- extract value at timepoint
            [~, indtmp] = min(abs(xdecode-timetoplot));
            
            ydecode = ydecode(:,indtmp);
            ydecode_neg = ydecode_neg(:, indtmp);
            ydecode_pos = ydecode_pos(:, indtmp);
            
            % -------------- PLOT
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            hsplots = [hsplots hsplot];
            title([birdname '-' thisbranch{1} ' [' loc{1} ']']);
            x = [1 2 3];
            y = [ydecode_neg ydecode ydecode_pos];
            plot(x, y, '-ok');
            xlim([0 4]);
            lt_plot_zeroline;
            ylim([-0.05 1]);
            
            % ================== OUTPUT
            % ------ LOCAL
            Yall = [Yall; y];
            Bregionall = [Bregionall; bregions];
            
            % ----- GLOBAL
            All_Y = [All_Y; y];
            All_Bregion = [All_Bregion; bregions];
            All_branchID = [All_branchID; bb*ones(size(bregions))];
            All_birdID = [All_birdID; i*ones(size(bregions))];
            
            
        end
        
        % ================== each branch
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title([birdname '-' thisbranch{1}]);
        % --- scatterplot
        title('RA(rd) and LMAN(k)');
        xlabel('pos minus neg control');
        ylabel('dat minus neg control');
        Ytmp = Yall - Yall(:,1);
        
        % --- lman
        indtmp = strcmp(Bregionall, 'LMAN');
        plot(Ytmp(indtmp, 3), Ytmp(indtmp,2), 'ok');
        % --- ra
        indtmp = strcmp(Bregionall, 'RA');
        plot(Ytmp(indtmp, 3), Ytmp(indtmp,2), 'or');
        
        xlim([-1 1]);
        ylim([-1 1]);
        line([0 1], [0 1]);
        lt_plot_zeroline;
        lt_plot_zeroline_vert;
        
    end
    
    % ================ PLOT FOR THIS BIRD, SUMMARY
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title([birdname '- ALL BRANCH']);
    title('RA(rd) and LMAN(k)');
    xlabel('pos minus neg control');
    ylabel('dat minus neg control');
    
    Ytmp = All_Y - All_Y(:,1);
    
    % --- lman
    indtmp = strcmp(All_Bregion, 'LMAN') & All_birdID==i;
    plot(Ytmp(indtmp, 3), Ytmp(indtmp,2), 'xk');
    % --- ra
    indtmp = strcmp(All_Bregion, 'RA') & All_birdID==i;
    plot(Ytmp(indtmp, 3), Ytmp(indtmp,2), 'xr');
    
    
    xlim([-1 1]);
    ylim([-1 1]);
    line([0 1], [0 1]);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
end

% ============== PLOT GLOBALLY
lt_figure; hold on;
% ----------------- 1) ALL DATA
lt_subplot(2,2,1); hold on
title('[all data] RA(rd) and LMAN(k)');
xlabel('pos minus neg control');
ylabel('dat minus neg control');
Ytmp = All_Y - All_Y(:,1);

% --- lman
indtmp = strcmp(All_Bregion, 'LMAN');
plot(Ytmp(indtmp, 3), Ytmp(indtmp,2), 'xk');
% --- ra
indtmp = strcmp(All_Bregion, 'RA');
plot(Ytmp(indtmp, 3), Ytmp(indtmp,2), 'xr');


xlim([-1 1]);
ylim([-1 1]);
line([0 1], [0 1]);
lt_plot_zeroline;
lt_plot_zeroline_vert;

% ---------------- 2)


%% quickly, analysis of covariance, are slopes different?

inds = All_birdID==2;

x = All_Y(inds,3)-All_Y(inds,1); % pos minus neg
y = All_Y(inds,2)-All_Y(inds,1); % dat minus neg
grpid = All_Bregion(inds);

aoctool(x, y, grpid)


%% #######################################################
%% #########################################################

%% ===== COLLECT ALL DATA INTO VECTOR
numbirds = length(DATSTRUCT_BYBRANCH.bird);
locationstoplot = BrainRegions;

All_X = [];
All_Y = [];
All_Yneg = [];
All_Ypos = [];
All_Bregion = {};
All_branchID = [];
All_birdID = [];
All_neurID = [];



for i=1:numbirds
    numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
    
    for bb=1:numbranches
        
        thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
        
        % ========= plot for this branch
        bregions = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion;
        xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode;
        ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode;
        ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg;
        ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos;
        neurID = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron;
        
        % ============= COLLECT
        All_Y = [All_Y; ydecode];
        All_Yneg = [All_Yneg; ydecode_neg];
        All_Ypos = [All_Ypos; ydecode_pos];
        All_Bregion = [All_Bregion; bregions];
        All_branchID = [All_branchID; bb*ones(size(bregions))];
        All_birdID = [All_birdID; i*ones(size(bregions))];
        All_neurID = [All_neurID; neurID];
        All_X = [All_X; xdecode];
    end
end

All_X= All_X(1,:);

%% =================== COLLECT DATA
if averageWithinNeuron==1
    
    figcount=1;
    subplotrows=4;
    subplotcols=2;
    fignums_alreadyused=[];
    hfigs=[];
    hsplots = [];
    
    
    All_Y_byneuron = [];
    All_Yneg_byneuron = [];
    All_Ypos_byneuron = [];
    All_IndsOrig_byneuron = {};
    All_Bregion_byneuron = {};
    All_birdID_byneuron = [];
    All_neurID_byneuron = [];
    
    All_UniqueID = mat2cell(num2str([All_birdID All_neurID]), ones(size(All_birdID,1),1));
    idlist = unique(All_UniqueID); % each individual unit
    
    for j=1:length(idlist)
        
        idthis = idlist{j};
        indstmp = find(strcmp(All_UniqueID, idthis));
        
        if length(indstmp)<minbranches
            disp('SKIPPED NEURON - not enouigh branches');
            continue
        end
        
        y = All_Y(indstmp,:);
        yneg = All_Yneg(indstmp,:);
        ypos = All_Ypos(indstmp,:);
        x = All_X(1,:);
        
        % ========== plot, overlay all traces for this neuron (across
        % branches)
        
        birdnum = unique(All_birdID(indstmp));
        birdname = ALLBRANCH.SummaryStruct.birds(birdnum).birdname;
        neurnum = unique(All_neurID(indstmp));
        bregion = unique(All_Bregion(indstmp));
        bregion = bregion{1};
        
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        if strcmp(bregion, 'RA')
            title([birdname '-n' num2str(neurnum) '[' bregion ']'], 'Color', 'r');
        elseif strcmp(bregion, 'LMAN')
            title([birdname '-n' num2str(neurnum) '[' bregion ']'], 'Color', 'g');
        end
        
        % ----- data
        plot(x, y', '-k');
        plot(x, mean(y,1), '-k', 'LineWidth', 2);
        
        % ----- neg
        plot(x, yneg', '-r');
        plot(x, mean(yneg,1), '-r', 'LineWidth', 2);
        
        % ----- pos
        plot(x, ypos', '-b');
        plot(x, mean(ypos,1), '-b', 'LineWidth', 2);
        
        axis tight
        lt_plot_zeroline_vert;
        
        
        % ================== COLLECT DATA
        All_Y_byneuron = [All_Y_byneuron; mean(y,1)];
        All_Yneg_byneuron = [All_Yneg_byneuron; mean(yneg,1)];
        All_Ypos_byneuron = [All_Ypos_byneuron; mean(ypos,1)];
        All_IndsOrig_byneuron = [All_IndsOrig_byneuron; indstmp];
        
        All_Bregion_byneuron =[All_Bregion_byneuron; bregion];
        All_birdID_byneuron = [All_birdID_byneuron; birdnum];
        All_neurID_byneuron = [All_neurID_byneuron; neurnum];
        
    end
    
    
    % =================== OVERWRITE DATA
    All_Y = All_Y_byneuron;
    All_Yneg = All_Yneg_byneuron;
    All_Ypos = All_Ypos_byneuron;
    All_Bregion = All_Bregion_byneuron;
    clear All_branchID; % since have averaged across branches
    All_birdID = All_birdID_byneuron;
    All_neurID = All_neurID_byneuron;
end


xdecode = All_X;