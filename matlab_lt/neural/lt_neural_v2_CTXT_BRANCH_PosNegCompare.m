function lt_neural_v2_CTXT_BRANCH_PosNegCompare(ALLBRANCH, BrainRegions, ...
    BirdToPlot, useDprime, ExptToPlot, analyfname, sortbypremotorpval, ...
    plotAllZscoreFR, averageWithinNeuron, minbranches, useCorrMetric)
apos =1;

timewindow = [-0.03 0.03]; % relative to syllable onset


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
        useDprime, timewindow);
    
else
    DATSTRUCT_BYBRANCH = lt_neural_v2_CTXT_BRANCH_OrgByBranchID(ALLBRANCH, '', '', ...
        useDprime, timewindow);
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

inds = All_birdID==1;

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

All_frcorr = [];
All_frcorr_neg = [];
All_frcorr_pos = [];

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
        
        frcorr_dat = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.FRcorr_dat;
        frcorr_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.FRcorr_neg;
        frcorr_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.FRcorr_pos;
        
        % ============= COLLECT
        All_Y = [All_Y; ydecode];
        All_Yneg = [All_Yneg; ydecode_neg];
        All_Ypos = [All_Ypos; ydecode_pos];
        All_Bregion = [All_Bregion; bregions];
        All_branchID = [All_branchID; bb*ones(size(bregions))];
        All_birdID = [All_birdID; i*ones(size(bregions))];
        All_neurID = [All_neurID; neurID];
        All_X = [All_X; xdecode];
       
        All_frcorr = [All_frcorr; frcorr_dat];
        All_frcorr_neg = [All_frcorr_neg; frcorr_neg];
        All_frcorr_pos = [All_frcorr_pos; frcorr_pos];

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
All_frcorr_bn = [];
All_frcorr_neg_bn = [];
All_frcorr_pos_bn = [];

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
        
        frcorr = All_frcorr(indstmp);
        frcorrNEG = All_frcorr_neg(indstmp);
        frcorrPOS = All_frcorr_pos(indstmp);
        
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
        
        All_frcorr_bn = [All_frcorr_bn; mean(frcorr)];
All_frcorr_neg_bn = [All_frcorr_neg_bn; mean(frcorrNEG)];
All_frcorr_pos_bn = [All_frcorr_pos_bn; mean(frcorrPOS)];

    end
    
    
    % =================== OVERWRITE DATA
    All_Y = All_Y_byneuron;
    All_Yneg = All_Yneg_byneuron;
    All_Ypos = All_Ypos_byneuron;
    All_Bregion = All_Bregion_byneuron;
    clear All_branchID; % since have averaged across branches
    All_birdID = All_birdID_byneuron;
    All_neurID = All_neurID_byneuron;
    
            All_frcorr = All_frcorr_bn;
        All_frcorr_neg = All_frcorr_neg_bn;
        All_frcorr_pos = All_frcorr_pos_bn;


end


xdecode = All_X;
%% ===================== PLOTS

for i=1:numbirds
    lt_figure; hold on;
    
    % ================= all data, LMAN vs. RA
    lt_subplot(2,2,1); hold on;
    title('LMAN');
    ylabel('vs. neg');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'LMAN');
    
    ydat = All_Y(inds,:);
    yneg = All_Yneg(inds,:);
    ypos = All_Ypos(inds,:);
    
    % -- subtract neg controls
    ydat = ydat-yneg;
    ypos = ypos-yneg;
    shadedErrorBar(xdecode, mean(ydat,1), std(ydat), {'Color', 'k'},1);
    shadedErrorBar(xdecode, mean(ypos,1), std(ypos), {'Color', 'b'},1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    % ================= all data, LMAN vs. RA
    lt_subplot(2,2,2); hold on;
    title('RA');
    ylabel('vs. neg');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'RA');
    
    ydat = All_Y(inds,:);
    yneg = All_Yneg(inds,:);
    ypos = All_Ypos(inds,:);
    
    % -- subtract neg controls
    ydat = ydat-yneg;
    ypos = ypos-yneg;
    shadedErrorBar(xdecode, mean(ydat,1), std(ydat), {'Color', 'k'},1);
    shadedErrorBar(xdecode, mean(ypos,1), std(ypos), {'Color', 'b'},1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5 SAME, WITHOUT
    % SUBTRACTING NEG CONTROLS
    % ================= all data, LMAN vs. RA
    lt_subplot(2,2,3); hold on;
    title('LMAN');
    ylabel('decode');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'LMAN');
    
    ydat = All_Y(inds,:);
    yneg = All_Yneg(inds,:);
    ypos = All_Ypos(inds,:);
    
    % -- subtract neg controls
    shadedErrorBar(xdecode, mean(ydat,1), std(ydat), {'Color', 'k'},1);
    shadedErrorBar(xdecode, mean(ypos,1), std(ypos), {'Color', 'b'},1);
    shadedErrorBar(xdecode, mean(yneg,1), std(yneg), {'Color', 'r'},1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    ylim([0 1]);
    
    % ================= all data, LMAN vs. RA
    lt_subplot(2,2,4); hold on;
    title('RA');
    ylabel('decode');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'RA');
    
    ydat = All_Y(inds,:);
    yneg = All_Yneg(inds,:);
    ypos = All_Ypos(inds,:);
    
    % -- subtract neg controls
    shadedErrorBar(xdecode, mean(ydat,1), std(ydat), {'Color', 'k'},1);
    shadedErrorBar(xdecode, mean(ypos,1), std(ypos), {'Color', 'b'},1);
    shadedErrorBar(xdecode, mean(yneg,1), std(yneg), {'Color', 'r'},1);
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    ylim([0 1]);
    
    
end

%% ========================= ONLY PREMOTOR WINDOW, MODEL TESTING

if useCorrMetric==1
% then use 1 minus corr of premotor FR

All_Y_point = 1-All_frcorr;
All_Yneg_point = 1-All_frcorr_neg;
All_Ypos_point = 1-All_frcorr_pos;

elseif useCorrMetric==0
    % then use mean decode in premotor window
assert(length(xdecode) == size(All_Y,2), 'asdfasd');
timebinstake = xdecode>=timewindow(1) & xdecode<timewindow(2);

All_Y_point = mean(All_Y(:, timebinstake), 2);
All_Yneg_point = mean(All_Yneg(:, timebinstake), 2);
All_Ypos_point = mean(All_Ypos(:, timebinstake), 2);
end

%% ===================== distribution of points

for i=1:numbirds
    lt_figure; hold on;
    
    % ================= all data, LMAN vs. RA
    hsplot = lt_subplot(2,2,1); hold on;
    hsplots = [hsplots hsplot];
    title('LMAN');
    xlabel('neg - dat - pos');
    ylabel('1 minus corr (premotor FR)');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'LMAN');
    
    ydat = All_Y_point(inds,:);
    yneg = All_Yneg_point(inds,:);
    ypos = All_Ypos_point(inds,:);
    x = [1 2 3];
    plot(x, [yneg ydat ypos], '-ok');
    
    
    
    % ================= all data, LMAN vs. RA
    hsplot = lt_subplot(2,2,2); hold on;
    hsplots = [hsplots hsplot];
    title('RA');
    xlabel('neg - dat - pos');
    ylabel('1 minus corr (premotor FR)');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'RA');
    
    ydat = All_Y_point(inds,:);
    yneg = All_Yneg_point(inds,:);
    ypos = All_Ypos_point(inds,:);
    x = [1 2 3];
    plot(x, [yneg ydat ypos], '-ok');
end


%% ===================== PLOTS
hsplots = [];

for i=1:numbirds
    lt_figure; hold on;
    
    % ================= all data, LMAN vs. RA
    hsplot = lt_subplot(2,2,1); hold on;
    hsplots = [hsplots hsplot];
    title('LMAN');
    xlabel('pos minus neg');
    ylabel('dat minus neg');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'LMAN');
    
    ydat = All_Y_point(inds,:);
    yneg = All_Yneg_point(inds,:);
    ypos = All_Ypos_point(inds,:);
    
    % -- subtract neg controls
    ydat = ydat-yneg;
    ypos = ypos-yneg;
    plot(ypos, ydat, 'ok');
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-0.3 1], [-0.3 1]);
    
    % ================= all data, LMAN vs. RA
    hsplot = lt_subplot(2,2,2); hold on;
    hsplots = [hsplots hsplot];
    title('RA');
    xlabel('pos minus neg');
    ylabel('dat minus neg');
    
    inds = All_birdID==i & strcmp(All_Bregion, 'RA');
    
    ydat = All_Y_point(inds,:);
    yneg = All_Yneg_point(inds,:);
    ypos = All_Ypos_point(inds,:);
    
    % -- subtract neg controls
    ydat = ydat-yneg;
    ypos = ypos-yneg;
    plot(ypos, ydat, 'ok');
    lt_plot_zeroline;
    lt_plot_zeroline_vert;
    line([-0.3 1], [-0.3 1]);
end

linkaxes(hsplots, 'xy');


%% ============= linear model

% dat(minus neg) ~ dat(pos) + [random effect of bird and branchID on slope]

ydat = double(All_Y_point - All_Yneg_point);
ypos = double(All_Ypos_point - All_Yneg_point);

if averageWithinNeuron==0
branchID = num2str([All_birdID All_branchID]);
tbl = table(ydat, ypos, branchID, All_Bregion, All_birdID, 'VariableNames', {'ydat', 'ypos', 'branchbird','bregion', 'birdID'});

% mdl = 'ydat ~ ypos + bregion:ypos + (-1 + ypos|branchbird) + (-1+bregion:ypos|branchbird)';  % full model
mdl = 'ydat ~ ypos + bregion:ypos + (-1 + ypos|branchbird)';  % reduced
% mdl = 'ydat ~ ypos + bregion:ypos + (ypos|birdID) + (-1+bregion:ypos|birdID)';

elseif averageWithinNeuron==1
tbl = table(ydat, ypos, All_Bregion, All_birdID, 'VariableNames', {'ydat', 'ypos', 'bregion', 'birdID'});

mdl = 'ydat ~ ypos + bregion:ypos + (-1 + ypos|birdID) + (-1+bregion:ypos|birdID)';  % full model
% mdl = 'ydat ~ ypos + bregion:ypos + (-1+bregion:ypos|birdID)';  % reduced
    
end
lme = fitlme(tbl, mdl)


% ------------------- plot distribution of estimates for a given input
birdtoget = 2;


%% ############################################################
%% ###########################################################
%% ################# WHITTLE DOWN DATA (EACH NEURON ONE DATAPOINT)


%% =============== PLOT MEAN BY BIRD AND BRAIN REGION


%%




%%
%% ===== COLLECT DATA INTO VECTOR
% numbirds = length(DATSTRUCT_BYBRANCH.bird);
% locationstoplot = BrainRegions;
%
%
% for i=1:numbirds
%     birdname = ALLBRANCH.SummaryStruct.birds(i).birdname;
%
%     numbranches = length(DATSTRUCT_BYBRANCH.bird(i).branchID);
%
%     AllXdecode = [];
%     AllYdecode = [];
%     AllYdecode_neg = [];
%     AllYdecode_pos = [];
%     AllNeurInds = [];
%     AllBregions = {};
%
%     for bb=1:numbranches
%
%         thisbranch = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).regexpstr;
%
%         % ========== COLLECT
%         xdecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.xdecode(1,:);
%         ydecode = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode;
%         ydecode_neg = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_neg;
%         ydecode_pos = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.ydecode_pos;
%         neurInds = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.IndsOriginal_neuron;
%         bregions = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.brainregion;
%
%         sylcontours_mean = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_mean;
%         sylcontours_x = DATSTRUCT_BYBRANCH.bird(i).branchID(bb).DAT.sylcontour_x;
%         sylcontours_mean = mean(sylcontours_mean,1);
%         sylcontours_x = mean(sylcontours_x,1);
%
%         % ===========
%         AllXdecode = [AllXdecode; xdecode];
%         AllYdecode = [AllYdecode; ydecode];
%         AllYdecode_neg = [AllYdecode_neg; ydecode_neg];
%         AllYdecode_pos = [AllYdecode_pos; ydecode_pos];
%         AllNeurInds = [AllNeurInds; neurInds];
%         AllBregions = [AllBregions; bregions];
%     end
%
%     Xdecode = AllXdecode(1,:);
%
%     % ================= GO THRU EACH NEURON AND PLOT
%     for loc = locationstoplot
%
%         neuronstoplot = AllNeurInds(strcmp(AllBregions, loc));
%         neuronstoplot = unique(neuronstoplot); % list of neurons that have data for this bregion
%
%         % ===========
%         for n=1:length(neuronstoplot)
%             nn = neuronstoplot(n);
%             [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
%             hsplots = [hsplots hsplot];
%             title([birdname '-n' num2str(nn) ' [' loc{1} ']']);
%
%
%             indstmp = AllNeurInds == nn;
%
%             x = Xdecode;
%             ymat = AllYdecode(indstmp, :);
%             yneg = AllYdecode_neg(indstmp, :);
%             ypos = AllYdecode_pos(indstmp, :);
%
%             plot(x, ymat, '-k');
%             plot(x, mean(ymat,1), '-k', 'LineWidth', 2);
%
%             plot(x, yneg, '-r');
%             plot(x, mean(yneg,1), '-r', 'LineWidth', 2);
%
%             plot(x, ypos, '-b');
%             plot(x, mean(ypos,1), '-b', 'LineWidth', 2);
%
%             line([0 0], ylim, 'LineStyle', '--', 'Color', 'y');
%         end
%     end
% end
