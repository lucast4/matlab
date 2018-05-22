function lt_neural_v2_CTXT_BRANCH_PlotOne(analyfname, birdtoget, branchtoget, ...
    neurtoget)

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
birdnum = find(strcmp({SummaryStruct.birds.birdname}, birdtoget));
numbranches = length(DATSTRUCT_BYBRANCH.bird(birdnum).branchID);

for j=1:numbranches

    if ~any(strcmp(DATSTRUCT_BYBRANCH.bird(birdnum).branchID(j).regexpstr{1}, branchtoget))
        continue
    end
    
    DAT = DATSTRUCT_BYBRANCH.bird(birdnum).branchID(j).DAT;
    nn = find(DAT.IndsOriginal_neuron == neurtoget);
    
    assert(length(nn)==1, 'cannot find this combo of bird, branch, neur...');
    
sts = lt_neural_ConfMatStats(DAT.PREMOTORDECODE_struct(nn).ConfMatAll_DAT{1});
y = sts.(decodestat);

sts = lt_neural_ConfMatStats(DAT.PREMOTORDECODE_struct(nn).ConfMatAll_NEG);
yneg = [sts.(decodestat)];


% ==================== PLOT
lt_figure; hold on;
lt_plot_histogram(yneg);
line([y y], ylim);
    

% ================== PLOT SHUFFLE DISTRIBUTION WITH ACTUAL DAT
yCI = prctile(yneg, [2.5 97.5]);
ymed = median(yneg);
lt_figure; hold on;
x = 1;
patch([x-0.5 x+0.5 x+0.5 x-0.5], [yCI(1) yCI(1) yCI(2) yCI(2)], [0.8 0.8 0.8]);
line([x-0.5 x+0.5], [ymed ymed]);

lt_plot(x, y);
xlim([0 2]);    
    % =============== plot decode timecourse
    if (1)
        lt_figure; hold on;
        y = DAT.ydecode(nn,:);
        x = DAT.xdecode;
        yneg = DAT.ydecode_neg(nn,:);
        plot(x, y, '-k');
        plot(x, yneg, '-r');
    end
end