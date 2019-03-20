function OUTDAT = lt_neural_v2_ANALY_FRsmooth_NormFR(OUTDAT, prewind, ignorecontext, ...
    epochtoplot, meanver, SwitchStruct, usepercent, nbasetime, ...
    nbasetime_ignoreswitch1, prctile_divs)

%% lt 3/16/19 - for a given neuron, scale all FRsmooth so that on aveage is similar pre and post learning
disp('NOTE: for a given neuron, will scale each timebin and sylalble by identical factor for all WN trials');
disp('(even though gain factor will be computed using only WN epoch trials)');


% ignorecontext = 1; % then will take mean over contexts, so that each unique syl is one datapoint
% meanver = 0; % 0: first take avearge, then take ratio; 1: first take ratios, then take average.

if ignorecontext==0
    disp('NOT YET CODED')
    pause;
end
%%
[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTDAT.All_birdnum, OUTDAT.All_exptnum, ...
    OUTDAT.All_swnum, OUTDAT.All_neurnum});

% -- convert to scalar (mean FR)
frt = OUTDAT.All_FRsmooth_t{1,1};
indst = frt>prewind(1) & frt<prewind(2);

% -- keep track of all inds that have been done
indsDone = [];
for i=1:length(indsgrpU)
    indsthis = indsgrp==indsgrpU(i);
    indsDone = [indsDone; find(indsthis)];
    
    % == base
    indsbase = OUTDAT.AllBase_indsepoch(indsthis);
    frmat = OUTDAT.All_FRsmooth(indsthis,1);
    frscal_base = [];
    for j=1:length(frmat)
        tmp = mean(mean(frmat{j}(indst, indsbase{j}),1),2);
        frscal_base = [frscal_base; tmp];
    end
    
    % == wn
    indsWN = OUTDAT.AllWN_indsepoch(indsthis);
    frmat = OUTDAT.All_FRsmooth(indsthis,2);
    frscal_wn = [];
    for j=1:length(frmat)
        indstmp = indsWN{j}(epochtoplot):indsWN{j}(epochtoplot+1)-1;
        tmp = mean(mean(frmat{j}(indst,indstmp),1),2);
        frscal_wn = [frscal_wn; tmp];
    end
    
    if meanver==0
        % -- first mean, then average
        gainfact = mean(frscal_base)/mean(frscal_wn);
    elseif meanver==1
        gainfact = mean(frscal_base./frscal_wn);
    end
    
    % ============== RESACLE ALL WN FR TO MATCH BASELINE
    frmat = OUTDAT.All_FRsmooth(indsthis,2);
    frmat = cellfun(@(x)x.*gainfact, frmat, 'UniformOutput', 0);
    OUTDAT.All_FRsmooth(indsthis,2) = frmat;
end

assert(length(unique(indsDone))==length(OUTDAT.All_birdnum));

%% ============= RECOMPUTE THINGS FOR WN RELATIVE TO BASELINE

OUTDAT = lt_neural_v2_ANALY_FRsmooth_MinusBase(OUTDAT, SwitchStruct,...
    usepercent, nbasetime, nbasetime_ignoreswitch1, prctile_divs);
