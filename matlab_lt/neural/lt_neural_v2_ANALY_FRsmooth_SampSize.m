function lt_neural_v2_ANALY_FRsmooth_BaseMinu(OUTDAT, SwitchStruct, ...
    epochtoplot, analytype, doShuff, syltypesneeded, premotorwind, nshuffs, ...
    minmotifs)
%%

dattype = 'neuron'; % for breaking down shuffle into lower level data.
% dattype = , bird, 'switch', 'expt', 'neuron'
%     nshuffs = 500;

%% ####################### GETTING DEVIATION FROM BASELINE SMOOTHED FR

% prctile_divs = [33 66 100]; % percentiles to divide up data by
% epochtoplot = 3; % i.e. out of the epochs decided by prctile_divs
% prctile_divs = [50 100]; % percentiles to divide up data by
% epochtoplot = 2; % i.e. out of the epochs decided by prctile_divs

% usepercent = 0; % if 1, then gets fr percent deviation from baseline. otherwise uses hz diff
% nbasetime = 60; % 60 minutes bnefore first WN trial is min time for baseline
% nbasetime = []; % 60 minutes bnefore first WN trial is min time for baseline

% analytype = 'AllMinusAll_FRsmooth';
% % AllMinusAll_FRsmooth
% % AllOnlyMinusDiff_FRsmooth
% % AllMinusAllMinusDiff_FRsmooth
% % AllOnlyMinusBase_FRsmooth
%
%
% syltypesneeded = [1 1 1] means needs minimum 1 targ, 1 same, 1 diff.
% [1 0 1] means doesnt care if has same type

%% ###################### LIMIT TO NEURONS THAT CONTAIN ALL SYL TYPES?

% =-==== go thru all switches. if bad then throw ou
maxneur = max(OUTDAT.All_neurnum);
numbirds = length(SwitchStruct.bird);

% ======================== SHUFFLE SYL TYPE? % within each neuron
indstokeep = [];
for i=1:numbirds
    numexpts = length(SwitchStruct.bird(i).exptnum);
    for ii=1:numexpts
        numswitch = length(SwitchStruct.bird(i).exptnum(ii).switchlist);
        for ss = 1:numswitch
            % ################### SECOND, PLOT SMOOTH FR FOR ALL SYLS, SUBTRACT BASE
            for nn=1:maxneur
                
                indsthis = find(OUTDAT.All_birdnum==i & OUTDAT.All_exptnum==ii & OUTDAT.All_swnum==ss ...
                    & OUTDAT.All_neurnum==nn);
                
                if isempty(indsthis)
                    continue
                end
                
                % ==== count how many of each syl type there exists
                numtarg = sum(OUTDAT.All_istarg(indsthis)==1);
                numsame = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==1);
                numdiff = sum(OUTDAT.All_istarg(indsthis)==0 & OUTDAT.All_issame(indsthis)==0);
                
                if all([numtarg numsame numdiff] >= syltypesneeded)
                    % then keep
                    disp([numtarg numsame numdiff]);
                    disp(indsthis)
                    indstokeep = [indstokeep; indsthis];
                end
                
            end
        end
    end
end

disp(['Keeping ' num2str(length(indstokeep)) '/' num2str(length(OUTDAT.All_birdnum)) ' datapoints, passes syltypes required criterion']);

OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);


%% =================== ONLY KEEP SWITCHES THAT HAVE A MINIMUM NUMBER OF MOTIFS

[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTDAT.All_birdnum, OUTDAT.All_exptnum, OUTDAT.All_swnum, OUTDAT.All_neurnum});
indstoremove = [];
for i=indsgrpU'
    
    indsthis = indsgrp==i;
    
    nmot = length(unique(OUTDAT.All_motifnum(indsthis)));
    if nmot < minmotifs
        indstoremove = [indstoremove; find(indsthis)];
    end
    
end

indstokeep = ~ismember(1:length(OUTDAT.All_birdnum), indstoremove');
OUTDAT = lt_structure_subsample_all_fields(OUTDAT, indstokeep, 1);

%% ####################### SUMMARY - collect, for each subtract global mean

% ======= ACTYAL DAT (NOT SHUFF)
shuffSylType =0;
plotOn = 1;
[OUTDAT_tmp] = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, shuffSylType, ...
    epochtoplot, plotOn);



