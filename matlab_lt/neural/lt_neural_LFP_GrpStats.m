function [allbnum, allenum, allswnum, allDat, ...
    allswitch_bnum, allswitch_enum, allswitch_swnum, allswitch_dat] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget)
%% lt 11/5/18 - makes groups and takes average (splitting into target, same, diff)

% fieldtoget = 'Spec1Mean_Base';

%% EXTRACT (DATAPOINT = CHANNEL PAIR) - IE.. AVERAGES OVER SIMILAR TYPE SYLLABLES
% ======= for each channel pair, get mean for targ, nontarg
% -- gets motifs for each channel pair
indsgrp = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch, OUTSTRUCT.chanpair});
indsgrp_unique = unique(indsgrp);

allbnum = [];
allenum = [];
allswnum = [];
allDat = nan([size(OUTSTRUCT.(fieldtoget){1},1), size(OUTSTRUCT.(fieldtoget){1},2), 3, length(indsgrp_unique)]); % [t, ff, targ-same-diff, cases]

for i=1:length(indsgrp_unique)
    j=indsgrp_unique(i);
    
    % for this channel pair ...
    % ----- info...
    allbnum = [allbnum; unique(OUTSTRUCT.bnum(indsgrp==j))];
    allenum = [allenum; unique(OUTSTRUCT.enum(indsgrp==j))];
    allswnum = [allswnum; unique(OUTSTRUCT.switch(indsgrp==j))];
    
    % ---------- target
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    datthis = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoget)(indsthis)), 3);
    allDat(:,:,1,i) = datthis;
    
    
    % ---------- same-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    if any(indsthis)
        datthis = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoget)(indsthis)), 3);
        allDat(:,:,2,i) = datthis;
    end
    
    % ---------- diff-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    if any(indsthis)
        datthis = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoget)(indsthis)), 3);
        allDat(:,:,3,i) = datthis;
    end
    
end

if nanmean(allDat(:))>1
    keyboard
end

%% ============ further group by switch? (i.e. average over all chanel pairs)
% ==================== PLOT, ONE FOR EACH BIRD (overlay all switches)0
indsgrp_switch = lt_tools_grp2idx({allbnum, allenum, allswnum});
indsgrp_switch_unique = unique(indsgrp_switch);

allswitch_bnum = [];
allswitch_enum = [];
allswitch_swnum =[];
%     allswitch_dat = nan(length(tbins), length(ffbins), 3, length(indsgrp_switch_unique)); % [t, ff, targ-same-diff, cases]
allswitch_dat = nan([size(allDat,1), size(allDat,2), 3, length(indsgrp_switch_unique)]);

for i=1:length(indsgrp_switch_unique)
    j=indsgrp_switch_unique(i);
    
    indsthis = indsgrp_switch==j; % all channel pairs for this switch
    % ---------------------- INFO
    allswitch_bnum = [allswitch_bnum; unique(allbnum(indsthis))];
    allswitch_enum = [allswitch_enum; unique(allenum(indsthis))];
    allswitch_swnum =[allswitch_swnum; unique(allswnum(indsthis))];
    
    % ================= target
    sylind = 1;
    cohmat = nanmean(squeeze(allDat(:,:,sylind,indsthis)), 3);
    allswitch_dat(:,:,sylind, i) = cohmat;
    
    % ================= same
    sylind = 2;
    cohmat = nanmean(squeeze(allDat(:,:,sylind,indsthis)), 3);
    allswitch_dat(:,:,sylind, i) = cohmat;
    
    % ================= target
    sylind = 3;
    cohmat = nanmean(squeeze(allDat(:,:,sylind,indsthis)), 3);
    allswitch_dat(:,:,sylind, i) = cohmat;
    
end

