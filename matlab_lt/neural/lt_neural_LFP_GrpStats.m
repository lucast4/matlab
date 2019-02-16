function [allbnum, allenum, allswnum, allDat, ...
    allswitch_bnum, allswitch_enum, allswitch_swnum, allswitch_dat, allNmotifs, allswitch_Nmotifs] = ...
    lt_neural_LFP_GrpStats(OUTSTRUCT, fieldtoget)
%% lt 2/7/19 - added output of sample size (i.e number of motifs)

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
if iscell(OUTSTRUCT.(fieldtoget))
    allDat = nan([size(OUTSTRUCT.(fieldtoget){1},1), size(OUTSTRUCT.(fieldtoget){1},2), 3, length(indsgrp_unique)]); % [t, ff, targ-same-diff, cases]
else
    allDat = nan([1, size(OUTSTRUCT.(fieldtoget),2), 3, length(indsgrp_unique)]); % [numrows(always 1), numcols, targ-same-diff, cases]
end
allNmotifs = nan([1, 3, length(indsgrp_unique)]); % [numrows(always 1), numcols, targ-same-diff, cases]

for i=1:length(indsgrp_unique)
    j=indsgrp_unique(i);
    
    % for this channel pair ...
    % ----- info...
    allbnum = [allbnum; unique(OUTSTRUCT.bnum(indsgrp==j))];
    allenum = [allenum; unique(OUTSTRUCT.enum(indsgrp==j))];
    allswnum = [allswnum; unique(OUTSTRUCT.switch(indsgrp==j))];
    
    % ---------- target
    indsthis = indsgrp==j & OUTSTRUCT.istarg==1;
    if iscell(OUTSTRUCT.(fieldtoget))
        datthis = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoget)(indsthis)), 3);
    else
        datthis= nanmean(OUTSTRUCT.(fieldtoget)(indsthis,:),1);
    end
    allDat(:,:,1,i) = datthis;
    allNmotifs(1,1,i) = sum(indsthis);
    
    % ---------- same-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==1;
    if any(indsthis)
        if iscell(OUTSTRUCT.(fieldtoget))
            datthis = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoget)(indsthis)), 3);
        else
            datthis= nanmean(OUTSTRUCT.(fieldtoget)(indsthis,:),1);
        end
        allDat(:,:,2,i) = datthis;
        allNmotifs(1,2,i) = sum(indsthis);
    else
        allNmotifs(1,2,i) = 0;
    end
    
    % ---------- diff-type
    indsthis = indsgrp==j & OUTSTRUCT.istarg==0 & OUTSTRUCT.issame==0;
    if any(indsthis)
        if iscell(OUTSTRUCT.(fieldtoget))
            datthis = nanmean(lt_neural_Coher_Cell2Mat(OUTSTRUCT.(fieldtoget)(indsthis)), 3);
        else
            datthis= nanmean(OUTSTRUCT.(fieldtoget)(indsthis,:),1);
        end
        allDat(:,:,3,i) = datthis;
        allNmotifs(1,3,i) = sum(indsthis);
    else
        allNmotifs(1,3,i) = 0;
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
allswitch_Nmotifs = nan(1, 3, length(indsgrp_switch_unique));

for i=1:length(indsgrp_switch_unique)
    %     j=indsgrp_switch_unique(i);
    
    indsthis = indsgrp_switch==indsgrp_switch_unique(i); % all channel pairs for this switch
    % ---------------------- INFO
    allswitch_bnum = [allswitch_bnum; unique(allbnum(indsthis))];
    allswitch_enum = [allswitch_enum; unique(allenum(indsthis))];
    allswitch_swnum =[allswitch_swnum; unique(allswnum(indsthis))];
    
    % ================= target
    sylind = 1;
    
    %     cohmat = nanmean(squeeze(allDat(:,:,sylind,indsthis)), 3);
    cohmat = squeeze(nanmean(allDat(:,:,sylind,indsthis), 4));
    allswitch_dat(:,:,sylind, i) = cohmat;
    allswitch_Nmotifs(1,sylind,i) = mean(allNmotifs(1,sylind,indsthis));
    assert(length(unique(squeeze(allNmotifs(1,sylind,indsthis))))==1);
    
    % ================= same
    sylind = 2;
    %     cohmat = nanmean(squeeze(allDat(:,:,sylind,indsthis)), 3);
    cohmat = squeeze(nanmean(allDat(:,:,sylind,indsthis), 4));
    allswitch_dat(:,:,sylind, i) = cohmat;
    allswitch_Nmotifs(1,sylind,i) = mean(allNmotifs(1,sylind,indsthis));
    if ~isnan(mean(allNmotifs(1,sylind,indsthis)))
        assert(length(unique(squeeze(allNmotifs(1,sylind,indsthis))))==1);
    end
    
    % ================= target
    sylind = 3;
    %     cohmat = nanmean(squeeze(allDat(:,:,sylind,indsthis)), 3);
    cohmat = squeeze(nanmean(allDat(:,:,sylind,indsthis), 4));
    allswitch_dat(:,:,sylind, i) = cohmat;
    allswitch_Nmotifs(1,sylind,i) = mean(allNmotifs(1,sylind,indsthis));
    if ~isnan(mean(allNmotifs(1,sylind,indsthis)))
        assert(length(unique(squeeze(allNmotifs(1,sylind,indsthis))))==1);
    end
end

