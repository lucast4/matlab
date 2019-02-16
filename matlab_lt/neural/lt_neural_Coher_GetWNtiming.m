function [OUTSTRUCT, allwntimes] = lt_neural_Coher_GetWNtiming(OUTSTRUCT, onlyusegoodtargsyls, ...
    useallwn, prctiletouse, SwitchStruct, SwitchCohStruct, PARAMS)
%% lt 2/2019 - extract WN time (using some percentile of onsets)

%%
[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});

allwntimes = []; % 2.5 PERCENTILEs
% all_midtimes = []; % m edian
WNonset = nan(length(OUTSTRUCT.bnum),1);
for i=1:length(indsgrpU)
    disp(i);
    
    % === what are target syls?
    indstmp = indsgrp==indsgrpU(i) & OUTSTRUCT.istarg==1;
    bnum = unique(OUTSTRUCT.bnum(indstmp));
    enum = unique(OUTSTRUCT.enum(indstmp));
    sw = unique(OUTSTRUCT.switch(indstmp));
    
    motiflist = unique(OUTSTRUCT.motifnum(indstmp));
    bname = SwitchStruct.bird(bnum).birdname;
    ename = SwitchStruct.bird(bnum).exptnum(enum).exptname;
    
    % ===
    all_pctiles = [];
    for mm=1:length(motiflist)
        dat = SwitchCohStruct.bird(bnum).exptnum(enum).switchlist(sw).motifnum(motiflist(mm));
        
        if onlyusegoodtargsyls==1
            sylbad = lt_neural_QUICK_LearnRemoveBadSyl(bname, ename, sw, dat.motifname);
            if sylbad==1
                continue
            end
        end
        
        wntimes = dat.WNhittimes_min - PARAMS.motif_predur;
        if useallwn==1
            indsgood = [dat.indsbase_epoch dat.indsWN]; 
        else
            indsgood = [dat.indsbase_epoch dat.indsWN_epoch];
        end
        
        % === only look at WN times within baseline or WN epochs
        wntimes = wntimes(indsgood);
        wntimes = wntimes(~isnan(wntimes)); % only look at hit trials
        
        % ============ GET PERCENTILES FOR TIMES
        all_pctiles = [all_pctiles; prctile(wntimes, prctiletouse)];
    end
    if any(isnan(all_pctiles))
        keyboard
    end
    
    % ================ PUT THIS BACK INTO OUTSTRUCT
    WNonset(indsgrp==indsgrpU(i)) = mean(all_pctiles);
    allwntimes = [allwntimes; mean(all_pctiles)];
end
OUTSTRUCT.WNonset = WNonset;


