function [OUTSTRUCT, PARAMS] = lt_neural_Coher_RealignbyWN(OUTSTRUCT, SwitchCohStruct, SwitchStruct, PARAMS)

%%
onlyusegoodtargsyls = 1;
prctiletouse = 2.5;
[indsgrp, indsgrpU] = lt_tools_grp2idx({OUTSTRUCT.bnum, OUTSTRUCT.enum, OUTSTRUCT.switch});

%% ======== GET TIMING OF WN ONSETS
% ======== for each switch, figure out timing for this WN file
% --- if multiple targets, then will take the mean of value for each
% target.

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
        indsgood = [dat.indsbase_epoch dat.indsWN_epoch];
        
        % === only look at WN times within baseline or WN epochs
        wntimes = wntimes(indsgood);
        wntimes = wntimes(~isnan(wntimes)); % only look at hit trials
        
        % ============ GET PERCENTILES FOR TIMES
        all_pctiles = [all_pctiles; prctile(wntimes, prctiletouse)];
    end
    
    
    % ================ PUT THIS BACK INTO OUTSTRUCT
    WNonset(indsgrp==indsgrpU(i)) = mean(all_pctiles);
end
OUTSTRUCT.WNonset = WNonset;



%% =========== FOR EACH CASE, SHIFT COHGRAMS SO THAT ALIGNED TO WN ONSET

fieldtoget = 'CohMean_Base';
 [ytmp, xtmp] = lt_neural_Coher_RealignbyWN_sub1(OUTSTRUCT, PARAMS, fieldtoget);
% -- add back to OUTSTRUCT.
OUTSTRUCT.([fieldtoget '_orig']) = OUTSTRUCT.(fieldtoget);
OUTSTRUCT.(fieldtoget) = ytmp;

fieldtoget = 'CohMean_WN';
 [ytmp, xtmp] = lt_neural_Coher_RealignbyWN_sub1(OUTSTRUCT, PARAMS, fieldtoget);
% -- add back to OUTSTRUCT.
OUTSTRUCT.([fieldtoget '_orig']) = OUTSTRUCT.(fieldtoget);
OUTSTRUCT.(fieldtoget) = ytmp;



PARAMS.tbins_BeforeAlignWN = PARAMS.tbins;
PARAMS.tbins = xtmp;
PARAMS.didRealignToWN=1;
