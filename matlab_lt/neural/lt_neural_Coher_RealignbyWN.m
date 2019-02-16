function [OUTSTRUCT, PARAMS] = lt_neural_Coher_RealignbyWN(OUTSTRUCT, ...
    SwitchCohStruct, SwitchStruct, PARAMS, useallwn)

%%
onlyusegoodtargsyls = 1;
prctiletouse = 2.5;

%% ======== GET TIMING OF WN ONSETS
% ======== for each switch, figure out timing for this WN file
% --- if multiple targets, then will take the mean of value for each
% target.

[OUTSTRUCT, allwntimes] = lt_neural_Coher_GetWNtiming(OUTSTRUCT, onlyusegoodtargsyls, ...
    useallwn, prctiletouse, SwitchStruct, SwitchCohStruct, PARAMS);

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
