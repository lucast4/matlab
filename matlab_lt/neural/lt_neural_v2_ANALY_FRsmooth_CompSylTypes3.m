function [rhomean_dat, rhomean_shuff] = lt_neural_v2_ANALY_FRsmooth_CompSylTypes3(OUTDAT, SwitchStruct, ...
    MOTIFSTATS_Compiled, nshuff, usediffFromBase)

%%
disp('params need to move up..')
pause;
% shuffSylType = 0;
epochtoplot = 3;
timewind = [-0.08 0.04];
% usediffFromBase = 0; % if 0, then also norms to global drift.
syltypesneeded = [1 1 1];

%% FIRST DO ACTUAL DATA (NO SHUFF)

% ======== get distances
All_RhoTargSame = lt_neural_v2_ANALY_FRsmooth_CompSylTypes2(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, timewind, usediffFromBase, syltypesneeded);

rhomean_dat = mean(cellfun(@mean, All_RhoTargSame));
% rhomean_dat = mean(cell2mat(All_RhoTargSame));


%% SECOND, DO MULTIPLE SHUFFLES
% nshuffs = 1000;
rhomean_shuff = [];
for sh = 1:nshuff
    disp(['shuff no: ' num2str(sh)]);
    % ======== do shuffle
    plotOn = 0;
    shuffonlynontargs =1;
    shuffSylType =1;
    OUTDAT_shuff = lt_neural_v2_ANALY_FRsmooth_Comps(OUTDAT, SwitchStruct, shuffSylType, ...
        epochtoplot, plotOn, shuffonlynontargs);
    
    % ======== get distances
    plotOn =0;
    All_RhoTargSame = lt_neural_v2_ANALY_FRsmooth_CompSylTypes2(OUTDAT_shuff, MOTIFSTATS_Compiled, ...
        SwitchStruct, timewind, usediffFromBase, syltypesneeded, plotOn);
    
    rhomean_shuff = [rhomean_shuff; mean(cellfun(@mean, All_RhoTargSame))];
end