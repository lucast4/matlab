function [rhomean_dat, rhomean_shuff, rhomean_dat_diff, rhomean_shuff_diff] =...
    lt_neural_v2_ANALY_FRsmooth_CompSylTypes3(OUTDAT, SwitchStruct, ...
    MOTIFSTATS_Compiled, nshuff, usediffFromBase, plotdifftype, epochtoplot)

%%
disp('params need to move up.. (I.E. out of hard coding here and use as argumetns)')
pause;
% shuffSylType = 0;
epochtoplot = 3;
timewind = [-0.08 0.03];
% usediffFromBase = 0; % if 0, then also norms to global drift.
syltypesneeded = [1 1 1];

%% FIRST DO ACTUAL DATA (NO SHUFF)

% ======== get distances
plotOn=1;
[All_RhoTargSame, All_RhoTargDiff] = lt_neural_v2_ANALY_FRsmooth_CompSylTypes2(OUTDAT, MOTIFSTATS_Compiled, ...
    SwitchStruct, timewind, usediffFromBase, syltypesneeded, plotOn, plotdifftype, epochtoplot);

rhomean_dat = mean(cellfun(@mean, All_RhoTargSame));
rhomean_dat_diff = mean(cellfun(@mean, All_RhoTargDiff));
% rhomean_dat = mean(cell2mat(All_RhoTargSame));


%% SECOND, DO MULTIPLE SHUFFLES
% nshuffs = 1000;
rhomean_shuff = [];
rhomean_shuff_diff = [];
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
    [All_RhoTargSame, All_RhoTargDiff] = lt_neural_v2_ANALY_FRsmooth_CompSylTypes2(OUTDAT_shuff, MOTIFSTATS_Compiled, ...
        SwitchStruct, timewind, usediffFromBase, syltypesneeded, plotOn, plotdifftype, epochtoplot);
    
    rhomean_shuff = [rhomean_shuff; mean(cellfun(@mean, All_RhoTargSame))];
    rhomean_shuff_diff = [rhomean_shuff_diff; mean(cellfun(@mean, All_RhoTargDiff))];
end