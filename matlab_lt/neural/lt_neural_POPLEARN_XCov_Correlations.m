function lt_neural_POPLEARN_XCov_Correlations(Yscalar, OUTSTRUCT, SwitchCohStruct, SwitchStruct, ...
    MOTIFSTATS_Compiled, PARAMS, clim, onlygoodexpt, expttype, plotlevel)
%% 2/1/19 - Plots each scalar, each experiement...

% clim = [-0.1 0.1];
useglobalmotifname = 0; % then all expt across a bird will be aligned... motif names will not necessarily be correct for each experiement though;
minDiffN = 3; % will only noramlize to diff ytpe if diff type has N this or larger.

%% ======== filter experiments
if onlygoodexpt==1
    % ===== filter outstruct
    [OUTSTRUCT, indstokeep] = lt_neural_Coher_QUICK_FilterOUTSTRUCT(OUTSTRUCT, SwitchStruct, expttype);
    % ===== filter the scalar array
    Yscalar = Yscalar(indstokeep);
end

% ====== SANITY CHECK.
assert(length(Yscalar) == length(OUTSTRUCT.bnum));


