function [ccRealAll, ccShiftAll, xlags_sec, neurpairID] = lt_neural_POPLEARN_FlipXcovDat(XcovDat, bregionwanted)
%% ask for brain region pair, extracts that pair (and flips so that is in order
% of the input (e.g. RA-LMAN is flipped to LMAN-RA)

% e.g., input is data from MOTIFSTATS_pop:
% XcovDat = MOTIFSTATS_pop.birds(i).exptnum(ii).DAT.setnum(ss).motif(mm).XCov_neurpair;
% bregionwanted = {'LMAN', 'RA'}; %

% % === output
% neurpairID, index after flatten XcovDat
%% ====== find the neuron pairs that match desired brain
% regions - note if need to fliplr the xcov
bregionlist = {XcovDat(:).bregtmp};
indspairs = [];
flippair = [];

for bb = 1:length(bregionlist)
    if all(strcmp(bregionlist{bb}, bregionwanted))
        % -- keep and don't flip
        indspairs = [indspairs bb];
        flippair = [flippair 0];
    elseif all(strcmp(fliplr(bregionlist{bb}), bregionwanted))
        % -- keep, and flip
        indspairs = [indspairs bb];
        flippair = [flippair 1];
    end
end
flippair = logical(flippair);

if isempty(indspairs)
    % then no neuron pairs with these brain regions
    ccRealAll = [];
    ccShiftAll = [];
    xlags_sec = [];
    return
end

% =================================== COLLECT XCOV TRACES
ccRealAll = {XcovDat(indspairs).ccRealAll};
ccShiftAll = {XcovDat(indspairs).ccShiftAll};
xlags_sec = XcovDat(indspairs(1)).x;

% -- flip any if needed
ccRealAll(flippair) = cellfun(@fliplr, ccRealAll(flippair), 'UniformOutput', 0);
ccShiftAll(flippair) = cellfun(@fliplr, ccShiftAll(flippair), 'UniformOutput', 0);

%%

neurpairID = indspairs;
