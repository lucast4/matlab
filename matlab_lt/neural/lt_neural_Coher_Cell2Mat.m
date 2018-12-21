function cohmat = lt_neural_Coher_Cell2Mat(CohCell)
% NOTE: tkaes min dimensions across all cells... (doens't matter if they
% are all same size)

assert(any(size(CohCell)==1), 'problem, can only have one channel ...');
% each cell is time x ff. multiple cells, one for each trial
% e.g. CohCell =
%
%   2527×1 cell array
%
%     [20×19 double]
%     [20×19 double]
%     [20×19 double]
%     [20×19 double]
%     [20×19 double]
%     [20×19 double]
%     [20×19 double]

%%

if isempty(CohCell)
    cohmat = [];
    return
end

mint = min(cellfun(@(x)size(x,1), CohCell, 'UniformOutput', 1));
minf = min(cellfun(@(x)size(x,2), CohCell, 'UniformOutput', 1));
ntrials = length(CohCell);
cohmat = nan(mint, minf, ntrials);

for tt=1:ntrials
    cohmat(:,:,tt) = CohCell{tt}(1:mint, 1:minf);
    
%     if all(size(cohmat(:,:,tt))==size(CohCell{tt}))
%         cohmat(:,:,tt) = CohCell{tt};
%     else
%         disp('skipped! wrong size');
%     end
end


