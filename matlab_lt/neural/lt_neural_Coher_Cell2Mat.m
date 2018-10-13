function cohmat = lt_neural_Coher_Cell2Mat(CohCell)

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
ntrials = length(CohCell);
cohexampl = CohCell{1};
cohmat = nan(size(cohexampl,1), size(cohexampl,2), ntrials);
for tt=1:ntrials
    if all(size(cohmat(:,:,tt))==size(CohCell{tt}))
        cohmat(:,:,tt) = CohCell{tt};
    else
        disp('skipped! wrong size');
    end
end


