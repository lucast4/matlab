function [datsort, ind1, ind2] = lt_tools_sort(dat)
%% lt 3/25/19 - same as matlab sort, except also outputs the rank of each element

% e.g.
% dat = [1 10 3 4 0];
% [datsort, ind1, ind2] = lt_tools_sort(dat);
% 
% then 
% datsort = dat(ind1),
% and ind2 is the rank of each element (i.e. [2 5 3 4 1])
% i.e. if want to know, for each index, what where it ranks (1 means
% smallest ...);


[datsort, ind1] = sort(dat);

ind2 = nan(size(ind1));
ind2(ind1) = 1:length(ind1);

