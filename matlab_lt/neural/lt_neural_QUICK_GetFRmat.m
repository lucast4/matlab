%% extracts fr mat that is clipped (predur and postdur relative to alignment)
function [frmat, t] = lt_neural_QUICK_GetFRmat(frmat, t, motif_predur, pretime, posttime)

% INPUTS:
% frmat, tbins x trials
% t, timebins for frmat
% motif_predur, sec, time of alignment
% pretime, how much time to keep preceding alignemnt, sec
% posttime, same, but for post alignment

% OUTPUT:
% frmat, same as input, but clipped in timebase.
% t, new timebase, changed so now with 0 at alginment point

%%
% --- which inds
twind = motif_predur + [-pretime posttime];
indskeep = t>=twind(1) & t<twind(2);

% --- get fr mat
frmat = frmat(indskeep, :);

% --- get, and center, timebase
t = t(indskeep) - motif_predur;