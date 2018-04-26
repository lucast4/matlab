function x = lt_neural_QUICK_XfromFRmat(frmat)

%% lt 4/24/18 - give timebase (x) from a frmat (time x trials), 
% assumes 1ms timebase

%% 
% starts from 0

n = size(frmat, 1);
x = 0:n-1;
x = 0.001*x;
