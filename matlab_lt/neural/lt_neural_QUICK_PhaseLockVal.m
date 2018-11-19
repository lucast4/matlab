function [plv, mean_cart, mean_pol, bootstats] = lt_neural_QUICK_PhaseLockVal(theta, doboot)
%% lt 11/12/18 - Operates on vector of angles.
% INPUT:
%  theta is [-pi, pi], as would get from phi in coherence
%   doboot = 1; then gets bootstrap stats
% calculations;
% NOTE: can be any angle from -inf to inf
% NOTE: theta of 0 is x,y = (1,0);
% NOTE: goes counterclockwise.

% OUTPUT:
%  plv is norm of sum of all angles (assume unity)
% mean_pol is angle of that mean vector.
% mean_cart is mean vector, in [x,y];

nboot = 15; % for getting standard deviatio...
%%

if ~exist('doboot', 'var')
    doboot =0;
end

%%

% ==== DATA
[x,y] = pol2cart(theta, ones(size(theta))); % assume unit vectors.
mean_cart = [mean(x) mean(y)];
[mean_pol, plv] = cart2pol(mean_cart(1), mean_cart(2));

% ==== BOOSTRAP
bootstats = struct;
if doboot==1
    polallshuff =[];
    plvallshuff = [];
    for n=1:nboot
        % --- get bootstrap CI for mean angle and plv
        indshuff = randi(length(theta), length(theta), 1);
        theta_shuff = theta(indshuff);
        
        % --- rerun using Resampled
        [x,y] = pol2cart(theta_shuff, ones(size(theta_shuff))); % assume unit vectors.
        [mean_pol_shuff, plv_shuff] = cart2pol(mean(x), mean(y));
        
        polallshuff =[polallshuff; mean_pol_shuff];
        plvallshuff = [plvallshuff; plv_shuff];
        
    end
    
    bootstats.plv.std = std(plvallshuff);
end

