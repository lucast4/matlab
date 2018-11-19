%% lt 11/6/18 - to simulate coherence calculations
function [coh, powercorr, power_a, power_b] = lt_tools_cohsandbox(phifrac, n, plotraw, mixcoeff)
% phifrac = 1, then entire uniform distribution,... [i.e fraction of entire
% 2*pi); range: 0 to 1

% n = trials;
% n = 100;

% mixcoeff = 1, then a and b are identical, 0, then entirely independnet.\
if ~exist('mixcoeff', 'var')
    mixcoeff = 0;
end
%%
a = rand(n,1);
b = rand(n,1);
c = rand(n,1);
% -- mix
a = (1-mixcoeff)*a + mixcoeff*c;
b = (1-mixcoeff)*b + mixcoeff*c;
% - power correlation
powercorr = corr(a,b);
power_a = sum(a.^2);
power_b = sum(b.^2);

phi = phifrac*(2*pi*rand(n,1));

if plotraw==1
figure; hold on;
subplot(2,2,1); hold on;
title('a and b');
plot(a, '-or');
plot(b, '-ob');

subplot(2,2,2); hold on;
rose(phi);
title('phi');
end

% === get coherence
coh = abs(sum((a.*b).*exp(i*phi),1)./n);

end