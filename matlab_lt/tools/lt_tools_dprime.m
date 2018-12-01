function dprime = lt_tools_dprime(x1, x2)
%% lt 11/27/18 - gets dprime of x1 minus x2.

dprime = (mean(x1)-mean(x2))/(sqrt(0.5*(var(x1)+var(x2))));
