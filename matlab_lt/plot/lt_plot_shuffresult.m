function lt_plot_shuffresult(ydat, yshuff)
%% lt 10/3/18 - plots histogram of shuffle result, plus data, and p val
% ydat scalar
% yshuff vector

nshuff = length(yshuff);
%%
lt_plot_histogram(yshuff);
line([ydat ydat], ylim, 'Color', 'r');
p = (sum(yshuff>=ydat)+1)/(nshuff+1);
lt_plot_pvalue(p, 'shuff>dat', 1);

