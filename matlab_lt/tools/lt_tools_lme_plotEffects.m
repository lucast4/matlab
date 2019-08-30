function lt_tools_lme_plotEffects(lme, newplot)
%% lt 5/26/19 - plot fixed effects
if ~exist('newplot', 'var')
    newplot = 1;
end

%% ======================
fixednames = lme.CoefficientNames;
fixedeffCI = lme.coefCI;
fixedeff = lme.fixedEffects;
fixedeffSE = lme.Coefficients.SE;

if newplot==1
lt_figure; hold on;
end
x = 1:length(fixednames);
errorbar(x, fixedeff,  fixedeff-fixedeffCI(:,1), -fixedeff+fixedeffCI(:,2), 'LineStyle', 'none');
errorbar(x, fixedeff,  fixedeff-fixedeffCI(:,1), -fixedeff+fixedeffCI(:,2), 'LineStyle', 'none');
% --- also plot SEM
errorbar(x, fixedeff,  fixedeffSE, 'LineStyle', 'none', 'Color', 'r');
plot(x, fixedeff, 'ok');
xlim([0 x(end)+1]);
lt_plot_zeroline;
set(gca, 'XTick', x, 'XTickLabel', fixednames);
rotateXLabels(gca, 45);

% --------- plot p values
YLIM = ylim;
pvals = lme.Coefficients.pValue;
for x=1:length(pvals)
    if pvals(x)<0.1
   lt_plot_text(x-0.3, YLIM(2), ['p=' num2str(pvals(x))], 'r', 10)
    end
end

