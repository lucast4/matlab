function lt_tools_lme_plotsummary(lme)
%% lt 6/2/19 - various summary plots

lt_figure; hold on;

%% ======= fixed effects
lt_subplot(3,4,1); hold on;
lt_tools_lme_plotEffects(lme, 0);

%% =============== PLOT RANDOM EFFECTS
[B, Bnames, Stats] = randomEffects(lme);

lt_subplot(3,4,2:4); hold on;
ylabel('random effects');
Xnames = {};
for x=1:length(Bnames.Name)
    y = B(x);
    yCI = [Stats.Lower(x) Stats.Upper(x)];
    p = Stats.pValue(x);
    
    errorbar(x, y,  y-yCI(1), -y+yCI(2), 'LineStyle', 'none', 'Color', 'k');
    plot(x, y, 'ok');
    
    % --- pvale
    if p<0.1
        lt_plot_text(x, y, ['p=' num2str(p)], 'r', 8);
    end
    % --- collect xlabels
    Xnames = [Xnames [Bnames.Name{x} '-' Bnames.Group{x} '-' Bnames.Level{x}]];
end
set(gca, 'XTick', 1:length(Xnames), 'XTickLabel', Xnames);

lt_plot_zeroline;
rotateXLabels(gca, 45);

%% ======= residuals
lt_subplot(3,4,5); hold on;
plotResiduals(lme, 'fitted');

