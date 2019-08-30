fig = figure; hold on;
set(gcf, 'Color', 'k')
set(gca, 'Color', 'k');
set(gca, 'Visible', 'off');
set(gca, 'Xtick', []);

box off

% rectangle('Position', [1 1 1 1], 'FaceColor', 'w');
% rectangle('Position', [1 1 1 1], 'Curvature', 1, 'FaceColor', 'w');
%
% rectangle('Position', [2 2 1 1], 'Curvature', 0, 'FaceColor', 'w');
% rectangle('Position', [2.5 2.5 1 1], 'Curvature', 1, 'FaceColor', 'w');

rectangle('Position', [3 3 1 1], 'Curvature', 0, 'FaceColor', 'w', 'EdgeColor', 'w');
rectangle('Position', [3.5 3.5 1 1], 'Curvature', 1, 'FaceColor', 'w', 'EdgeColor', 'w');

rectangle('Position', [5 5 1 1], 'Curvature', 0, 'FaceColor', 'w', 'EdgeColor', 'k');
rectangle('Position', [5.5 5.5 1 1], 'Curvature', 1, 'FaceColor', 'w', 'EdgeColor', 'k');

xlim([0 10]);
ylim([0 10]);


%%

fig = figure; hold on;
set(gcf, 'Color', 'k')
set(gca, 'Color', 'k');
set(gca, 'Visible', 'off');
set(gca, 'Xtick', []);
xlim([0 10]);
ylim([0 10]);

box off

c = [8 1];
w = 1;
h = 1;
iscirc = 0;
edgeon = 0;

lt_draw_insert(c, w, h, iscirc, edgeon)
lt_draw_insert(c+0.5, w, h, iscirc, edgeon)

%% =============== DRAW A RANDOM SET OF OBJECTS
close all;

plotEasy=1;
studytime = 1; % 0.8
nStim = 5;
nObj = 5; % max num obj
[allplots, allplots_mask] = lt_draw_PlotStim(plotEasy, ...
    studytime, nStim, nObj);

disp('NEXT BLOCK!!!');

plotEasy=0;
lt_draw_PlotStim(plotEasy, ...
    studytime, nStim, nObj, allplots, allplots_mask);

%% =============== [FIRST GENERATE STIMULI]

close all;

plotEasy=1;
studytime = 0.9; % 0.8
nStim = 10;
nObj = 5; % max num obj
genMode =1;
[allplots, allplots_mask] = lt_draw_PlotStim(plotEasy, ...
    studytime, nStim, nObj, [], [], genMode);


% ============== PLOT EACH STIM 2X (easy and hard), randomized
close all;
easyVec = [zeros(1, length(allplots)), ones(1, length(allplots))];
trialInds = [1:length(allplots) 1:length(allplots)];

indrand = randperm(2*length(allplots));
easyVec = easyVec(indrand);
trialInds = trialInds(indrand);


for i=1:length(trialInds)
    indthis = trialInds(i);
    plotEasy=easyVec(i);

% lt_draw_PlotStim(plotEasy, ...
%     studytime, [], [], allplots(indthis), allplots_mask(indthis));
lt_draw_PlotStim(plotEasy, ...
    studytime, [], [], allplots(indthis), allplots_mask(indthis));
end

%%