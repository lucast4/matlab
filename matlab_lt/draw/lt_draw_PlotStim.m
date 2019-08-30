function [allplots, allplots_mask] = lt_draw_PlotStim(plotEasy, ...
    studytime, nStim, nObj, allplots, allplots_mask, genMode)
% genMode use this to generate stimuli, will close all plots.
if ~exist('genMode', 'var')
    genMode = 0;
end

if genMode==1
    studytime=0;
    dodraw =0;
elseif genMode==0
    assert(exist('studytime', 'var')==1);
    dodraw =1;
end

%% if have allplots, then will run using old stimuli.
% the other params will still matter

%% plots sequence of stimuli.
% saves output.

% plotEasy=1;
% studytime = 2; % 0.8

%% optional - input sequence of stimuli to plot

if ~exist('allplots', 'var')
    newstim =1;
elseif ~isempty(allplots)
    newstim = 0;
    assert(exist('allplots_mask','var')==1);
    nStim = length(allplots);
else
    newstim =1;
end

%% ============================ PLOT STIM
if newstim==1
    allplots = {};
    allplots_mask = {};
end

for j=1:nStim
    disp(['Stim #' num2str(j)]);
    figure; hold on;
    N = nObj;
    
    %     % ==========
    %     h1 = subplot(2,1,1); hold on;
    %     set(gcf, 'Color', 'k')
    %     set(gca, 'Color', 'k');
    %     set(gca, 'Visible', 'off');
    %     set(gca, 'Xtick', []);
    %     xlim([0 10]);
    %     ylim([0 10]);
    
    h2 = subplot(2,1,2); hold on;
    title(num2str(j), 'Color', 'w');
    set(gcf, 'Color', 'k')
    set(gca, 'Color', 'k');
    set(gca, 'Visible', 'off');
    set(gca, 'Xtick', []);
    xlim([0 10]);
    ylim([0 10]);
    
    
    if newstim==1
        htmp = {};
        for i=1:N
            c = [3*rand+4 3*rand+4];
            w = 2;
            h = 2;
            iscirc = randi(3)-1;
            
            if plotEasy==1
                subplot(h2);
                edgeon = 1;
                hfig = lt_draw_insert(c, w, h, iscirc, edgeon);
            else
                subplot(h2)
                edgeon = 0;
                hfig = lt_draw_insert(c, w, h, iscirc, edgeon);
            end
            htmp = [htmp; hfig];
        end
        allplots = [allplots; {htmp}];
    elseif newstim==0
        % ===== PLOT PREVIOUS PLOT
        N = length(allplots{j});
        for i=1:N
            tmp = allplots{j}{i};
            
            if plotEasy==1
                subplot(h2);
                edgeon = 1;
                lt_draw_insert([], [], [], [],edgeon, tmp);
            else
                subplot(h2)
                edgeon = 0;
                lt_draw_insert([], [], [], [],edgeon, tmp);
            end
        end
    end
    pause(studytime);
    % close all;
    
    
    %% %%%%%%%%%%%%%%%%%%% RANDOM STIMULI TO REMOVE AFTEREFFECTS.
    figure; hold on;
    N = randi(10)+40;
    
    %     h1 = subplot(2,1,1); hold on;
    %     set(gcf, 'Color', 'k')
    %     set(gca, 'Color', 'k');
    %     set(gca, 'Visible', 'off');
    %     set(gca, 'Xtick', []);
    %     xlim([0 10]);
    %     ylim([0 10]);
    
    h2 = subplot(2,1,2); hold on;
    set(gcf, 'Color', 'k')
    set(gca, 'Color', 'k');
    set(gca, 'Visible', 'off');
    set(gca, 'Xtick', []);
    xlim([0 10]);
    ylim([0 10]);
    
    if newstim==1
        htmp = {};
        for i=1:N
            c = [6*rand+2 6*rand+2];
            w = 0.5;
            h = 0.5;
            iscirc = 1;
            
            subplot(h2);
            edgeon = 1;
            hfig = lt_draw_insert(c, w, h, iscirc, edgeon);
            htmp = [htmp; hfig];
        end
        allplots_mask = [allplots_mask; {htmp}];
    else
        
        % ===== PLOT PREVIOUS PLOT
        N = length(allplots_mask{j});
        for i=1:N
            tmp = allplots_mask{j}{i};
            
            subplot(h2);
            edgeon = 1;
            lt_draw_insert([], [], [], [],edgeon, tmp);
        end
        
    end
    if dodraw==1
        pause;
    end
    
end
