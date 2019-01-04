function lt_plot_colormap(type)
%% lt 1/3/19 - for heatmaps
% type = 'centered'; % then centered at 0



%%
if strcmp(type, 'centered')
%         cmap = [[linspace(1, 0, 10)' linspace(0, 0, 10)' linspace(0, 0, 10)']; ...
%     [linspace(0, 0, 10)' linspace(0, 0, 10)' linspace(0, 1, 10)']];
% 
%     n = 10;
%     cmap = [[linspace(0, 1, n)' linspace(1, 1, n)' linspace(1,1, n)']; ...
%     [linspace(1,1, n)' linspace(1, 1, n)' linspace(1, 0, n)']];

    n = 25;
    cmap = [0.8*[linspace(1, 1, n)' linspace(0.2, 1, n)'  linspace(0.2, 1, n)']; ...
    0.8*[linspace(1, 0.2, n)' linspace(1, 0.2, n)' linspace(1, 1, n)' ]];

colormap(cmap);
end
