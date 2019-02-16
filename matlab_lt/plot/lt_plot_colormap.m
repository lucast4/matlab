function lt_plot_colormap(type)
%% lt 1/3/19 - for heatmaps
% type = 'centered'; % then centered at 0
% type = 'pval'; % white to black

n = 25; % number of gradations,.

%%
if strcmp(type, 'centered')
    %         cmap = [[linspace(1, 0, 10)' linspace(0, 0, 10)' linspace(0, 0, 10)']; ...
    %     [linspace(0, 0, 10)' linspace(0, 0, 10)' linspace(0, 1, 10)']];
    %
    %     n = 10;
    %     cmap = [[linspace(0, 1, n)' linspace(1, 1, n)' linspace(1,1, n)']; ...
    %     [linspace(1,1, n)' linspace(1, 1, n)' linspace(1, 0, n)']];
    
    
    cmap = [0.8*[linspace(1, 1, n)' linspace(0.2, 1, n)'  linspace(0.2, 1, n)']; ...
        0.8*[linspace(1, 0.2, n)' linspace(1, 0.2, n)' linspace(1, 1, n)' ]];
    
    cmap = flipud(cmap);
elseif strcmp(type, 'pval')
    
    cmap = [linspace(0, 1, n)' linspace(0, 1, n)'  linspace(0, 1, n)'];
end

colormap(gca, cmap);
