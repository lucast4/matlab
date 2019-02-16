function     lt_neural_PLOT_rasterline(spktimes, yval, plotcol, plotdot, ...
    linesize)

if ~exist('plotdot', 'var')
    plotdot = [];
end

if ~exist('linesize', 'var')
    linesize = 0.8;
end

%% plots a single line of raster (lt, 11/4/17)
% spktimes; array, in sec
% yval: position to plot (will go from -1 to ...)


if ~isempty(spktimes)
if plotdot==1
   plot(spktimes, yval, '.', 'Color', plotcol); 
else
    for ttt =1:length(spktimes)
        
        line([spktimes(ttt) spktimes(ttt)], [yval-linesize/2 yval+linesize/2], ...
            'Color', plotcol, 'LineWidth', 1);
    end
end
end