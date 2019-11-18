function lt_plot_distributionPlot(Y, ver)
if ~exist('ver', 'var')
    ver=2;
end


x = 1:size(Y,2);

if ver==1
    distributionPlot(Y,'histOpt',0, 'showMM', 6);
elseif ver==2
    binsize = 0.025;
    tmp = round(100*max(abs(Y(:))))/100;
    tmp = tmp+(binsize-mod(tmp, binsize));
    % tmp = tmp+binsize/2;
    xcenter = -tmp:binsize:tmp;

    distributionPlot(Y,'xValues', x, 'showMM', 4, 'addSpread', 0, 'color', ...
        'b', 'histOpt', 0, 'divFactor', xcenter);
end