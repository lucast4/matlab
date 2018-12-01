function lt_neural_QUICK_PlotSylPatches(ons, offs, y, linemode)
%% lt 11/21/18 - like raster, but patches for syls


%        ons =
%
%   Columns 1 through 4
%
%                          0         0.117333333333335         0.297166666666662         0.489533333333334
%
%   Columns 5 through 8
%
%          0.573299999999996         0.735099999999996          1.07296666666667          1.20693333333333
%
%   Columns 9 through 12
%
%           1.36206666666666          1.51903333333333          1.66276666666666                    1.7418




%           offs =
%
%   Columns 1 through 4
%
%           0.03393333333333         0.194099999999999         0.354700000000001         0.552799999999998
%
%   Columns 5 through 8
%
%          0.656733333333335         0.939966666666663                    1.1325          1.27713333333333
%
%   Columns 9 through 12
%
%                     1.4266                    1.5873                    1.7214                    1.8084

% y = 1 % line to plot.

% linemode=1; then instead of patches plots lines.

if ~exist('linemode', 'var')
    linemode=0;
end
%% % RUN
for kk=1:length(ons)
    
    if linemode==1
        plot([ons(kk), offs(kk)], [y y], '-k', 'LineWidth', 2);
    else
    X = [ons(kk) offs(kk) offs(kk) ons(kk)];
    Y = [y-0.4 y-0.4 y+0.4 y+0.4];
    h=patch(X, Y, 'w', 'FaceColor', [0 0 0], 'EdgeColor', [0.8 0.8 0.8], ...
        'FaceAlpha', 0.1);
    set(h, 'FaceAlpha', 0.1);
    end
end
