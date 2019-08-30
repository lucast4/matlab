function hfig = lt_draw_insert(c, w, h, iscirc, edgeon, handle)
%% lt 5/29/19 - inserts a new drawing
% NOTE: hfig is a struct, not a handle. did this so would not empty when
% close figure.

% === option 1: enter params
% c = [5 4]; %center (x y)
% w = 3;
% h = 4;
%
% iscirc = 1; 0=rect, 2=triangle
% edgeon = 1;

% === option 2: enter handle (rectangle object)
% first 4 params must be empty.
% onluy that matter are: edgeon and handle.
% leave edgeon empty to use whatever is in handle.


%%
if edgeon==1
    ec = 'k';
else
    ec = 'w';
end
if isempty(c)
    assert(exist('handle', 'var')==1);
    
    if strcmp(handle.Type, 'rectangle')
    hfig = rectangle('Position', handle.Position, 'Curvature', ...
        handle.Curvature, 'FaceColor', handle.FaceColor, ...
        'EdgeColor', ec, 'LineWidth', ...
        handle.LineWidth, 'LineStyle', handle.LineStyle);
    elseif strcmp(handle.Type, 'patch')
        hfig = fill(handle.Vertices(:,1), handle.Vertices(:,2), ...
            handle.FaceColor, 'EdgeColor', ec);
    end
else
    if iscirc==2
        % ==== draw a triangle
    hfig = fill([c(1)-w/2 c(1)+w/2 c(1)], [c(2)-h/2 c(2)-h/2 c(2)+h/2], ...
        'w', 'EdgeColor', ec);
        
%     hfig = rectangle('Position', [c(1)-w/2 c(2)-h/2 w h], 'Curvature', [iscirc iscirc], ...
%         'FaceColor', 'w', 'EdgeColor', ec);    
    else
    
    hfig = rectangle('Position', [c(1)-w/2 c(2)-h/2 w h], 'Curvature', [iscirc iscirc], ...
        'FaceColor', 'w', 'EdgeColor', ec);
    end
end

% === extract fields from fig, so is not handle anymore.
tmp = struct;
fnames = fieldnames(hfig);
for j=1:length(fnames)
    tmp.(fnames{j}) = hfig.(fnames{j});
end
hfig = tmp;