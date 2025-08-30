function fh = show_flow_field(data, opt)
% function fh = show_flow_field(data, opt)
% Tstamp = 1 x time
% raster = cell x time x trial
% ID = cell x trial

def.min_sample = 100;
def.xlim = [];
def.ylim = [];
def.zlim = [];
def.trj_wid = 4;
def.trj_color = [1 0 0; 0 0 1];
def.input_col = [0 0 0];
def.legend = {};
def.xlabel = 'X';
def.ylabel = 'Y';

opt = safeStructAssign(def, opt);

U = data.U;
V = data.V;
Cnt = data.Cnt;
Xgrid = data.Xgrid;
Ygrid = data.Ygrid;
if isfield(data, 'mean_trajectory')
    mean_trajectory = data.mean_trajectory;
end
if isfield(data, 'W')
    InputWeight = data.W;
end

ncond = size(mean_trajectory, 2);

U(Cnt(:) < opt.min_sample) = 0;
V(Cnt(:) < opt.min_sample) = 0;

nonZeroMask = (U ~= 0) | (V ~= 0);
[rowIdx, colIdx] = find(nonZeroMask);

rowMin = min(rowIdx);
rowMax = max(rowIdx);
colMin = min(colIdx);
colMax = max(colIdx);

Xgrid = Xgrid(rowMin:rowMax, colMin:colMax);
Ygrid = Ygrid(rowMin:rowMax, colMin:colMax);
U = U(rowMin:rowMax, colMin:colMax);
V = V(rowMin:rowMax, colMin:colMax);

% Plot
fh = figure;
hold on;
quiver(Xgrid, Ygrid, U, V, 'AutoScale', 'on', 'LineWidth', 1.5, 'color', 'k');

ch = nan(ncond,1);
if isfield(data, 'mean_trajectory')
    for c = 1:ncond
        ch(c) = plot(mean_trajectory(:, c, 1), mean_trajectory(:, c, 2), 'color', opt.trj_color(c, :), 'linestyle', '-', 'linew', opt.trj_wid);
    end
end
if ~isempty(opt.legend)
    legend(ch, opt.legend, 'autoupdate', 'off');
end
ax = gca;
xlims = ax.XLim;
ylims = ax.YLim;

if isfield(data, 'W')
    InputWeight = InputWeight / max(abs(InputWeight(:))) * (xlims(2) - xlims(1)) / 7;
    hw = (xlims(2) - xlims(1)) / 40;
    for c= 1:size(InputWeight, 2)
        drawAxisArrow(gca, 0, 0, InputWeight(1, c), InputWeight(2, c), 'LineWidth', 4, 'HeadWidth', hw, 'HeadLength', hw, 'Color', opt.input_col(c,:));
        drawAxisArrow(gca, 0, 0, -InputWeight(1, c), -InputWeight(2, c), 'LineWidth', 4, 'HeadWidth', hw, 'HeadLength', hw, 'Color', opt.input_col(c,:));
    end
end
set(gca, 'xlim', xlims, 'ylim', ylims');
xlabel(opt.xlabel);
ylabel(opt.ylabel);
axis equal;

end


function h = drawAxisArrow(ax, x_start, y_start, x_end, y_end, varargin)
% DRAWAXISARROW Draws an arrow in axis coordinates with triangular head
%
%   h = drawAxisArrow(ax, x_start, y_start, x_end, y_end)
%   h = drawAxisArrow(ax, x_start, y_start, x_end, y_end, 'PropertyName', PropertyValue, ...)
%
%   Inputs:
%       ax        - Target axes handle
%       x_start   - Starting x-coordinate
%       y_start   - Starting y-coordinate
%       x_end     - Ending x-coordinate
%       y_end     - Ending y-coordinate
%
%   Optional Properties:
%       'Color'      - Arrow color (default: 'k' black)
%       'LineWidth'  - Shaft width (default: 2)
%       'HeadWidth'  - Width of arrowhead (default: 0.5, in axis units)
%       'HeadLength' - Length of arrowhead (default: 0.5, in axis units)
%       'LineStyle'  - Line style (default: '-')
%
%   Output:
%       h - Handle to the patch object containing the arrow

    % Default parameters
    params.Color = 'k';
    params.LineWidth = 2;
    params.HeadWidth = 0.5;
    params.HeadLength = 0.5;
    params.LineStyle = '-';
    
    % Parse optional parameters
    if ~isempty(varargin)
        for i = 1:2:length(varargin)
            if isfield(params, varargin{i})
                params.(varargin{i}) = varargin{i+1};
            end
        end
    end
    
    % Calculate arrow direction vector
    dx = x_end - x_start;
    dy = y_end - y_start;
    arrow_length = sqrt(dx^2 + dy^2);
    
    % Normalize direction vector
    if arrow_length > 0
        ux = dx / arrow_length;
        uy = dy / arrow_length;
    else
        ux = 0;
        uy = 1; % Default vertical if zero length
    end
    
    % Calculate perpendicular vector (for head width)
    px = -uy;
    py = ux;
    
    % Adjust head dimensions if they're too large for the arrow
    if params.HeadLength > 0.8 * arrow_length
        params.HeadLength = 0.8 * arrow_length;
    end
    
    % Calculate shaft end point (start of arrowhead)
    shaft_end_x = x_end - ux * params.HeadLength;
    shaft_end_y = y_end - uy * params.HeadLength;
    
    % Points for arrow shaft (rectangle)
    shaft_width = params.LineWidth/70; % Scale line width to axis units
    shaft_points_x = [
        x_start + px * shaft_width
        shaft_end_x + px * shaft_width
        shaft_end_x - px * shaft_width
        x_start - px * shaft_width
    ];
    
    shaft_points_y = [
        y_start + py * shaft_width
        shaft_end_y + py * shaft_width
        shaft_end_y - py * shaft_width
        y_start - py * shaft_width
    ];
    
    % Points for arrow head (triangle)
    head_points_x = [
        shaft_end_x + px * params.HeadWidth/2
        x_end
        shaft_end_x - px * params.HeadWidth/2
    ];
    
    head_points_y = [
        shaft_end_y + py * params.HeadWidth/2
        y_end
        shaft_end_y - py * params.HeadWidth/2
    ];
    
    % Combine all points
    all_x = [shaft_points_x(1:2); head_points_x; shaft_points_x(3:4)];
    all_y = [shaft_points_y(1:2); head_points_y; shaft_points_y(3:4)];
    
    % Draw the arrow as a single patch object
    h = patch(ax, all_x, all_y, params.Color, ...
              'EdgeColor', params.Color, ...
              'LineWidth', params.LineWidth, ...
              'LineStyle', params.LineStyle);
end


