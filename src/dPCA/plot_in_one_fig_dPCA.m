function [fh_1, fh_2] = plot_in_one_fig_dPCA(fh_idv, layout, figsize, task_ylim, choice_ylim)

if nargin<4
    choice_ylim = [-60 60];
end

if nargin<4
    task_ylim = [-40 40];
end

if nargin<3
    figsize = [300 300];
end

if nargin<2
    layout = [2 ceil(length(fh_idv)/2)];
end

fh_1 = figure('Position', [50 100 figsize(1) figsize(2)]);
fh_2 = figure('Position', [50 100 figsize(1) figsize(2)]);

for i = 1:length(fh_idv)
    axesList = findobj(fh_idv{i}, 'Type', 'axes');
    positions = arrayfun(@(ax) ax.Position(2), axesList);
    [~, idx] = sort(positions, 'descend');
    axesList = axesList(idx);

    figure(fh_1) % task axis
    subplot(layout(1),layout(2),i)
    copyobj(axesList(2).Children, gca);
    hold on; plot(xlim, [0 0], ':', 'Color', [.5 .5 .5])
    title(['session ', num2str(i)])
    format_panel(gca, 'ylim', task_ylim, 'xlabel', 'Time (ms)', 'ylabel', {'Pair signal'})
    xtickangle(0)

    figure(fh_2) % choice axis
    subplot(layout(1),layout(2),i)
    copyobj(axesList(1).Children, gca);
    hold on; plot(xlim, [0 0], ':', 'Color', [.5 .5 .5])
    title(['session ', num2str(i)])
    format_panel(gca, 'ylim', choice_ylim, 'xlabel', 'Time (ms)', 'ylabel', {'Choice signal'})
    xtickangle(0)
end

end