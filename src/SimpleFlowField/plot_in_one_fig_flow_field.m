function fh_all = plot_in_one_fig_flow_field(fh_idv, layout, figsize)

if nargin<3
    figsize = [300 300];
end

if nargin<2
    layout = [2 ceil(length(fh_idv)/2)];
end

fh_all = figure('Position', [50 100 figsize(1) figsize(2)]);
for i = 1:length(fh_idv)
    old_ax = findobj(fh_idv{i}, 'Type', 'axes');
    new_ax = subplot(layout(1),layout(2),i);
    copyobj(old_ax.Children, new_ax);

    title(['session ', num2str(i)])
    xlabel('PC 1')
    ylabel('PC 2')
end

end