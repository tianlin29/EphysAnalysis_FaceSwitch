function fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, layout, figsize)

if nargin<3
    figsize = [300 300];
end

if nargin<2
    layout = [2 ceil(length(fh_idv)/2)];
end

opt.xlim = log([3 100]/100);
opt.xtick = log([3 6 12 24 48 96]'/100);
opt.xticklabel = {'0', '6', '', '24', '', '96'};

fh_all = figure('Position', [50 100 figsize(1) figsize(2)]);
nax = min(length(fh_idv), layout(1)*layout(2));
for i = 1:nax
    old_ax = get(fh_idv{i}, 'Children');
    new_ax = subplot(layout(1),layout(2),i);
    copyobj(old_ax.Children, new_ax);

    title(['session ', num2str(i)])
    xtickangle(0)
end
format_panel(gcf, ...
    'ylim', [.5 1], 'ytick', .5:.125:1, ...
    'xlim', opt.xlim, 'xtick', opt.xtick, 'xticklabel', opt.xticklabel, ...
    'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');

end