function fh_all = plot_in_one_fig_unsigned_choice_switch_cost(fh_idv, layout, figsize)

if nargin<3
    figsize = [300 300];
end

if nargin<2
    layout = [2 ceil(length(fh_idv)/2)];
end

opt.xlim = log([3 100]/100);
opt.xtick = log([3 6 12 24 48 96]'/100);
opt.xticklabel = {'0', '6', '', '24', '', '96'};

title_list = {{'Nick learnTask2', 'first half'}, 'Nick learnTask3', 'Nick learnTask4', 'Woody learnTask3', 'Woody learnTask4', 'Second half'};
fh_all = figure('Position', [50 100 figsize(1) figsize(2)]);
count = 0;
for i = 1:2
    for j = 1:size(fh_idv,2)
        count = count+1;

        old_ax = get(fh_idv{i,j}, 'Children');
        new_ax = subplot(layout(1),layout(2),count);
        copyobj(old_ax.Children, new_ax);

        if count<=6; title(title_list{count}); end
        xtickangle(0)
    end
end
format_panel(gcf, ...
    'ylim', [.5 1], 'ytick', .5:.125:1, ...
    'xlim', opt.xlim, 'xtick', opt.xtick, 'xticklabel', opt.xticklabel, ...
    'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');

end