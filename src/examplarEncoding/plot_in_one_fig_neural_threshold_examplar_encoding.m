function fh_all = plot_in_one_fig_neural_threshold_examplar_encoding(fh_idv, layout, figsize, normalize_threshold)

if nargin<4
    normalize_threshold = false;
end

if nargin<3
    figsize = [300 300];
end

if nargin<2
    layout = [2 ceil(length(fh_idv)/2)];
end

title_list = {'Monkey N', 'Monkey W'};
new_ax = nan(size(fh_idv));
fh_all = figure('Position', [50 100 figsize(1) figsize(2)]);
count = 0;
for i = 1:size(fh_idv,1)
    for j = 1:size(fh_idv,2)
        count = count + 1;
        if isempty(fh_idv{i,j})
            if i==1; title(title_list{j}); end
            continue;
        else
            old_ax = get(fh_idv{i,j}, 'Children');
            new_ax(i,j) = subplot(layout(1),layout(2),count);
            title(title_list{i})
            copyobj(old_ax(1).Children, new_ax(i,j));
        end
        format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Neural thershold')
        xtickangle(0)
    end
end

end