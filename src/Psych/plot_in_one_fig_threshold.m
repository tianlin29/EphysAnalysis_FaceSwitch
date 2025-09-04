function fh_all = plot_in_one_fig_threshold(fh_idv, layout, figsize)

if nargin<3
    figsize = [300 300];
end

if nargin<2
    layout = [2 ceil(length(fh_idv)/2)];
end

title_list = {'Nick', 'Woody'};
new_ax = nan(size(fh_idv));
fh_all = figure('Position', [50 100 figsize(1) figsize(2)]);
count = 0;
for i = 1:size(fh_idv,1)
    for j = 1:size(fh_idv,2)
        count = count + 1;
        
        if isempty(fh_idv{i,j})
            new_ax(i,j) = subplot(layout(1),layout(2),count);
            if i==1; title(title_list{j}); end
            continue;
        else
            old_ax = get(fh_idv{i,j}, 'Children');
            new_ax(i,j) = subplot(layout(1),layout(2),count);
            if i==1; title(title_list{j}); end
            copyobj(old_ax(end).Children, new_ax(i,j));
        end

        if i==1 && j==1
            xlabel('#Session')
            ylabel('Threshold')
        end
    end
end

end