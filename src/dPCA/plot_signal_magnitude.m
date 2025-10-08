function [fh_pair, fh_choice] = plot_signal_magnitude(monkey_list, experiment_list, InterimDir)

[task_signal, choice_signal] = deal(cell(2, 3)); % {monkey, experiment}
for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        [~, n_files] = get_file_path(monkey, experiment);
        for session_id = 1:n_files
            data = load(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d.mat', monkey, experiment, session_id))).data;
            dpc = data.dpc; % (dim, time, cond) 2task*2choice
            tstamp = data.tstamp; % (1, time)

            % smooth dPC score
            for d = 1:size(dpc, 1)
                for c = 1:size(dpc, 3)
                    dpc(d,:,c) = nanconv(dpc(d,:,c), fspecial('average', [1,100]), 'same');
                end
            end

            % calculate task signal difference in dPC score
            dpc_mn = [mean(dpc(1,:,[1 3]), 3)', mean(dpc(1,:,[2 4]), 3)'];
            I = tstamp>200 & tstamp<600;
            task_signal{monkey_id, exp_id}(session_id) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));

            % calculate choice signal difference in dPC score
            dpc_mn = [mean(dpc(2,:,[1 2]), 3)', mean(dpc(2,:,[3 4]), 3)'];
            I = tstamp>200 & tstamp<600;
            choice_signal{monkey_id, exp_id}(session_id) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));
        end
    end
end

color_list = [44 145 224;
    58 191 153;
    240 169 58]/255;

% plot pair signal version 2 (plot 2 monkeys together)
normalized_task_signal = cellfun(@(x,y) x./y, task_signal, choice_signal, 'uni', 0);
task_signal_to_plot = normalized_task_signal;

fh_pair = figure('Position', [50 100 300*1.2 200*1.2]);
for exp_id = 1:length(experiment_list)
    subplot(2,3,exp_id); hold on
    nses = max([length(task_signal_to_plot{1, exp_id}), length(task_signal_to_plot{2, exp_id})]);
    for monkey_id = 1:length(monkey_list)
        if exp_id==1 && monkey_id==2; continue; end
        
        if monkey_id==1
            plot(task_signal_to_plot{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7)
        else
            plot(task_signal_to_plot{monkey_id, exp_id}, '-', 'Color', color_list(exp_id,:), 'LineWidth', 1);
            scatter(1:length(task_signal_to_plot{monkey_id, exp_id}), task_signal_to_plot{monkey_id, exp_id}, 7, 'filled', 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', color_list(exp_id,:))
        end
    end
    format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Pair signal', ...
        'xlim', [0.5 nses+0.5], ...
        'xtick', 1:2:nses, ...
        'ylim', [0 1])
    xtickangle(0)
    title(sprintf('Pair 1 vs. %d', exp_id+1))
end
axPosition = gca().Position; % 0.7320    0.1307    0.1730    0.3204

% plot choice signal
fh_choice = figure('Position', [50 100 300*1.2 200*1.2]);
ax = subplot(2,3,6); hold on
for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        nses = length(choice_signal{monkey_id, exp_id});
        if monkey_id==1
            plot(choice_signal{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7);
        else
            plot(choice_signal{monkey_id, exp_id}, '-', 'Color', color_list(exp_id,:), 'LineWidth', 1);
            scatter(1:nses, choice_signal{monkey_id, exp_id}, 7, 'filled', 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', color_list(exp_id,:))
        end
    end
end
format_panel(fh_choice, 'xlabel', '#Session', 'ylabel', 'Choice signal', ...
    'xlim', [0.5 15+0.5], ...
    'xtick', 1:2:15, ...
    'ylim', [20 120])
xtickangle(0)
ax.Position = axPosition;

end


