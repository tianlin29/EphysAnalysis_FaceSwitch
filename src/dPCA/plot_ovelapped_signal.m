function [fh, fh_avg, fh_first_second, fh_pair1_2] = plot_ovelapped_signal(monkey, experiment, n_files, InterimDir, list)

% average signal
[task_signal] = deal(nan(801, 2, n_files, 2)); % (time, 2 cond, session, first/second)
for t = 1:2
    for n = 1:n_files
        data = load(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%s.mat', monkey, experiment, n, list{t}))).data;
        dpc = data.dpc; % (dim, time, cond) 2task*2choice
        tstamp = data.tstamp; % (1, time)

        % smooth dPC score
        for d = 1:size(dpc, 1)
            for c = 1:size(dpc, 3)
                dpc(d,:,c) = nanconv(dpc(d,:,c), fspecial('average', [1,100]), 'same');
            end
        end

        % calculate task signal difference in dPC score
        dpc_mn = [mean(dpc(1,:,[1 3]), 3)', mean(dpc(1,:,[2 4]), 3)']; % task signal of task 1 and 2
        mn = mean(dpc_mn, 1);
        if mn(1)>mn(2) % task 1 is negative
            dpc_mn = -dpc_mn;
        end
        task_signal(:,:,n,t) = dpc_mn;
    end
end

% prepare plotting options
switch experiment
    case {'learnTask2', 'passiveLong'}
        opt.plot.color = [0 0 0; 44 145 224]/255;
    case 'learnTask3'
        opt.plot.color = [0 0 0; 58 191 153]/255;
    case 'learnTask4'
        opt.plot.color = [0 0 0; 240 169 58]/255;
    case 'faceColor'
        opt.plot.color = [0 0 0; 1 0 0];
end
opt.plot.linewidth = [0.5 2];
opt.plot.linestyle = {'-', '--'};

% plot overlapped signal of 2task*first/second
fh = figure('Position', [50 100 900 900]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for t = 1:2
        for c = 1:2
            plot(tstamp, task_signal(:,c,n,t), 'Color', opt.plot.color(c,:), 'LineWidth', opt.plot.linewidth(t))
        end
    end
    plot(xlim, [0 0], ':', 'Color', 'black')
    title(sprintf('session %d', n))
end
format_panel(gcf, 'xlabel', 'Time (ms)', 'ylabel', 'Pair signal', 'ylim', [-60 60])

% test whether change of signal is significant
difference = task_signal(:,2,:,:) - task_signal(:,1,:,:); % pair2 - pair1

p = nan(n_files, 1);
for n = 1:n_files
    p(n) = run_ttest_paired(difference(:,:,n,1)', '>', difference(:,:,n,2)');
end

% plot average of significant sessions
% if strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2')
%     I = logical([1 1 0 0 1 1 1 0 1 1 1 1 1 1 0]);
% elseif strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask3')
%     I = logical([0 1 1 0 0 0 0 0 0 0]);
% elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask3')
%     I = logical([1 0 1 0 1 1 0 1 0 1 1 0 0]);
% end
I = true(size(1,n_files));
task_signal_mn = mean(task_signal(:,:,I,:), 3);

% opt.less_timepoint
t_idx = round(linspace(1, length(tstamp), 50));
ttt = tstamp(t_idx);
x = task_signal_mn(t_idx,:,:,:);

fh_avg = figure('Position', [50 100 110 110]);
hold on
for t = 1:2
    for c = 1:2
        plot(ttt, x(:,c,1,t), 'Color', opt.plot.color(c,:), 'LineWidth', opt.plot.linewidth(t))
    end
end
plot(xlim, [0 0], ':', 'Color', 'black')
title(sprintf('%s %s', monkey, experiment))
format_panel(gcf, 'xlabel', 'Time (ms)', 'ylabel', 'Pair signal', 'ylim', [-60 60], 'xlim', [tstamp(1) tstamp(end)])
xtickangle(0)

% plot difference between first and second half trials
fh_first_second = figure('Position', [50 100 900 900]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for t = 1:2
        plot(tstamp, difference(:,:,n,t), 'Color', 'black', 'LineWidth', opt.plot.linewidth(t))
    end
    plot(xlim, [0 0], ':', 'Color', 'black')
    if p(n)<0.05
        title(sprintf('session %d (*)', n), 'Color', 'red')
    else
        title(sprintf('session %d (ns)', n), 'Color', 'black')
    end
end
format_panel(gcf, 'xlabel', 'Time (ms)', 'ylabel', 'Pair signal', 'ylim', [-10 80])

% plot difference of pair 1 and 2
fh_pair1_2 = figure('Position', [50 100 900 900]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot(tstamp, (-task_signal(:,1,n,2)) - (-task_signal(:,1,n,1)), 'Color', opt.plot.color(1,:))
    plot(tstamp, task_signal(:,2,n,2) - task_signal(:,2,n,1), 'Color', opt.plot.color(2,:))
    plot(xlim, [0 0], ':', 'Color', 'black')
    title(sprintf('session %d', n))
end
format_panel(gcf, 'xlabel', 'Time (ms)', 'ylabel', 'Difference of pair signal', 'ylim', [-40 20])

end