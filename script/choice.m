run('../Initialize.m');
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask3'; % learnTask2, learnTask3, learnTask4, faceColor
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'choice'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'choice'); mkdir(InterimDir);

%%
D = combineTrialData(fullfile(PreprocDir, sprintf('%s_%s.mat', monkey, experiment)));

%% block design
session = cellfun(@(x) x.session_id, D);
cond = cellfun(@(x) x.cond, D);
fh = figure;
for n = 1:max(session)
    subplot(5,5,n); hold on
    cond_ = cond(session==n)';
    h = imagesc(cond_, [1 2]);
    switch experiment
        case {'learnTask2', 'faceColor'}
            map = [0 0 0; 44 145 224]/255;
        case 'learnTask3'
            map = [0 0 0; 58 191 153]/255;
        case 'learnTask4'
            map = [0 0 0; 240 169 58]/255;
    end
    colormap(map)
    set(h, 'AlphaData', 0.4*(ones(size(cond_))));
    title(sprintf('session %d', n))
end
format_panel(fh, 'axis', 'normal', ...
    'xtick', 0:200:1000, 'xlabel', '#Trial');
print(fh, '-dpdf', fullfile(FigDir, sprintf('block_design_%s_%s.pdf', monkey, experiment)));

%% block design for paper
figure;
subplot(2,1,1); hold on
opt.color = [0 0 0; 44 145 224; 58 191 153; 240 169 58]/255;
x = [0 5; 2 3; 3 4; 4 5];
y = [1 1; 2 2; 3 3; 4 4];
for i = 1:4
    plot(x(i,:), y(i,:), 'Color', opt.color(i,:))
end
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Stimulus set', 'axis', 'normal', 'ylim', [0 5])
set(gca, 'YDir', 'reverse', 'XTick', [], 'YTick', []);

subplot(2,1,2)
fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;
task_set = fnd.getp('task_set'); task_set = task_set(1,:);
imagesc(task_set, [1 2]); alpha(0.4);
map = [0 0 0; 44 145 224]/255;
colormap(map)
format_panel(gca, 'axis', 'normal', 'xlim', [0 length(task_set)], 'xtick', 0:400:1200, 'xlabel', '#Trial');
yticks([])
set(gcf, 'Position', [100 100 120 150])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('block_design.pdf')));

%% accuracy
monkey = 'Nick';
correct_rate = nan(4,2); % (pair, before/after)
for n = 1:3
    experiment = sprintf('learnTask%d', n+1);
    [~, n_files] = get_file_path(monkey, experiment);
    for i = 1:2
        if i==1; date_id = 1; else; date_id = n_files; end
        fnd = load(get_file_path(monkey, experiment, date_id, 'FND_sorted')).fnd;
        fnd = fnd.extract_trial(~isnan(fnd.getp('targ_cho')));
        targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
        targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:);
        task_set = fnd.getp('task_set'); task_set = task_set(1,:);
        correct = targ_cho==targ_cor;

        if n==1 && i==1
            correct_rate(n,2) = mean(correct(task_set==1));
        end
        correct_rate(n+1,i) = mean(correct(task_set==2));
    end
end

% plot
figure('Position', [100 100 220 120]);
b = bar(correct_rate, 1, 'FaceColor', 'flat', 'EdgeColor', 'none');
color_list = [55 159 48; 240 67 50]/255;
for k = 1:2
    b(k).CData = color_list(k,:);
end
format_panel(gca, 'xtick', 1:4, 'xticklabel', {'Pair 1', 'Pair 2', 'Pair 3', 'Pair 4'}, ...
    'ylabel', 'Performance', 'ylim', [0.4 1])
title('Data')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('accuracy_data_RNN_model.pdf')))

%% unsigned choice
session = cellfun(@(x) x.session_id, D);
cond = cellfun(@(x) x.cond, D);
coh = cellfun(@(x) x.coh, D);
targ_cor = cellfun(@(x) x.targ_cor, D);
resp = cellfun(@(x) x.resp, D);
cor = resp==targ_cor;

clear opt
% select data
opt.session_list = 1:n_files;
% process data
opt.log = true;
opt.constant = true;
opt.verbose = false;
% plot
switch experiment
    case {'learnTask2', 'faceColor'}
        opt.color = [0 0 0; 44 145 224]/255;
    case 'learnTask3'
        opt.color = [0 0 0; 58 191 153]/255;
    case 'learnTask4'
        opt.color = [0 0 0; 240 169 58]/255;
end
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond, coh, cor, session, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [5 5], [500 500]*1.4);
save(fullfile(InterimDir, sprintf('unsigned_choice_%s_%s.mat', monkey, experiment)), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('unsigned_choice_summary_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_choice_%s_%s.pdf', monkey, experiment)));

%% plot threshold of all sessions together
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

fh_summary = cell(2, 3); % {monkey, experiment}
for exp_id = 1:length(experiment_list)
    for monkey_id = 1:length(monkey_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end
        D = combineTrialData(fullfile(PreprocDir, sprintf('%s_%s.mat', monkey, experiment)));
        [~, n_files] = get_file_path(monkey, experiment);

        session = cellfun(@(x) x.session_id, D);
        cond = cellfun(@(x) x.cond, D);
        coh = cellfun(@(x) x.coh, D);
        targ_cor = cellfun(@(x) x.targ_cor, D);
        resp = cellfun(@(x) x.resp, D);
        cor = resp==targ_cor;

        clear opt
        % select data
        opt.session_list = 1:n_files;
        % process data
        opt.log = true;
        opt.constant = false;
        opt.verbose = false;
        % plot
        switch experiment
            case 'learnTask2'
                opt.color = [0 0 0; 44 145 224; 44 145 224]/255;
            case 'learnTask3'
                opt.color = [0 0 0; 58 191 153; 58 191 153]/255;
            case 'learnTask4'
                opt.color = [0 0 0; 240 169 58; 240 169 58]/255;
        end
        opt.average = false; % do not average across learning sessions

        [~, ~, fh_summary{monkey_id, exp_id}, ~] = run_unsigned_choice_2cond(cond, coh, cor, session, opt);
    end
end

fh_all = plot_in_one_fig_threshold(fh_summary, [2 3], [300 200]*1.2);
print(fh_all, '-dpdf', fullfile(FigDir, 'threshold.pdf'));

%% [stat] accuracy of unsigned choice
stat = load(fullfile(InterimDir, sprintf('unsigned_choice_%s_%s.mat', monkey, experiment))).stat;
nses = length(stat);
acc1 = cellfun(@(x) x.acc1, stat);
acc2 = cellfun(@(x) x.acc2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', acc1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Accuracy of pair 1 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Accuracy of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
else
    error('check')
end
[rho, p_val] = corr((1:nses)', acc2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Accuracy of pair 2 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Accuracy of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
else
    error('check')
end

%% [stat] threshold of unsigned choice
stat = load(fullfile(InterimDir, sprintf('unsigned_choice_%s_%s.mat', monkey, experiment))).stat;
nses = length(stat);
thres1 = cellfun(@(x) x.thres1, stat);
thres2 = cellfun(@(x) x.thres2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', thres1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Threshold of pair 1 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Threshold of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Threshold of pair 1 decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end
[rho, p_val] = corr((1:nses)', thres2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Threshold of pair 2 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Threshold of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Threshold of pair 2 decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

%% choice
session = cellfun(@(x) x.session_id, D);
cond = cellfun(@(x) x.cond, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);

clear opt
% select data
opt.session_list = 1:n_files;
% process data
opt.constant = false;
opt.verbose = false;
% plot
switch experiment
    case {'learnTask2', 'faceColor'}
        opt.color = [0 0 0; 44 145 224]/255;
    case 'learnTask3'
        opt.color = [0 0 0; 58 191 153]/255;
    case 'learnTask4'
        opt.color = [0 0 0; 240 169 58]/255;
end
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_choice_2cond(cond, coh, resp, session, opt);
fh_all = plot_in_one_fig_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
save(fullfile(InterimDir, sprintf('choice_%s_%s.mat', monkey, experiment)), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('choice_summary_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('choice_%s_%s.pdf', monkey, experiment)));

%% plot choice for the first session of each pair
file_path = {'\\10.10.49.250\rawdata\FaceSwitch_Monkey\Nick\20250528\Nick20250528_01.mat', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Nick\20250704\Nick20250704_01.mat', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Nick\20250721\Nick20250721_01.mat';
    '', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Woody\20231103\Woody20231103_01.mat', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Woody\20240110\Woody20240110_01.mat'};
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

for monkey_idx = 1:2
    for experiment_idx = 1:3
        monkey = monkey_list{monkey_idx};
        experiment = experiment_list{experiment_idx};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        D = load(file_path{monkey_idx, experiment_idx}).trial_data;
        D = D(~isnan(cellfun(@(x) x.response, D))); % responded trial only

        session = ones(size(D));
        cond = cellfun(@(x) x.curr_task_set, D);
        coh = cellfun(@(x) x.morph_level(1), D);
        resp = cellfun(@(x) 3-x.response, D);

        clear opt
        % select data
        opt.session_list = 1;
        % process data
        opt.constant = false;
        opt.verbose = false;
        % plot
        switch experiment_idx
            case 1
                opt.color = [0 0 0; 44 145 224]/255;
            case 2
                opt.color = [0 0 0; 58 191 153]/255;
            case 3
                opt.color = [0 0 0; 240 169 58]/255;
        end
        opt.average = false; % do not average across learning sessions
        [~, fh_idv, ~, stat] = run_choice_2cond(cond, coh, resp, session, opt);
        print(fh_idv, '-dpdf', fullfile(FigDir, sprintf('choice_of_the_first_session_%s_%s.pdf', monkey, experiment)));
    end
end

%% choice switch cost
session = cellfun(@(x) x.session_id, D);
n_after_switch = cellfun(@(x) x.n_after_switch, D);
cond_switch = nan(size(n_after_switch)); idx = 30;
cond_switch(n_after_switch>=idx) = 1; cond_switch(n_after_switch<idx) = 2; 
coh = cellfun(@(x) x.coh, D);
targ_cor = cellfun(@(x) x.targ_cor, D);
resp = cellfun(@(x) x.resp, D);
cor = resp==targ_cor;

clear opt
% select data
opt.session_list = 1:n_files;
% process data
opt.log = true;
opt.constant = false;
opt.verbose = false;
% plot
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond_switch, coh, cor, session, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [5 5], [500 500]*1.5);
save(fullfile(InterimDir, sprintf('choice_switch_cost_%s_%s.mat', monkey, experiment)), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('summary_choice_switch_cost_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('choice_switch_cost_%s_%s.pdf', monkey, experiment)));

%% [stat] switch cost in glmfit
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

p_summary = cell(2, 3); % {monkey, experiment}
fh = figure('Position', [50 100 600 300]);
for exp_id = 1:length(experiment_list)
    for monkey_id = 1:length(monkey_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        stat = load(fullfile(InterimDir, sprintf('choice_switch_cost_%s_%s.mat', monkey, experiment))).stat;
        p_summary{monkey_id, exp_id} = cellfun(@(x) x.stat_glmfit.p(2), stat);

        subplot(2, 3, (monkey_id-1)*3 + exp_id); hold on
        plot(p_summary{monkey_id, exp_id}, '.-', 'Color', 'black')
        plot(xlim, [0.05 0.05], '-', 'Color', 'red')
        title(sprintf('%s %s', monkey, experiment))
    end
end
format_panel(fh, 'xlabel', '#Session', 'ylabel', 'p-value', 'ytick', [0 0.05 0.5])
print(fh, '-dpdf', fullfile(FigDir, sprintf('p_value_choice_switch_cost.pdf')));

%% choice switch cost 比较前半部分和后半部分
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

fh = cell(2, 5); % {first+half, monkey*experiment}
count = 0;
for exp_id = 1:length(experiment_list)
    for monkey_id = 1:length(monkey_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end
        D = combineTrialData(fullfile(PreprocDir, sprintf('%s_%s.mat', monkey, experiment)));
        [~, n_files] = get_file_path(monkey, experiment);
        count = count + 1;

        session = cellfun(@(x) x.session_id, D);
        n_after_switch = cellfun(@(x) x.n_after_switch, D);
        cond_switch = nan(size(n_after_switch)); idx = 30;
        cond_switch(n_after_switch>=idx) = 1; cond_switch(n_after_switch<idx) = 2;
        coh = cellfun(@(x) x.coh, D);
        targ_cor = cellfun(@(x) x.targ_cor, D);
        resp = cellfun(@(x) x.resp, D);
        cor = resp==targ_cor;

        clear opt
        % process data
        opt.log = true;
        opt.constant = false;
        opt.verbose = false;
        % plot
        opt.average = true;
        % opt.legend = {'Non-switch', 'switch'};

        % select data
        opt.session_list = 1:round(n_files/2);
        [fh{1,count}, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond_switch, coh, cor, session, opt);

        opt.session_list = round(n_files/2):n_files;
        [fh{2,count}, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond_switch, coh, cor, session, opt);
    end
end

fh_all = plot_in_one_fig_unsigned_choice_switch_cost(fh, [2 5], [500 200]*1.5);
print(fh_all, '-dpdf', fullfile(FigDir, 'first_half_vs_second_half_sessions.pdf'));

%% 线性回归模型
targ_cor = cellfun(@(x) x.targ_cor, D);
resp = cellfun(@(x) x.resp, D);
cor = double(resp==targ_cor);
n_after_switch = cellfun(@(x) x.n_after_switch, D);
session = cellfun(@(x) x.session_id, D);

tbl = table(cor, n_after_switch, session, VariableNames={'correct', 'n_after_switch', 'session'});
modelspec = 'correct ~ n_after_switch + session';
mdl = fitglm(tbl, modelspec, 'Distribution', 'binomial', 'Link', 'logit');

%% [test] accuracy vs. switch
session = cellfun(@(x) x.session_id, D);
cond = cellfun(@(x) x.cond, D);
n_after_switch = cellfun(@(x) x.n_after_switch, D);
targ_cor = cellfun(@(x) x.targ_cor, D);
resp = cellfun(@(x) x.resp, D);
cor = resp==targ_cor;

figure;
for n = 1:max(session)
    u_n_after_switch = 0:79;
    I1 = cond==1 & session==n;
    I2 = cond==2 & session==n;
    [p1, pse1] = calcGroupMean(cor(I1), n_after_switch(I1), u_n_after_switch, 'binary');
    [p2, pse2] = calcGroupMean(cor(I2), n_after_switch(I2), u_n_after_switch, 'binary');

    clear opt
    opt.color = [0 0 0; 1 0 0];

    subplot(5,5,n); hold on
    plot(u_n_after_switch, p1, '.', 'markers', 7, 'Color', opt.color(1,:));
    cerrorbar(u_n_after_switch, p1, pse1, 'Color', opt.color(1,:));
    plot(u_n_after_switch, p2, '.', 'markers', 7, 'Color', opt.color(2,:));
    cerrorbar(u_n_after_switch, p2, pse2, 'Color', opt.color(2,:));
    title(sprintf('session %d', n))
end







