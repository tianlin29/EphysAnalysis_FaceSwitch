run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, threeExemplar, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier', experiment); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% classifier (pair axis for threeExemplar)
file_list = 1:n_files;
n_files = length(file_list);

REPEAT = 10;
for r = 1:REPEAT
    data = cell(n_files, 2); % (session, condition)
    for n = 1:n_files
        fprintf('repreat %d, session %d\n', r, n)

        % set options
        clear opt;
        opt.epoch = 1;
        opt.tstamp = {[100,200,300,400,500,600]};
        opt.t_win = 200;
        opt.glm_type = 'lasso';

        % load neural data
        fnd = load(get_file_path(monkey, experiment, file_list(n), 'FND_sorted')).fnd;
        if subsession==1; fnd = fnd.extract_trial(fnd.getp('targ_cor')==fnd.getp('targ_cho')); end % correct trials only for the main task
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        task_set = fnd.getp('task_set');
        targ_cor = fnd.getp('targ_cor');

        % for choice 1, get pair axis
        ID = task_set; ID(targ_cor~=1) = NaN;
        opt.Kfold = 5;
        data{n,1} = PopulationClassifier(fnd, ID, opt);

        % for choice 2, get pair axis
        ID = task_set; ID(targ_cor~=2) = NaN;
        opt.Kfold = 5;
        data{n,2} = PopulationClassifier(fnd, ID, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r)), 'data', 'opt')
end

%% classifier (category axis for main task)
if strcmp(experiment, 'faceColor')
    file_list = [1:8, 18:23];
else
    file_list = 1:n_files;
end
n_files = length(file_list);

REPEAT = 10;
for r = 1:REPEAT
    data = cell(n_files, 3); % (session, condition)
    for n = 1:n_files
        fprintf('repreat %d, session %d\n', r, n)

        % set options
        clear opt;
        opt.epoch = 1;
        opt.tstamp = {[100,200,300,400,500,600]}; % 对于不同位置的时刻点，使用的是不同的随机种子，因此结果会不同
        opt.t_win = 200;
        opt.glm_type = 'lasso';

        % load neural data
        fnd = load(get_file_path(monkey, experiment, file_list(n), 'FND_sorted')).fnd;
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        task_set = fnd.getp('task_set');
        targ_cor = fnd.getp('targ_cor');

        % project pair 1 to pair 1
        ID = targ_cor; ID(task_set~=1) = NaN;
        opt.Kfold = 5;
        data{n,1} = PopulationClassifier(fnd, ID, opt);

        % project pair 2 to pair 2
        ID = targ_cor; ID(task_set~=2) = NaN;
        opt.Kfold = 5;
        data{n,2} = PopulationClassifier(fnd, ID, opt);

        % build pair 1 classifier
        ID = targ_cor; ID(task_set~=1) = NaN;
        opt.Kfold = 1;
        data_category_axis = PopulationClassifier(fnd, ID, opt);

        % project pair 2 to pair 1
        ID = targ_cor; ID(task_set~=2) = NaN;
        opt.Kfold = 1;
        opt.refBeta = data_category_axis.data_ind; % attach pair 1 classification result
        data{n,3} = PopulationClassifier(fnd, ID, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r)), 'data', 'opt')
end


%% classifier (category axis for passiveLong)
REPEAT = 10;
for r = 1:REPEAT
    data = cell(n_files, 3); % (session, condition)
    for n = 1:n_files
        fprintf('repreat %d, session %d\n', r, n)

        % set options
        clear opt;
        opt.epoch = 1;
        opt.tstamp = {[100,200,300,400,500,600]}; % 对于不同位置的时刻点，使用的是不同的随机种子，因此结果会不同
        opt.t_win = 200;
        opt.glm_type = 'lasso';

        % (1)
        % load main task
        fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 1)).fnd;
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        targ_cor = fnd.getp('targ_cor');

        % category axis of main task
        ID = targ_cor;
        opt.Kfold = 5;
        data{n,1} = PopulationClassifier(fnd, ID, opt);

        % (2)
        % load passive fixation
        fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 2)).fnd;
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        targ_cor = fnd.getp('targ_cor');

        % category axis of passive fixation
        ID = targ_cor;
        opt.Kfold = 5;
        data{n,2} = PopulationClassifier(fnd, ID, opt);

        % build category axis
        ID = targ_cor;
        opt.Kfold = 1;
        data_category_axis = PopulationClassifier(fnd, ID, opt);

        % (3)
        % load main task
        fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 1)).fnd;
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        targ_cor = fnd.getp('targ_cor');

        % cross-experiment projection
        ID = targ_cor;
        opt.Kfold = 1;
        opt.refBeta = data_category_axis.data_ind; % attach main task classification result
        data{n,3} = PopulationClassifier(fnd, ID, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r)), 'data', 'opt')
end

%% plot threshold vs. session/timebin
% load data
n_repeats = 10;
n_conds = 3;
[coh, cor, cond, session, repeat, opt_classifier] = merge_across_repetition(InterimDir, monkey, experiment, n_repeats, n_files, n_conds);

% plot
[fh_session, fh_timebin] = plot_threshold(coh, cor, cond, session, repeat, opt_classifier, experiment);
print(fh_session, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_session.pdf', monkey, experiment)))
print(fh_timebin, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_timebin.pdf', monkey, experiment)))

%% plot accuracy vs. session/timebin
% load data
n_repeats = 10;
n_conds = 3;
[coh, cor, cond, session, repeat, opt_classifier] = merge_across_repetition(InterimDir, monkey, experiment, n_repeats, n_files, n_conds);

% plot
[fh_session, fh_timebin] = plot_accuracy(coh, cor, cond, session, repeat, opt_classifier, experiment);
print(fh_session, '-dpdf', fullfile(FigDir, sprintf('neural_accuracy_%s_%s_session.pdf', monkey, experiment)))
print(fh_timebin, '-dpdf', fullfile(FigDir, sprintf('neural_accuracy_%s_%s_timebin.pdf', monkey, experiment)))

%% plot psychometric function
% load data
n_repeats = 10;
n_conds = 2;
[coh, cor, cond, session, repeat, opt_classifier] = merge_across_repetition(InterimDir, monkey, experiment, n_repeats, n_files, n_conds);

% plot
[fh_summary, fh_all_psych] = plot_psychometric_function(coh, cor, cond, session, repeat, opt_classifier, experiment);
for t = 1:length(fh_summary)
    print(fh_summary{t}, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_%d_%d_ms.pdf', ...
        monkey, experiment, opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2)));
    print(fh_all_psych{t}, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s_%d_%d_ms.pdf', ...
        monkey, experiment, opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2)));
end

%% poster format
% pyschometric function
set(gcf, 'Position', [50 100 700 700])
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s_%d_%d_ms.pdf', ...
    monkey, experiment, opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2)));

% threshold
set(gcf, 'Position', [300 300 268 103])
format_panel(gca, 'ylim', [0 7], 'xlim', [0.5 15.5], 'xtick', 1:2:15);
xtickangle(0)

print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_%d_%d_ms.pdf', ...
    monkey, experiment, opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2)));

%% check nunit
nunit = nan(n_files, 1);
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    r = fnd.FR({2, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    [nunit(n), ntime, ntrial] = size(fnd.data{1});
end



