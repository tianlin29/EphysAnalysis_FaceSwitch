run('../Initialize.m');
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask3'; % learnTask2, learnTask3, learnTask4, faceColor
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% check randomness in classification
% classification
n_files = 1;
data = cell(n_files, 3); % (session, condition)
for n = 1:n_files
    n
    % set options
    clear opt;
    opt.epoch = 1;
    opt.tstamp = {[200,300]}; % 对于不同位置的时刻点，使用的是不同的随机种子，因此结果会不同
    opt.t_win = 200;
    opt.glm_type = 'lasso';

    % load neural data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2); % extract epoch 2

    % extract task variables
    task_set = fnd.getp('task_set');
    targ_cor = fnd.getp('targ_cor');

    % project pair 1 to pair 1
    rng(1);
    ID = targ_cor; ID(task_set~=1) = NaN;
    opt.Kfold = 5;
    data{n,1} = PopulationClassifier(fnd, ID, opt);

    % project pair 2 to pair 2
    rng(1);
    ID = targ_cor; ID(task_set~=2) = NaN;
    opt.Kfold = 5;
    data{n,2} = PopulationClassifier(fnd, ID, opt);

    % build pair 1 classifier
    rng(1);
    ID = targ_cor; ID(task_set~=1) = NaN;
    opt.Kfold = 1;
    data_pair1_axis = PopulationClassifier(fnd, ID, opt);
        
    % project pair 2 to pair 1
    rng(1);
    ID = targ_cor; ID(task_set~=2) = NaN;
    opt.Kfold = 1;
    opt.refBeta = data_pair1_axis.data_ind; % attach pair 1 classification result
    data{n,3} = PopulationClassifier(fnd, ID, opt);
end
save(fullfile(InterimDir, sprintf('classifier_test.mat')), 'data', 'opt')

% load data
D = load(fullfile(InterimDir, sprintf('classifier_test.mat'))); % get data and opt
data = D.data;
opt_classifier = D.opt;

% merge data from all sessions and conditions
[coh, cor, cond, session] = deal([]);
for n = 1:n_files
    for c = 1:3
        ntrial = length(data{n,c}.param.morph_level);
        coh = [coh; data{n,c}.param.morph_level];
        cor = [cor; data{n,c}.Correct{1}];
        cond = [cond; c*ones(ntrial,1)];
        session = [session; n*ones(ntrial,1)];
    end
end

% plot psychometric function
n_time = length(opt_classifier.tstamp{1});
thres = nan(n_files, n_time, 3);
for t = 1:n_time
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
        case 'faceColor'
            opt.color = [0 0 0; 1 0 0; 1 0 0];
    end
    opt.linewidth = [0.5, 0.5, 1.5];
    opt.average = false; % do not average across learning sessions
    opt.normalize_threshold = false;

    [~, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(cond, coh, cor(:,t), session, opt);
end


%% classifier
data = cell(n_files, 3); % (session, condition)
for n = 1:n_files
    n
    % set options
    clear opt;
    opt.epoch = 1;
    opt.tstamp = {100:100:600}; % 对于不同位置的时刻点，使用的是不同的随机种子，因此结果会不同
    opt.t_win = 200;
    opt.glm_type = 'lasso';

    % load neural data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2); % extract epoch 2

    % extract task variables
    task_set = fnd.getp('task_set');
    targ_cor = fnd.getp('targ_cor');

    % project pair 1 to pair 1
    rng(1);
    ID = targ_cor; ID(task_set~=1) = NaN;
    data{n,1} = PopulationClassifier(fnd, ID, opt);

    % project pair 2 to pair 2
    rng(1);
    ID = targ_cor; ID(task_set~=2) = NaN;
    data{n,2} = PopulationClassifier(fnd, ID, opt);

    % project pair 2 to pair 1
    rng(1);
    ID = targ_cor; ID(task_set~=2) = NaN;
    opt.refBeta = data{n,1}.data_ind; % attach pair 1 classification result
    data{n,3} = PopulationClassifier(fnd, ID, opt);
end
save(fullfile(InterimDir, sprintf('classifier_%s_%s.mat', monkey, experiment)), 'data', 'opt')

%% plot
% load data
D = load(fullfile(InterimDir, sprintf('classifier_%s_%s.mat', monkey, experiment))); % get data and opt
data = D.data;
opt_classifier = D.opt;

% merge data from all sessions and conditions
[coh, cor, cond, session] = deal([]);
for n = 1:n_files
    for c = 1:3
        ntrial = length(data{n,c}.param.morph_level);
        coh = [coh; data{n,c}.param.morph_level];
        cor = [cor; data{n,c}.Correct{1}];
        cond = [cond; c*ones(ntrial,1)];
        session = [session; n*ones(ntrial,1)];
    end
end

% plot psychometric function
n_time = length(opt_classifier.tstamp{1});
thres = nan(n_files, n_time, 3);
for t = 1:n_time
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
        case 'faceColor'
            opt.color = [0 0 0; 1 0 0; 1 0 0];
    end
    opt.linewidth = [0.5, 0.5, 1.5];
    opt.average = false; % do not average across learning sessions
    opt.normalize_threshold = false;

    [~, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(cond, coh, cor(:,t), session, opt);
    thres_ = cellfun(@(x) x.thres, stat', 'uni', 0); thres_ = cell2mat(thres_); thres(:,t,:) = thres_; % save threshold
    fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [5 5], [500 500]*1.5);
    print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_%d_%d_ms.pdf', ...
        monkey, experiment, opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2)));
    print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s_%d_%d_ms.pdf', ...
        monkey, experiment, opt_classifier.tstamp{1}(t)-opt_classifier.t_win/2, opt_classifier.tstamp{1}(t)+opt_classifier.t_win/2)));
end

I = abs(thres)>nanmean(thres(:))+3*nanstd(thres(:)) | thres<0;
if opt.normalize_threshold
    thres = thres ./ thres(:,:,1);
end
thres(I) = NaN;
ylim_range = round([nanmin(thres(:))-0.1*nanmean(thres(:)), nanmax(thres(:))+0.1*nanmean(thres(:))], 1);

% plot how threshold changes with session
fh_session = figure('Position', [50 100 600 600]);
for n = 1:n_files
    subplot(5,5,n); hold on
    title(sprintf('Session %d', n))
    for c = 1:3
        plot(thres(n,:,c), '.-', 'markers', 7, 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
    end
    format_panel(gca, 'xlabel', '#Time bin', 'ylabel', 'Threshold', ...
        'xlim', [1-0.5 n_time+0.5], 'xtick', 1:2:n_time, 'ylim', ylim_range)
    xtickangle(0)
end
print(fh_session, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_session.pdf', monkey, experiment)))

% plot how threshold changes with time in a trial
fh_timebin = figure('Position', [50 100 600 600]);
for t = 1:n_time
    subplot(5,5,t); hold on
    title(sprintf('Time bin %d', t))
    for c = 1:3
        plot(thres(:,t,c), '.-', 'markers', 7, 'Color', opt.color(c,:), 'LineWidth', opt.linewidth(c));
    end
    format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Threshold', ...
        'xlim', [1-0.5 n_files+0.5], 'xtick', 1:2:n_files, 'ylim', ylim_range)
    xtickangle(0)
end
print(fh_timebin, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_timebin.pdf', monkey, experiment)))

