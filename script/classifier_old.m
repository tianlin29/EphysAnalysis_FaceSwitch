run('../Initialize.m');
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask3'; % learnTask2, learnTask3, learnTask4, faceColor
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);


%% classifier
[data_1to1, data_2to2, data_1to2] = deal(cell(1, n_files));
if strcmp(experiment, 'faceColor')
    file_list = [1:8, 18:22];
else
    file_list = 1:n_files;
end

for n = 1:length(file_list)
    n
    fnd = load(get_file_path(monkey, experiment, file_list(n), 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % set options
    clear opt;
    opt.epoch = 1;
    opt.glm_type = 'lasso';
    opt.tstamp = {100:100:600};
    opt.t_win = 200;

    % get ID
    morph_level = fnd.getp('morph_level')*100;
    targ_cor = fnd.getp('targ_cor');
    targ_cho = fnd.getp('targ_cho');
    task_set = fnd.getp('task_set');
    trc = trial_classifier('plus_cat', 1);
    ID = trc.category(targ_cor); % classify trials based on stimulus category

    % project pair 1 to pair 1
    rng(1)
    ID_ = ID;
    ID_(task_set~=1) = NaN;
    trial_classifier_result(ID_, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});
    data_1to1{n} = PopulationClassifier(fnd, ID_, opt);

    % project pair 2 to pair 2
    rng(1)
    ID_ = ID;
    ID_(task_set~=2) = NaN;
    trial_classifier_result(ID_, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});
    data_2to2{n} = PopulationClassifier(fnd, ID_, opt);

    % project pair 1 to pair 2
    rng(1)
    ID_ = ID;
    ID_(task_set~=2) = NaN;
    trial_classifier_result(ID_, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});
    opt.refBeta = data_1to1{n}.data_ind; % attach pair 1 classification result
    data_1to2{n} = PopulationClassifier(fnd, ID_, opt);
end
save(fullfile(InterimDir, sprintf('classifier_%s_%s.mat', monkey, experiment)), 'data_1to1', 'data_2to2', 'data_1to2')

%% plot classifier
n_files = length(file_list);

% transform classification result to roma data format
load(fullfile(InterimDir, sprintf('classifier_%s_%s.mat', monkey, experiment)));
correct = [];
morph_level = [];
task_set = [];
session = [];
for n = 1:n_files
    % pair 1
    correct_ = data_1to1{n}.Correct{1};
    morph_level_ = data_1to1{n}.param.morph_level;

    correct = [correct; correct_];
    morph_level = [morph_level; morph_level_];
    task_set = [task_set; ones(length(morph_level_),1)];
    session = [session; n*ones(length(morph_level_),1)];

    % new pair
    correct_ = data_2to2{n}.Correct{1};
    morph_level_ = data_2to2{n}.param.morph_level;

    correct = [correct; correct_];
    morph_level = [morph_level; morph_level_];
    task_set = [task_set; 2*ones(length(morph_level_),1)];
    session = [session; n*ones(length(morph_level_),1)];

    % new pair project to pair 1 axis
    correct_ = data_1to2{n}.Correct{1};
    morph_level_ = data_1to2{n}.param.morph_level;

    correct = [correct; correct_];
    morph_level = [morph_level; morph_level_];
    task_set = [task_set; 3*ones(length(morph_level_),1)];
    session = [session; n*ones(length(morph_level_),1)];
end

% plot neural unsigned psychometric function
% loop through timepoints
tstamp = {100:100:600};
t_win = 200;
thres = nan(3, length(file_list), length(tstamp{1}));
for t = 1:length(tstamp{1})
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

    [~, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(task_set, morph_level, correct(:,t), session, opt);
    thres_ = cellfun(@(x) x.thres(:), stat, 'uni', 0); thres_ = cell2mat(thres_); thres(:,:,t) = thres_; % all the thresholds
    fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
    print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_%d_%d_ms.pdf', monkey, experiment, tstamp{1}(t)-t_win/2, tstamp{1}(t)+t_win/2)));
    print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s_%d_%d_ms.pdf', monkey, experiment, tstamp{1}(t)-t_win/2, tstamp{1}(t)+t_win/2)));
end

%% [test] 解码正确率随trial内时间和session的变化
thres(thres<0) = NaN;

% 3d mesh
figure; hold on
[X, Y] = meshgrid(1:n_files, 1:length(tstamp{1}));
mesh(X, Y, squeeze(thres(1,:,:))', 'FaceColor', opt.color(1,:), 'FaceAlpha', .3);
mesh(X, Y, squeeze(thres(2,:,:))', 'FaceColor', opt.color(2,:), 'FaceAlpha', .3);
mesh(X, Y, squeeze(thres(3,:,:))', 'FaceColor', opt.color(3,:), 'FaceAlpha', .9);
xlabel('#Session')
ylabel('Time bin')
zlabel('Threshold')

% 2d plot
figure('Position', [20 700 800 150]); % time
for t = 1:length(tstamp{1})
    subplot(1,length(tstamp{1}),t); hold on
    plot(1:n_files, squeeze(thres(1,:,t)), '.-', 'markers', 7, 'Color', opt.color(1,:), 'LineWidth', opt.linewidth(1));
    plot(1:n_files, squeeze(thres(2,:,t)), '.-', 'markers', 7, 'Color', opt.color(2,:), 'LineWidth', opt.linewidth(2));
    plot(1:n_files, squeeze(thres(3,:,t)), '.-', 'markers', 7, 'Color', opt.color(3,:), 'LineWidth', opt.linewidth(3));
    title(sprintf('Time bin %d', t))
end
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Threshold', 'ylim', [0.25 2])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_timebin.pdf', monkey, experiment)));

figure('Position', [20 170 800 150*3]); % session
for n = 1:n_files
    subplot(3,6,n); hold on
    plot(1:length(tstamp{1}), squeeze(thres(1,n,:)), '.-', 'markers', 7, 'Color', opt.color(1,:), 'LineWidth', opt.linewidth(1));
    plot(1:length(tstamp{1}), squeeze(thres(2,n,:)), '.-', 'markers', 7, 'Color', opt.color(2,:), 'LineWidth', opt.linewidth(2));
    plot(1:length(tstamp{1}), squeeze(thres(3,n,:)), '.-', 'markers', 7, 'Color', opt.color(3,:), 'LineWidth', opt.linewidth(3));
    title(sprintf('Session %d', n))
end
format_panel(gcf, 'xlabel', 'Time bin', 'ylabel', 'Threshold', 'ylim', [0.25 2])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_session.pdf', monkey, experiment)));

%% [poster] decoding accuracy of an example session
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [2 3], [300 200]*1.2);
% set(gca, 'Position', [0.1300    0.6546    0.2134    0.2681])
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_example_session_%s_%s.pdf', monkey, experiment)));

%% [poster] plot everything together
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

        % pair-dependent decoder
        % transform classification result to roma data format
        load(fullfile(InterimDir, sprintf('classifier_%s_%s.mat', monkey, experiment)));
        correct = [];
        morph_level = [];
        task_set = [];
        session = [];
        [~, n_files] = get_file_path(monkey, experiment);
        for n = 1:n_files
            % pair 1
            correct_ = data_1to1{n}.Correct{1};
            morph_level_ = data_1to1{n}.param.morph_level;
            I_valid = ~isnan(correct_);

            correct = [correct; correct_(I_valid)];
            morph_level = [morph_level; morph_level_(I_valid)];
            task_set = [task_set; ones(sum(I_valid),1)];
            session = [session; n*ones(sum(I_valid),1)];

            % new pair
            correct_ = data_2to2{n}.Correct{1};
            morph_level_ = data_2to2{n}.param.morph_level;
            I_valid = ~isnan(correct_);

            correct = [correct; correct_(I_valid)];
            morph_level = [morph_level; morph_level_(I_valid)];
            task_set = [task_set; 2*ones(sum(I_valid),1)];
            session = [session; n*ones(sum(I_valid),1)];

            % new pair project to pair 1 axis
            correct_ = data_1to2{n}.Correct{1};
            morph_level_ = data_1to2{n}.param.morph_level;
            I_valid = ~isnan(correct_);

            correct = [correct; correct_(I_valid)];
            morph_level = [morph_level; morph_level_(I_valid)];
            task_set = [task_set; 3*ones(sum(I_valid),1)];
            session = [session; n*ones(sum(I_valid),1)];
        end

        % plot neural unsigned psychometric function
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
        opt.linewidth = [0.5, 0.5, 1.5];
        opt.average = false; % do not average across learning sessions
        % quick summary plot
        opt.normalize_threshold = true;

        [~, ~, fh_summary{monkey_id, exp_id}, ~] = run_unsigned_choice_3cond(task_set, morph_level, correct, session, opt);
    end
end

fh_all = plot_in_one_fig_neural_threshold(fh_summary, [2 3], [300 200]*1.2, opt.normalize_threshold);
print(fh_all, '-dpdf', fullfile(FigDir, 'neural_threshold_linear_classifier.pdf'));

%% compare decoding result and monkey's choice
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

[fh_correct, fh_wrong] = deal(cell(2, 3)); % {monkey, experiment}
for exp_id = 1:length(experiment_list)
    for monkey_id = 1:length(monkey_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        % pair-dependent decoder
        % transform classification result to roma data format
        load(fullfile(InterimDir, sprintf('classifier_%s_%s.mat', monkey, experiment)));
        correct = [];
        morph_level = [];
        task_set = [];
        cor_choice = [];
        session = [];
        [~, n_files] = get_file_path(monkey, experiment);
        for n = 1:n_files
            % pair 1
            correct_ = data_1to1{n}.Correct{1};
            morph_level_ = data_1to1{n}.param.morph_level;
            cor_choice_ = data_1to1{n}.param.targ_cor==data_1to1{n}.param.targ_cho;
            I_valid = ~isnan(correct_);

            correct = [correct; correct_(I_valid)];
            morph_level = [morph_level; morph_level_(I_valid)];
            task_set = [task_set; ones(sum(I_valid),1)];
            cor_choice = [cor_choice; cor_choice_(I_valid)];
            session = [session; n*ones(sum(I_valid),1)];

            % new pair
            correct_ = data_2to2{n}.Correct{1};
            morph_level_ = data_2to2{n}.param.morph_level;
            cor_choice_ = data_2to2{n}.param.targ_cor==data_2to2{n}.param.targ_cho;
            I_valid = ~isnan(correct_);

            correct = [correct; correct_(I_valid)];
            morph_level = [morph_level; morph_level_(I_valid)];
            task_set = [task_set; 2*ones(sum(I_valid),1)];
            cor_choice = [cor_choice; cor_choice_(I_valid)];
            session = [session; n*ones(sum(I_valid),1)];

            % new pair project to pair 1 axis
            correct_ = data_1to2{n}.Correct{1};
            morph_level_ = data_1to2{n}.param.morph_level;
            cor_choice_ = data_1to2{n}.param.targ_cor==data_1to2{n}.param.targ_cho;
            I_valid = ~isnan(correct_);

            correct = [correct; correct_(I_valid)];
            morph_level = [morph_level; morph_level_(I_valid)];
            task_set = [task_set; 3*ones(sum(I_valid),1)];
            cor_choice = [cor_choice; cor_choice_(I_valid)];
            session = [session; n*ones(sum(I_valid),1)];
        end

        [correct_trials, wrong_trials] = deal(nan(3, max(session)));
        for s = 1:max(session) % session
            for t = 1:3 % three decoding methods
                correct_tmp = correct(session==s & task_set==t);
                cor_choice_tmp = cor_choice(session==s & task_set==t);
                correct_trials(t,s) = sum(correct_tmp & cor_choice_tmp) / sum(cor_choice_tmp);
                wrong_trials(t,s) = sum(~correct_tmp & ~cor_choice_tmp) / sum(~cor_choice_tmp);
            end
        end

        clear opt;
        switch experiment
            case 'learnTask2'
                opt.color = [0 0 0; 44 145 224; 44 145 224]/255;
            case 'learnTask3'
                opt.color = [0 0 0; 58 191 153; 58 191 153]/255;
            case 'learnTask4'
                opt.color = [0 0 0; 240 169 58; 240 169 58]/255;
        end
        opt.linewidth = [0.5, 0.5, 1.5];

        fh_correct{monkey_id, exp_id} = figure; hold on
        for t = 1:3
            plot(correct_trials(t,:), 'Color', opt.color(t,:), 'LineWidth', opt.linewidth(t));
        end
        format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Correct prediction (%)')
        title('Correct trials')

        fh_wrong{monkey_id, exp_id} = figure; hold on
        for t = 1:3
            plot(wrong_trials(t,:), 'Color', opt.color(t,:), 'LineWidth', opt.linewidth(t));
        end
        format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Correct prediction (%)')
        title('Wrong trials')
    end
end

fh_correct_all = plot_in_one_fig_decoding_correlation(fh_correct, [2 3], [300 200]*1.2, [.6 1]);
print(fh_correct_all, '-dpdf', fullfile(FigDir, 'correct_trials_decoding_result.pdf'));

fh_wrong_all = plot_in_one_fig_decoding_correlation(fh_wrong, [2 3], [300 200]*1.2, [.3 .8]);
print(fh_wrong_all, '-dpdf', fullfile(FigDir, 'wrong_trials_decoding_result.pdf'));



