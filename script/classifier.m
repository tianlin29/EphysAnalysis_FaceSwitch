run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% classifier
[data_1to1, data_2to2, data_1to2] = deal(cell(1, n_files));
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % set options
    clear opt;
    opt.epoch = 1;
    opt.glm_type = 'lasso';
    opt.tstamp = {200};
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
    case 'faceColor'
        opt.color = [0 0 0; 1 0 0; 1 0 0];
end
opt.linewidth = [0.5, 0.5, 1.5];
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(task_set, morph_level, correct, session, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s.pdf', monkey, experiment)));

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



