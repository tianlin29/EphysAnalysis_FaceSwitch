run('../Initialize.m');
monkey = 'Woody';
experiment = 'learnTask3';
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
    opt.tstamp = {400};
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
        opt.color = [99 97 172; 254 160 64; 254 160 64]/255;
    case 'learnTask3'
        opt.color = [99 97 172; 242 128 128; 242 128 128]/255;
    case 'learnTask4'
        opt.color = [99 97 172; 178 34 34; 178 34 34]/255;
end
opt.linewidth = [0.5, 0.5, 1.5];
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(task_set, morph_level, correct, session, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s.pdf', monkey, experiment)));

%% plot everything together
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

fh_summary = cell(3, 2); % {experiment, monkey}
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
                opt.color = [99 97 172; 254 160 64; 254 160 64]/255;
            case 'learnTask3'
                opt.color = [99 97 172; 242 128 128; 242 128 128]/255;
            case 'learnTask4'
                opt.color = [99 97 172; 178 34 34; 178 34 34]/255;
        end
        opt.linewidth = [0.5, 0.5, 1.5];
        opt.average = false; % do not average across learning sessions

        [~, ~, fh_summary{exp_id, monkey_id}, ~] = run_unsigned_choice_3cond(task_set, morph_level, correct, session, opt);
    end
end

fh_all = plot_in_one_fig_neural_threshold(fh_summary, [3 2], [200 300]*1.5);
print(fh_all, '-dpdf', fullfile(FigDir, 'neural_threshold.pdf'));




