run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask4'; % learnTask2, learnTask3, learnTask4, faceColor
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'exemplarEncoding'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'exemplarEncoding'); mkdir(InterimDir);

%% examplar encoding
[session, cond, coh, cor, category, choice] = deal([]);
for n = 1:n_files
    n
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % examplar decoding
    clear opt;
    opt.examplar_high_coh = false;

    [nunit, ntime, ntrial] = size(fnd.data{1});
    task_set = fnd.getp('task_set'); task_set = task_set(1,:);
    targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
    targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:);
    morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:);
    trial = 1:ntrial;
    FR = fnd.FR({1, [100 300]}); % (unit, trial)

    [distance_to_pair_1, distance_to_pair_2] = deal(nan(2, ntrial));
    for tr = 1:ntrial
        I_pair1_human  = targ_cor==targ_cho & targ_cor==1 & task_set==1;
        I_pair1_monkey = targ_cor==targ_cho & targ_cor==2 & task_set==1;
        I_pair2_human  = targ_cor==targ_cho & targ_cor==1 & task_set==2;
        I_pair2_monkey = targ_cor==targ_cho & targ_cor==2 & task_set==2;

        if opt.examplar_high_coh
            I_pair1_human = I_pair1_human & abs(morph_level)>=0.6;
            I_pair1_monkey = I_pair1_monkey & abs(morph_level)>=0.6;
            I_pair2_human = I_pair2_human & abs(morph_level)>=0.6;
            I_pair2_monkey = I_pair2_monkey & abs(morph_level)>=0.6;
        end

        examplar_pair1_human  = mean(FR(:, I_pair1_human & trial~=tr), 2); % (unit, 1)
        examplar_pair1_monkey = mean(FR(:, I_pair1_monkey & trial~=tr), 2);
        examplar_pair2_human  = mean(FR(:, I_pair2_human & trial~=tr), 2);
        examplar_pair2_monkey = mean(FR(:, I_pair2_monkey & trial~=tr), 2);

        current_trial = FR(:, tr);
        distance_to_pair_1(:, tr) = [pdist([examplar_pair1_human'; current_trial'], 'euclidean');
            pdist([examplar_pair1_monkey'; current_trial'], 'euclidean')];
        distance_to_pair_2(:, tr) = [pdist([examplar_pair2_human'; current_trial'], 'euclidean');
            pdist([examplar_pair2_monkey'; current_trial'], 'euclidean')];
    end

    % 1to1
    category_decoding_1to1 = nan(1,ntrial);
    category_decoding_1to1(distance_to_pair_1(1,:)<distance_to_pair_1(2,:)) = 1;
    category_decoding_1to1(distance_to_pair_1(2,:)<distance_to_pair_1(1,:)) = 2;
    correct_ = category_decoding_1to1==targ_cor;

    I_valid = task_set==1;
    cor = [cor; correct_(I_valid)'];
    coh = [coh; morph_level(I_valid)'];
    cond = [cond; ones(1,sum(I_valid))'];
    session = [session; n*ones(1,sum(I_valid))'];
    category = [category; targ_cor(I_valid)'];
    choice = [choice; targ_cho(I_valid)'];

    % 2to2
    category_decoding_2to2 = nan(1,ntrial);
    category_decoding_2to2(distance_to_pair_2(1,:)<distance_to_pair_2(2,:)) = 1;
    category_decoding_2to2(distance_to_pair_2(2,:)<distance_to_pair_2(1,:)) = 2;
    correct_ = category_decoding_2to2==targ_cor;

    I_valid = task_set==2;
    cor = [cor; correct_(I_valid)'];
    coh = [coh; morph_level(I_valid)'];
    cond = [cond; 2*ones(1,sum(I_valid))'];
    session = [session; n*ones(1,sum(I_valid))'];
    category = [category; targ_cor(I_valid)'];
    choice = [choice; targ_cho(I_valid)'];

    % 2to1
    category_decoding_2to1 = nan(1,ntrial);
    category_decoding_2to1(distance_to_pair_1(1,:)<distance_to_pair_1(2,:)) = 1;
    category_decoding_2to1(distance_to_pair_1(2,:)<distance_to_pair_1(1,:)) = 2;
    correct_ = category_decoding_2to1==targ_cor;

    I_valid = task_set==2;
    cor = [cor; correct_(I_valid)'];
    coh = [coh; morph_level(I_valid)'];
    cond = [cond; 3*ones(1,sum(I_valid))'];
    session = [session; n*ones(1,sum(I_valid))'];
    category = [category; targ_cor(I_valid)'];
    choice = [choice; targ_cho(I_valid)'];
end

if opt.examplar_high_coh
    save(fullfile(InterimDir, sprintf('examplar_encoding_high_coh_%s_%s.mat', monkey, experiment)), 'session', 'cond', 'coh', 'cor', 'category', 'choice');
else
    save(fullfile(InterimDir, sprintf('examplar_encoding_%s_%s.mat', monkey, experiment)), 'session', 'cond', 'coh', 'cor', 'category', 'choice');
end
beep

%% plot decoding result
load(fullfile(InterimDir, sprintf('examplar_encoding_%s_%s.mat', monkey, experiment)));

clear opt;
% select data
opt.session_list = [1:8, 18:20, 23];
% process data
opt.log = true;
opt.constant = false;
opt.verbose = false;
% plot
switch experiment
    case {'learnTask2', 'faceColor'}
        opt.color = [0 0 0; 44 145 224; 44 145 224]/255;
    case 'learnTask3'
        opt.color = [0 0 0; 58 191 153; 58 191 153]/255;
    case 'learnTask4'
        opt.color = [0 0 0; 240 169 58; 240 169 58]/255;
end
opt.linewidth = [0.5, 0.5, 1.5];
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_3cond(cond, coh, cor, session, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [5 5], [500 500]*1.5, [1:8, 18:20, 23]);
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_neural_psychometric_func_%s_%s.pdf', monkey, experiment)));

%% plot everything together
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

        [~, n_files] = get_file_path(monkey, experiment);
        load(fullfile(InterimDir, sprintf('examplar_encoding_%s_%s.mat', monkey, experiment)));

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
        opt.normalize_threshold = false;

        [~, ~, fh_summary{monkey_id, exp_id}, ~] = run_unsigned_choice_3cond(cond, coh, cor, session, opt);
    end
end

fh_all = plot_in_one_fig_neural_threshold_examplar_encoding(fh_summary, [2 3], [300 200]*1.2, opt.normalize_threshold);
print(fh_all, '-dpdf', fullfile(FigDir, 'neural_threshold_exemplar_encoding_high_coh.pdf'));

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

        load(fullfile(InterimDir, sprintf('examplar_encoding_%s_%s.mat', monkey, experiment)));
        cor_choice = choice==category;

        [correct_trials, wrong_trials] = deal(nan(3, max(session)));
        for s = 1:max(session) % session
            for t = 1:3 % three decoding methods
                correct_tmp = cor(session==s & cond==t);
                cor_choice_tmp = cor_choice(session==s & cond==t);
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









