run('../Initialize.m');
monkey = 'Nick';
experiment = 'faceColor';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

%% PSTH (face task)
DETREND = false;
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND')).fnd;

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100; morph_level(task_set==2) = NaN;
    color_level = fnd.getp('color_level')*100; color_level(task_set==1) = NaN;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    if DETREND
        trc = trial_classifier('stim_group', {[0 12], [12 24], [24 48], [48 Inf]}, 'plus_cat', 1, 'include_0coh', true);
    else
        trc = trial_classifier('plus_cat', 1, 'include_0coh', true);
    end
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    ID(task_set~=1) = NaN; % one task only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot PSTH
    psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], DETREND);

    clear opt;
    opt.cutoff = fnd.cutoff();
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.event_label = {fnd.alignto.event};
    opt.legend = arrayfun(@(x) num2str(x, '%.0f'), mean_coh, 'uni', 0);
    opt.legend_pos = [0.88 0.5 0.07 0.3];
    fh = showPopPSTH(fnd.tstamp, psth_data, opt);
    if DETREND
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_coh_face_task_detrend_session_%d.pdf', n)));
    else
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_coh_face_task_session_%d.pdf', n)));
    end
end

%% PSTH (color task)
DETREND = true;
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND')).fnd;

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100; morph_level(task_set==2) = NaN;
    color_level = fnd.getp('color_level')*100; color_level(task_set==1) = NaN;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('plus_cat', 1, 'include_0coh', true);
    [ID, mean_coh] = trc.stim_choice(color_level, targ_cho, targ_cor); % correct trials only
    ID(task_set~=2) = NaN; % one task only
    trial_classifier_result(ID, {'color_level', 'targ_cor', 'targ_cho', 'task_set'}, {color_level, targ_cor, targ_cho, task_set});

    % plot PSTH
    psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], DETREND);

    clear opt;
    opt.cutoff = fnd.cutoff();
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.event_label = {fnd.alignto.event};
    opt.legend = arrayfun(@(x) num2str(x, '%.0f'), mean_coh, 'uni', 0);
    opt.legend_pos = [0.88 0.5 0.07 0.3];
    fh = showPopPSTH(fnd.tstamp, psth_data, opt);
    if DETREND
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_coh_color_task_detrend_session_%d.pdf', n)));
    else
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_coh_color_task_session_%d.pdf', n)));
    end
end

%% PCA (face task)
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND')).fnd;

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100; morph_level(task_set==2) = NaN;
    color_level = fnd.getp('color_level')*100; color_level(task_set==1) = NaN;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('plus_cat', 1, 'include_0coh', true);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    ID(task_set~=1) = NaN; % one task only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot PCA
    clear opt
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.epoch = 2;
    opt.PC_range = [250 600];
    opt.PSTH_conv = {'boxcar', 100};
    opt.PC_kernel = {'bartlett', 200};

    fh = showPCA_trajectory(fnd, ID, opt);
    saveas(fh, fullfile(FigDir, sprintf('PCA_trajectory_face_task_session_%d.fig', n)));
    print(fh, '-dpdf', fullfile(FigDir, sprintf('PCA_trajectory_face_task_session_%d.pdf', n)));
end

%% PCA (color task)
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND')).fnd;

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100; morph_level(task_set==2) = NaN;
    color_level = fnd.getp('color_level')*100; color_level(task_set==1) = NaN;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('plus_cat', 1, 'include_0coh', true);
    [ID, mean_coh] = trc.stim_choice(color_level, targ_cho, targ_cor); % correct trials only
    ID(task_set~=2) = NaN; % one task only
    trial_classifier_result(ID, {'color_level', 'targ_cor', 'targ_cho', 'task_set'}, {color_level, targ_cor, targ_cho, task_set});

    % plot PCA
    clear opt
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.epoch = 2;
    opt.PC_range = [250 600];
    opt.PSTH_conv = {'boxcar', 100};
    opt.PC_kernel = {'bartlett', 200};

    fh = showPCA_trajectory(fnd, ID, opt);
    saveas(fh, fullfile(FigDir, sprintf('PCA_trajectory_color_task_session_%d.fig', n)));
    print(fh, '-dpdf', fullfile(FigDir, sprintf('PCA_trajectory_color_task_session_%d.pdf', n)));
end

%% PCA (coh * task)
% color level
% day 1: [96 80]
% day 2: [96 80 60 40 32 24]
% day 3: [96 80 60 40 32 24	18]

for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND')).fnd;

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100; morph_level(task_set==2) = NaN;
    color_level = fnd.getp('color_level')*100; color_level(task_set==1) = NaN;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    ID = task_set;
    ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:));
    ID(targ_cho~=targ_cho) = NaN;
    trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set});

    % plot PCA
    clear opt
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
    opt.epoch = 2;
    opt.PC_range = [250 600];
    opt.PSTH_conv = {'boxcar', 100};
    opt.PC_kernel = {'bartlett', 200};

    fh = showPCA_trajectory(fnd, ID, opt);
    saveas(fh, fullfile(FigDir, sprintf('PCA_trajectory_task-choice_session_%d.fig', n)));
    print(fh, '-dpdf', fullfile(FigDir, sprintf('PCA_trajectory_task-choice_session_%d.pdf', n)));
end


