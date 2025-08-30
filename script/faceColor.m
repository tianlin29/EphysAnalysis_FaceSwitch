run('../Initialize.m');
monkey = 'Nick';
experiment = 'faceColor';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'faceColor'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'faceColor'); mkdir(InterimDir);

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

%% PSTH of single units (color task)
n = 3;
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
psth = fnd.PSTH(ID, {'boxcar', 100}, [], [], 2);

clear opt;
opt.epoch = 2;
opt.unitID = fnd.unitID;
opt.plot = set_plot_opt('vik', max(ID(:)));

fh = showSinglePSTH(fnd.tstamp, psth, opt);

for f = 1:length(fh)
    if f==1
        exportgraphics(fh{f}, fullfile(FigDir, 'PSTH_single_unit_color_task.pdf'), 'ContentType', 'vector');
    else
        exportgraphics(fh{f}, fullfile(FigDir, 'PSTH_single_unit_color_task.pdf'), 'ContentType', 'vector', 'Append', true);
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

%% compare two manifolds in a PC space
for n = 3:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'fnd')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor')); % correct trials only

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100;
    color_level = fnd.getp('color_level')*100;
    coh = nan(size(task_set)); coh(task_set==1) = morph_level(task_set==1); coh(task_set==2) = color_level(task_set==2); 
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('stim_group', {[0 24], [24 40], [40 60], [60 80], [80 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(coh, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho', 'task_set'}, {coh, targ_cor, targ_cho, task_set});

    ID_1 = ID; ID_1(task_set~=1) = NaN;
    ID_2 = ID; ID_2(task_set~=2) = NaN;

    % plot manifold
    data_1 = fnd.PSTH(ID_1, {'boxcar', 100});
    data_2 = fnd.PSTH(ID_2, {'boxcar', 100});

    opt.mean_coh = mean_coh;
    opt.plot1 = set_plot_opt('red', max(ID(:)));
    opt.plot2 = set_plot_opt('blue', max(ID(:)));
    opt.PC_range = [250 600];
    opt.Time = [100 200 300 400 500 600];
    opt.epoch = 1;

    opt.roughness = 5e-4;
    opt.dim_sign = [1 -1 1];
    opt.dim = [3 1 2];
    opt.view = [-161 -29];
    opt.PC_target = 'concatenate'; % data1, data2, average, concatenate
    fh = comparePCA(fnd.tstamp, fnd.tstamp, data_1, data_2, opt);

    % change format
    format_panel(gcf, 'fig_size', [1270 250])
    a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',6)
    h = get(gca,'xlabel'); set(h, 'FontSize', 7); h = get(gca,'ylabel'); set(h, 'FontSize', 7); h = get(gca,'zlabel'); set(h, 'FontSize', 7);
    print(fh, '-dpdf', fullfile(FigDir, sprintf('compare_manifold_pair_%s_%s_session%d.pdf', monkey, experiment, n)));
end

%% show view
[currentAzimuth, currentElevation] = view;
disp(['Azimuth: ', num2str(currentAzimuth)]);
disp(['Elevation: ', num2str(currentElevation)]);
