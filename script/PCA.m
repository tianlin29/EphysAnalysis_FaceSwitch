run('../Initialize.m');
monkey = 'Woody';
experiment = 'learnTask3';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'PCA'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'PCA'); mkdir(InterimDir);

%% check trial number and unit number
[nunit, ntask1, ntask2] = deal(nan(n_files, 1));
for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    % select unit
    r = fnd.FR({2, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor')); % correct trials only

    [nunit(n), ~, ~] = size(fnd.data{1});
    task_set = fnd.getp('task_set'); task_set = task_set(1,:);
    ntask1(n) = sum(task_set==1);
    ntask2(n) = sum(task_set==2);
end

tbl = table(nunit, ntask1, ntask2);

%% manifold (coh)
TASK_ID = 2;

for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID); % one task only

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('stim_group', {[0 12], [12 24], [24 40], [40 60], [60 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot manifold
    data = fnd.PSTH(ID, {'boxcar', 100});

    opt.mean_coh = mean_coh;
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.PC_range = [250 600];
    opt.Time = [100 200 300 400 500 600];
    opt.epoch = 1;

    opt.roughness = 5e-4;
    opt.dim_sign = [1 -1 1];
    opt.dim = [3 1 2];
    opt.view = [-145 34];
    fh = showPCA(fnd.tstamp, data, opt);

    % change format
    format_panel(gcf, 'fig_size', [1270 250])
    a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',6)
    h = get(gca,'xlabel'); set(h, 'FontSize', 7); h = get(gca,'ylabel'); set(h, 'FontSize', 7); h = get(gca,'zlabel'); set(h, 'FontSize', 7);
    print(fh, '-dpdf', fullfile(FigDir, sprintf('manifold_pair%d_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)));
end

%% show view
[currentAzimuth, currentElevation] = view;
disp(['Azimuth: ', num2str(currentAzimuth)]);
disp(['Elevation: ', num2str(currentElevation)]);

%% PCA (coh)
TASK_ID = 2;

for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID); % one task only

    if TASK_ID==2
        if strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==1
            fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=60);
        elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==6
            fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=18);
        elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[7 8 9 10])
            fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=12);
        elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[11 12 13 14 15])
            fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=6);
        end

        if strcmp(monkey, 'Woody')
            fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=18);
        end
    end

    if TASK_ID==1
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=18);
    end

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot PCA
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.epoch = 1;
    opt.PC_range = [250 600];
    opt.PSTH_conv = {'boxcar', 100};
    opt.PC_kernel = {'bartlett', 200};

    fh = showPCA_trajectory(fnd, ID, opt);
    % saveas(fh, fullfile(FigDir, sprintf('PCA_coh_pair%d_%s_%s_session%d.fig', TASK_ID, monkey, experiment, n)));
    print(fh, '-dpdf', fullfile(FigDir, sprintf('PCA_coh_pair%d_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)));
end

%% PCA (choice*pair)
for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
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
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    ID = task_set;
    ID(targ_cho==2) = ID(targ_cho==2)+max(ID(:));
    trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set});

    % plot PCA
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
    opt.epoch = 1;
    opt.PC_range = [250 600];
    opt.PSTH_conv = {'boxcar', 100};
    opt.PC_kernel = {'bartlett', 200};

    fh = showPCA_trajectory(fnd, ID, opt);
    % saveas(fh, fullfile(FigDir, sprintf('PCA_choice-pair_%s_%s_session%d.fig', monkey, experiment, n)));
    print(fh, '-dpdf', fullfile(FigDir, sprintf('PCA_choice-pair_%s_%s_session%d.pdf', monkey, experiment, n)));
end

%% single trial trajectories in PC space
for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
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
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    ID = targ_cho;
    trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set});

    % get PCA coef
    data = fnd.PSTH(ID, {'gaussian', 20});

    clear opt;
    opt.plot = set_plot_opt('vik', 2);
    opt.PC_range = [250 600];
    opt.epoch = 1;
    [coef, detPSTH] = getPC_coef(fnd.tstamp, data, opt);

    % show single trial trajectory
    clear opt;
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
    fh = show_single_tr_trajectory(fnd.tstamp{1}, fnd.raster{1}, ID, coef, detPSTH, opt);
    % saveas(fh, fullfile(FigDir, sprintf('single_trial_trajectory_choice_%s_%s_session%d.fig', monkey, experiment, n)));
    print(fh, '-dpdf', fullfile(FigDir, sprintf('single_trial_trajectory_choice_%s_%s_session%d.pdf', monkey, experiment, n)));
end

%% compare two manifolds in a PC space
for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'fnd_sorted')).fnd;
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
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('stim_group', {[0 24], [24 40], [40 60], [60 80], [80 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

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
    opt.view = [-48 26];
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
