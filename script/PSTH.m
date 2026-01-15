run('../Initialize.m');
monkey = 'Woody'; % Nick, Woody
experiment = 'noCue'; % learnTask2, learnTask2_increaseSwitchFreq, learnTask3, learnTask4
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'PSTH'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'PSTH'); mkdir(InterimDir);

%% save raster for RNN fitting
fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;

% select unit
r = fnd.FR({2, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% select trial
fnd = fnd.extract_trial(fnd.getp('task_set')==1); % one task only
fnd = fnd.extract_trial(fnd.getp('targ_cor')==fnd.getp('targ_cho'));

raster = fnd.raster(2); raster = raster{1}; % (unit, time, trial)
save('D:\Engine_D2\latent_dynamics\raster.mat', 'raster');

targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:);
save('D:\Engine_D2\latent_dynamics\targ_cho.mat', 'targ_cho');


%% PSTH (coh)
TASK_ID = 1;
DETREND = true;

for n = 1:9
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    % select unit
    r = fnd.FR({2, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID); % one task only

    if strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==1
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=60);
    elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==6
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=18);
    elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[7 8 9 10])
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=12);
    elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[11 12 13 14 15])
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=6);
    end

    if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask3')
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=6);
    end

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot PSTH
    psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], DETREND);

    clear opt
    opt.cutoff = fnd.cutoff();
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.legend = arrayfun(@(x) num2str(x, '%.0f'), mean_coh, 'uni', 0);
    opt.legend_pos = [0.01 0.25 0.1 0.5];
    opt.event_label = {fnd.alignto.event};
    fh = showPopPSTH(fnd.tstamp, psth_data, opt);
    if DETREND
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_coh_pair%d_detrended_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)));
    else
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_coh_pair%d_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)));
    end
end

%% PSTH (choice)
TASK_ID = 1;
DETREND = true;

for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    % select unit
    r = fnd.FR({2, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID); % one task only

    if strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==1
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=60);
    elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==6
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=18);
    elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[7 8 9 10])
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=12);
    elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[11 12 13 14 15])
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=6);
    end

    if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask3')
        fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=6);
    end

    % get ID
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('stim_group', {[0 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot PSTH
    psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], DETREND);

    clear opt
    opt.cutoff = fnd.cutoff();
    opt.plot = set_plot_opt('roma', max(ID(:)));
    opt.legend = {'Choice 1', 'Choice 2'};
    opt.legend_pos = [0.01 0.25 0.1 0.5];
    opt.event_label = {fnd.alignto.event};
    fh = showPopPSTH(fnd.tstamp, psth_data, opt);
    if DETREND
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_choice_pair%d_detrended_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)));
    else
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_choice_pair%d_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)));
    end
end


%% PSTH (choice*pair)
DETREND = true;

for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    % select unit
    r = fnd.FR({2, [100 500]});
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

    % plot PSTH
    psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], DETREND);

    clear opt
    opt.cutoff = fnd.cutoff();
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
    opt.legend = {'Pair1 Choice1', 'Pair2 Choice1', 'Pair1 Choice2', 'Pair2 Choice2'};
    opt.legend_pos = [0.01 0.25 0.1 0.5];
    opt.event_label = {fnd.alignto.event};
    fh = showPopPSTH(fnd.tstamp, psth_data, opt);
    if DETREND
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_choice-pair_detrended_%s_%s_session%d.pdf', monkey, experiment, n)));
    else
        print(fh, '-dpdf', fullfile(FigDir, sprintf('PSTH_choice-pair_%s_%s_session%d.pdf', monkey, experiment, n)));
    end
end

%% example units (coh)
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted_preprocessed')).fnd;
unitID = [1 1 1;
          1 17 1];
trc = trial_classifier('stim_group', {[0 24], [24 48], [48 Inf]}, 'plus_cat', 1);

clear opt;
opt.epoch = [1 2 4];
opt.less_timepoints = 1;
opt.cutoff = [{[1 551]}, {[51 451]}, {[1 696]}, {[51 401]}];
opt.plot = set_plot_opt('vik', 6);
opt.plot.linewidth = 0.75 * ones(6,1);
opt.event_label = {'Stim on', 'Stim off', 'FP off', 'Resp'};
for n = 1:size(unitID, 1)
    fnd1 = fnd.set_unit_criteria('customID', unitID(n,:));
    ID = trc.stim_choice(fnd1.getp('morph_level')*100, fnd1.getp('targ_cho')); % sort by coh
    psth_data = fnd1.PSTH(ID, {'boxcar', 100});
    fh = showPopPSTH(fnd1.tstamp, psth_data, opt);
    title(sprintf('channel%d unit%d', unitID(n,2), unitID(n,3)))
    set(findobj(gcf, 'Type', 'axes'), 'XTickLabelRotation', 90);
    allAxes = findobj(gcf, 'type', 'axes');
    set(allAxes(3), 'YColor', 'none');
    set(gcf, 'pos', [100 100 245 85])  
    print(fh, '-dpdf', fullfile(FigDir, ['PSTH_coh' sprintf('_ch%d_u%d', unitID(n,2), unitID(n,3)) '.pdf']));
end

%% example units (pair)
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted_preprocessed')).fnd;
unitID = [1 1 2;
          1 10 1;
          1 53 3;
          1 54 1;
          1 80 2;
          1 86 1;
          1 93 2];
trc = trial_classifier('stim_group', {[0 Inf]}, 'plus_cat', 1);

clear opt;
opt.epoch = [1 2 4];
opt.less_timepoints = 1;
opt.cutoff = [{[1 551]}, {[51 451]}, {[1 696]}, {[51 401]}];
opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
switch experiment
    case {'learnTask2', 'faceColor', 'threeExemplar', 'noCue', 'learnTask2_increaseSwitchFreq'}
        new_color_set = [0 0 0; 44 145 224; 0 0 0; 44 145 224]/255;
    case 'learnTask3'
        new_color_set = [0 0 0; 58 191 153; 0 0 0; 58 191 153]/255;
    case 'learnTask4'
        new_color_set = [0 0 0; 240 169 58; 0 0 0; 240 169 58]/255;
end
opt.plot.color = new_color_set;
opt.plot.facecolor = new_color_set;
opt.plot.edgecolor = new_color_set;
opt.plot.linewidth = 0.75 * ones(6,1);
opt.event_label = {'Stim on', 'Stim off', 'FP off', 'Resp'};
for n = 1:size(unitID, 1)
    fnd1 = fnd.set_unit_criteria('customID', unitID(n,:));
    targ_cho = fnd1.getp('targ_cho');
    task_set = fnd1.getp('task_set');
    ID = task_set; ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:));
    psth_data = fnd1.PSTH(ID, {'boxcar', 100});
    fh = showPopPSTH(fnd1.tstamp, psth_data, opt);
    title(sprintf('channel%d unit%d', unitID(n,2), unitID(n,3)))
    set(findobj(gcf, 'Type', 'axes'), 'XTickLabelRotation', 90);
    allAxes = findobj(gcf, 'type', 'axes');
    set(allAxes(3), 'YColor', 'none');
    set(gcf, 'pos', [100 100 245 85])  
    print(fh, '-dpdf', fullfile(FigDir, ['PSTH_pair' sprintf('_ch%d_u%d', unitID(n,2), unitID(n,3)) '.pdf']));
end

%% PSTH of single units (coh*pair)
TASK_ID = 2;

for n = 10:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = extract_epoch(fnd, 2);

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

    trc = trial_classifier('stim_group', {[0 24], [24 48], [48 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % plot PSTH
    psth = fnd.PSTH(ID, {'boxcar', 100});

    clear opt;
    opt.unitID = fnd.unitID;
    opt.plot = set_plot_opt('vik', max(ID(:)));

    fh = showSinglePSTH(fnd.tstamp, psth, opt);

    for f = 1:length(fh)
        if f==1
            exportgraphics(fh{f}, fullfile(FigDir, sprintf('PSTH_singleUnit_coh_pair%d_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)), 'ContentType', 'vector');
        else
            exportgraphics(fh{f}, fullfile(FigDir, sprintf('PSTH_singleUnit_coh_pair%d_%s_%s_session%d.pdf', TASK_ID, monkey, experiment, n)), 'ContentType', 'vector', 'Append', true);
        end
    end
end

%% PSTH of single units (choice*pair)
for n = n_files:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = extract_epoch(fnd, 2);

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

    % plot PSTH
    psth = fnd.PSTH(ID, {'boxcar', 100});

    clear opt;
    opt.unitID = fnd.unitID;
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);

    fh = showSinglePSTH(fnd.tstamp, psth, opt);

    for f = 1:length(fh)
        if f==1
            exportgraphics(fh{f}, fullfile(FigDir, sprintf('PSTH_singleUnit_choice-pair_%s_%s_session%d.pdf', monkey, experiment, n)), 'ContentType', 'vector');
        else
            exportgraphics(fh{f}, fullfile(FigDir, sprintf('PSTH_singleUnit_choice-pair_%s_%s_session%d.pdf', monkey, experiment, n)), 'ContentType', 'vector', 'Append', true);
        end
    end
end
