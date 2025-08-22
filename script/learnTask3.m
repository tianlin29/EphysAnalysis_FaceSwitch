run('../Initialize.m');
monkey = 'Nick';
experiment = 'learnTask3';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

%% PSTH of single units (pair 1)
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

% get ID
task_set = fnd.getp('task_set');
morph_level = fnd.getp('morph_level')*100;
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');

trc = trial_classifier('plus_cat', 1, 'include_0coh', true);
[ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
ID(task_set~=1) = NaN; % one task only
trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

% plot PSTH
psth = fnd.PSTH(ID, {'boxcar', 100}, [], [], 2);

clear opt;
opt.epoch = 2;
opt.unitID = fnd.unitID;
opt.plot = set_plot_opt('vik', max(ID(:)));

fh = showSinglePSTH(fnd.tstamp, psth, opt);

for f = 1:length(fh)
    if f==1
        exportgraphics(fh{f}, fullfile(FigDir, 'PSTH_single_unit_pair1.pdf'), 'ContentType', 'vector');
    else
        exportgraphics(fh{f}, fullfile(FigDir, 'PSTH_single_unit_pair1.pdf'), 'ContentType', 'vector', 'Append', true);
    end
end

%% comparison of PSTH between two choices [useless?]
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

% preprocess
trialFlg = fnd.getp('task_set')==1;
fnd = fnd.extract_trial(trialFlg);

% get ID
trc = trial_classifier('stim_group', {[0 12], [12 24], [24 48], [48 Inf]}, 'plus_cat', 1);

[~, mean_coh] = trc.stim_category(fnd.getp('morph_level')*100, fnd.getp('targ_cor'));
ID1 = trc.stim_plus_choice(fnd.getp('morph_level')*100, fnd.getp('targ_cho'), fnd.getp('targ_cor'));
ID2 = trc.stim_minus_choice(fnd.getp('morph_level')*100, fnd.getp('targ_cho'), fnd.getp('targ_cor'));

% plot PSTH
psth_data1 = fnd.PSTH(ID1, {'boxcar', 100}, [], [], 2);
psth_data2 = fnd.PSTH(ID2, {'boxcar', 100}, [], [], 2);

clear opt;
opt.cutoff = fnd.cutoff();
opt.plot = set_plot_opt('roma', 2);
opt.epoch = 2;
opt.title = arrayfun(@(x)sprintf('coh=%1.0f', x), mean_coh, 'uni', 0);
fh = showPopPSTH_choice_dependent(fnd.tstamp, psth_data1, psth_data2, opt);
print(fh, '-dpdf', fullfile(FigDir, 'PopPSTH_choice.pdf'));

%% population average PSTH (align preference) [useless?]
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

% preprocess
trialFlg = fnd.getp('task_set')==1;
fnd = fnd.extract_trial(trialFlg);

% get preference
    % get stimulus preference of each unit
    % if the target is targ_cho instead of targ_cor, it becomes choice
    % preference
clear opt;
opt.epoch = 2;
opt.t_range = [100 200];
opt.require_significance = false;
[stim_pref, pval] = get_unit_preference(fnd, fnd.getp('targ_cor'), opt);

% plot PSTH
trc = trial_classifier('stim_group', {[0 15], [15 30], [30 60], [60 Inf]}, 'plus_cat', 1);
ID = trc.stim_choice(fnd.getp('morph_level')*100, fnd.getp('targ_cho')); % correct trials only

psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], [], 2);

clear opt;
cutoff = fnd.cutoff();
opt.cutoff = cutoff(2);
opt.plot = set_plot_opt('vik', max(ID(:)));
opt.flip_preference = stim_pref==2; % units prefering category 2 will be reversed
fh = showPopPSTH(fnd.tstamp(2), psth_data(2), opt);
print(fh, '-dpdf', fullfile(FigDir, 'PopPSTH_aligned.pdf'));

%% flow field (task 1 or 2)
TASK_ID = 1;
VERBOSE = false;
PREPROCESS_FND = false;

fh_idv = cell(n_files, 1);
for n = 1:1
    % load file
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    if PREPROCESS_FND
        % remove low FR units
        FR = fnd.FR({2, [100 400]}); % (unit, trial)
        I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);
        fprintf('%d low FR units are removed.\n', sum(~I))

        % remove units changed over time
        I = ~isnan(fnd.misc.SNR);
        fnd = fnd.set_unit_criteria('custom', I);
        fprintf('%d units changed over time.\n', sum(~I))
    end

    % get PC space to plot
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    ID = targ_cho; ID(targ_cho~=targ_cor) = NaN; % correct trials only
    ID(task_set~=TASK_ID) = NaN; % one task only
    if VERBOSE; trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set}); end

    clear opt;
    opt.plot = set_plot_opt('roma', 2);
    opt.PC_range = [250 600];
    opt.conv_kernel = fspecial('average', [1 100]);
    opt.epoch = 2;

    data = fnd.PSTH(ID, {'gaussian', 20}, [], [], opt.epoch);
    [coef, detPSTH] = getPC_coef(fnd.tstamp, data, opt);

    % calculate flow field
    clear opt;
    opt.PC_coef = coef;
    opt.detPSTH = detPSTH;
    opt.dim = [1 2];
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
    opt.bin_num = 20;

    data = calc_flow_field(fnd.tstamp{2}, fnd.raster{2}, ID, opt);

    % show flow field (PC 1-2)
    clear pltopt;
    pltopt.min_sample = 1000;
    pltopt.trj_color = [
        0.00  0.45  0.74;
        0.85  0.33  0.10];
    pltopt.legend = {'Choice 1', 'Choice 2'};
    pltopt.xlabel = 'PC 1';
    pltopt.ylabel = 'PC 2';
    pltopt.title = sprintf('Task %d', TASK_ID);

    fh_idv{n} = show_flow_field(data, pltopt);
    print(fh_idv{n}, '-dpdf', fullfile(FigDir, sprintf('flow_field_PC1-2_pair%d_session%d.pdf', TASK_ID, n)));
end

fh_all = plot_in_one_fig_flow_field(fh_idv, [5 3], [300 500*1.2]*2);
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('flow_field_PC1-2_pair%d.pdf', TASK_ID)));



