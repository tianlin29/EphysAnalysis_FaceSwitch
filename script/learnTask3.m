run('../Initialize.m');
monkey = 'Woody';
experiment = 'learnTask3';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

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



