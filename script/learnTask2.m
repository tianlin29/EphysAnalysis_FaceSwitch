run('../Initialize.m');
monkey = 'Nick';
experiment = 'learnTask2';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

task_color = [99 97 172; 254 160 64]/255;
subj_list.Nick = {'Nick20250528', 'Nick20250529', 'Nick20250530', 'Nick20250531', 'Nick20250601', ...
    'Nick20250604', 'Nick20250605', 'Nick20250606', 'Nick20250607', 'Nick20250608', ...
    'Nick20250609', 'Nick20250610', 'Nick20250611', 'Nick20250612', 'Nick20250613'}; % 15 sessions

%%
D = combineTrialData(fullfile(PreprocDir, sprintf('%s_%s.mat', monkey, experiment)));

%% unsigned choice
date = cellfun(@(x) x.date, D, 'uni', 0);
cond = cellfun(@(x) x.cond, D);
coh = cellfun(@(x) x.coh, D);
targ_cor = cellfun(@(x) x.targ_cor, D);
resp = cellfun(@(x) x.resp, D);
cor = resp==targ_cor;

clear opt
% select data
opt.subj_list = subj_list.(monkey);
% opt.subj_list = opt.subj_list(end-2:end-1);
% process data
opt.log = true;
opt.constant = false;
opt.verbose = false;
% plot
opt.color = task_color;
opt.average = false; % do not average across learning sessions
opt.legend = {'Pair 1', 'Pair 2'};

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond, coh, cor, date, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
save(fullfile(InterimDir, 'unsigned_choice.mat'), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, 'unsigned_choice_summary.pdf'));
print(fh_all, '-dpdf', fullfile(FigDir, 'unsigned_choice.pdf'));

%% [stat] unsigned choice
stat = load(fullfile(InterimDir, 'unsigned_choice.mat')).stat;
nses = length(stat);
acc1 = cellfun(@(x) x.acc1, stat);
acc2 = cellfun(@(x) x.acc2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', acc1', 'Type', 'Spearman');
fprintf('Accuracy of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
[rho, p_val] = corr((1:nses)', acc2', 'Type', 'Spearman');
fprintf('Accuracy of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))

%% flow field (task 1 or 2)
TASK_ID = 1;
VERBOSE = false;

fh_idv = cell(n_files, 1);
for n = 1:n_files
    % load file
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

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








