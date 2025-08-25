run('../Initialize.m');
monkey = 'Nick';
experiment = 'learnTask4';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

task_color = [99 97 172; 178 34 34]/255;
subj_list.Woody = {'Woody20240110', 'Woody20240111', 'Woody20240112', 'Woody20240113', 'Woody20240115'}; % 5 sessions
subj_list.Nick = {'Nick20250721', 'Nick20250722', 'Nick20250723', 'Nick20250724', 'Nick20250725', ...
    'Nick20250726', 'Nick20250727', 'Nick20250728'}; % 8 sessions

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
% process data
opt.log = true;
opt.constant = false;
opt.verbose = false;
% plot
opt.color = task_color;
opt.average = false; % do not average across learning sessions
opt.legend = {'Pair 1', 'Pair 4'};

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond, coh, cor, date, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
save(fullfile(InterimDir, 'unsigned_choice.mat'), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, 'unsigned_choice_summary.pdf'));
print(fh_all, '-dpdf', fullfile(FigDir, 'unsigned_choice.pdf'));

%% [stat] accuracy of unsigned choice
stat = load(fullfile(InterimDir, 'unsigned_choice.mat')).stat;
nses = length(stat);
acc1 = cellfun(@(x) x.acc1, stat);
acc2 = cellfun(@(x) x.acc2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', acc1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Accuracy of pair 1 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Accuracy of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
else
    error('check')
end
[rho, p_val] = corr((1:nses)', acc2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Accuracy of pair 2 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Accuracy of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
else
    error('check')
end

%% [stat] threshold of unsigned choice
stat = load(fullfile(InterimDir, 'unsigned_choice.mat')).stat;
nses = length(stat);
thres1 = cellfun(@(x) x.thres1, stat);
thres2 = cellfun(@(x) x.thres2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', thres1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Threshold of pair 1 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Threshold of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Threshold of pair 1 decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end
[rho, p_val] = corr((1:nses)', thres2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Threshold of pair 2 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Threshold of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Threshold of pair 2 decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

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



