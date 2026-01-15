run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Woody'; % Nick, Woody
experiment = 'passiveLong'; % learnTask2, learnTask3, learnTask4, faceColor, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'dPCA', 'passiveLong'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'dPCA', 'passiveLong'); mkdir(InterimDir);

%% dPCA of passiveLong main task
type = 'main';
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 1)).fnd;
    classifier = 'pair-choice';
    step = [1 1 1 0]; % regularization, dPCA, projection, plotting
    run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n, type);
end

fh_proj = cell(n_files, 1);
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 1)).fnd;
    classifier = 'pair-choice';
    step = [0 0 0 1]; % regularization, dPCA, projection, plotting
    fh_proj{n} = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n, type);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [5 5], [500 500]*1.5, [-40 40], [-80 80]);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_pair_%s_%s_%s.pdf', monkey, experiment, type)));
print(fh_2, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s_%s.pdf', monkey, experiment, type)));

%% dPCA of passive fixation
type = 'passive';
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 2)).fnd;
    classifier = 'pair-choice';
    step = [1 1 1 0]; % regularization, dPCA, projection, plotting
    run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n, type);
end

fh_proj = cell(n_files, 1);
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted', 2)).fnd;
    classifier = 'pair-choice';
    step = [0 0 0 1]; % regularization, dPCA, projection, plotting
    fh_proj{n} = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n, type);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [5 5], [500 500]*1.5, [-40 40], [-80 80]);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_pair_%s_%s_%s.pdf', monkey, experiment, type)));
print(fh_2, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s_%s.pdf', monkey, experiment, type)));

%% plot overlapped pair signal
list = {'main', 'passive'};
[fh, fh_avg, fh_first_second, fh_pair1_2] = plot_ovelapped_signal(monkey, experiment, n_files, InterimDir, list);

print(fh, '-dpdf', fullfile(FigDir, sprintf('pair_signal_%s_%s.pdf', monkey, experiment)));
print(fh_avg, '-dpdf', fullfile(FigDir, sprintf('average_pair_signal_%s_%s.pdf', monkey, experiment)));
print(fh_first_second, '-dpdf', fullfile(FigDir, sprintf('difference_main_vs_passive_%s_%s.pdf', monkey, experiment)));
print(fh_pair1_2, '-dpdf', fullfile(FigDir, sprintf('difference_pair1_vs_2_%s_%s.pdf', monkey, experiment)));


%% functions
function fh_proj = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, session_id, type)

%% 1. load data/setup parameters
% select epoch
fnd = fnd.extract_epoch(2);

% select unit
r = fnd.FR({1, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% select trial
if strcmp(type, 'main')
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
end

% classifier
coh = fnd.getp('morph_level')*100;
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');
task_set = fnd.getp('task_set');

switch classifier
    case 'pair-choice'
        ID_dpca = task_set; ID_dpca(targ_cor==2) = ID_dpca(targ_cor==2) + max(ID_dpca(:)); % classifier for dPCA
        % trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});
        ID = ID_dpca; % classifier for projection

        param_combination = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        param_name = {'Pair', 'Choice', 'Condition-independent', 'S/C Interaction'};
        target = {{'Pair', 1}, {'Choice', 1}};
end

% options
regularization = true;
tbin = 20; % ms
epoch = 1;

%% 2. determine regularization parameter for dPCA
if step(1) && regularization
    clear opt;
    opt.tbin = tbin; % ms
    opt.t_range = [0 600];
    opt.param_combination = param_combination;
    opt.param_name = param_name;

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_optimize_regularization(r, fnd.tstamp{epoch}, opt);
    print(fh, '-dpdf', fullfile(FigDir, sprintf('step1_regularization_%s_%s_session%d_%s.pdf', monkey, experiment, session_id, type)));
    save(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d_%s.mat', monkey, experiment, session_id, type)), 'data')
end

%% 3. run main dPCA
if step(2)
    clear opt;
    opt.detrend = false;
    opt.tbin = tbin;
    opt.t_range = [0 600];
    opt.param_combination = param_combination;
    opt.param_name = param_name;
    opt.show_figure = false;
    if regularization; opt.regularization = load(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d_%s.mat', monkey, experiment, session_id, type))).data; end

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_stim_choice(r, fnd.tstamp{epoch}, opt);
    if ~isempty(fh); saveas(fh, fullfile(FigDir, sprintf('step2_dPCA_%s_%s_session%d_%s.fig', monkey, experiment, session_id, type)), 'fig'); end
    save(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d_%s.mat', monkey, experiment, session_id, type)), 'data');
end

%% 4. compute neural data along dPCA axes
if step(3)
    clear opt;
    opt.epoch = epoch;
    opt.target = target;
    opt.coefficient_data = load(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d_%s.mat', monkey, experiment, session_id, type))).data;

    data = gen_popresp_dPCA_axis(fnd, ID, opt);
    data.cutoff = [find(fnd.tstamp{opt.epoch}==0), find(fnd.tstamp{opt.epoch}==600)];
    save(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%s.mat', monkey, experiment, session_id, type)), 'data');
end

%% 5. show neural data along dPCA axes
fh_proj = [];
if step(4)
    clear opt;
    opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
    opt.plot = set_plot_opt_2cond('roma', 'roma', max(ID(:))/2);
    switch experiment
        case {'learnTask2', 'passiveLong'}
            opt.plot.color = [0 0 0; 44 145 224; 0 0 0; 44 145 224]/255;
        case 'learnTask3'
            opt.plot.color = [0 0 0; 58 191 153; 0 0 0; 58 191 153]/255;
        case 'learnTask4'
            opt.plot.color = [0 0 0; 240 169 58; 0 0 0; 240 169 58]/255;
        case 'faceColor'
            opt.plot.color = [0 0 0; 1 0 0; 0 0 0; 1 0 0];
    end

    fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%s.mat', monkey, experiment, session_id, type)), opt);
    format_panel(fh_proj, 'ylim', [-55 55])
    print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('step3_projection_%s_%s_session%d_%s.pdf', monkey, experiment, session_id, type)));
end

end
