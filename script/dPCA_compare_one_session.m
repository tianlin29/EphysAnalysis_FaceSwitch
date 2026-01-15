%% initialize dPCA
run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask3'; % learnTask2, learnTask3, learnTask4, faceColor, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'dPCA', 'compare_one_session'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'dPCA', 'compare_one_session'); mkdir(InterimDir);

%% dPCA
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    step = [1 1 1 0]; % regularization, dPCA, projection, plotting
    run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n);
end

%% plot all sessions
[fh_proj_1, fh_proj_2] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    step = [0 0 0 1]; % regularization, dPCA, projection, plotting
    [fh_proj_1{n}, fh_proj_2{n}] = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n);
end

[fh_1_pair, fh_1_choice] = plot_in_one_fig_dPCA(fh_proj_1, [5 5], [500 500]*1.5, [-50 50], [-50 50]);
[fh_2_pair, fh_2_choice] = plot_in_one_fig_dPCA(fh_proj_2, [5 5], [500 500]*1.5, [-50 50], [-50 50]);

print(fh_1_pair, '-dpdf', fullfile(FigDir, sprintf('dPCA_pair_%s_%s_first.pdf', monkey, experiment)));
print(fh_2_pair, '-dpdf', fullfile(FigDir, sprintf('dPCA_pair_%s_%s_second.pdf', monkey, experiment)));
print(fh_1_choice, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s_first.pdf', monkey, experiment)));
print(fh_2_choice, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s_second.pdf', monkey, experiment)));

%% plot overlapped pair signal
list = {'first', 'second'};
[fh, fh_avg, fh_first_second, fh_pair1_2] = plot_ovelapped_signal(monkey, experiment, n_files, InterimDir, list);

print(fh, '-dpdf', fullfile(FigDir, sprintf('pair_signal_%s_%s.pdf', monkey, experiment)));
print(fh_avg, '-dpdf', fullfile(FigDir, sprintf('average_pair_signal_%s_%s.pdf', monkey, experiment)));
print(fh_first_second, '-dpdf', fullfile(FigDir, sprintf('difference_first_vs_second_%s_%s.pdf', monkey, experiment)));
print(fh_pair1_2, '-dpdf', fullfile(FigDir, sprintf('difference_pair1_vs_2_%s_%s.pdf', monkey, experiment)));



%% functions
function [fh_proj_1, fh_proj_2] = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, session_id)

%% 1. load data/setup parameters
% select epoch
fnd = fnd.extract_epoch(2);

% select unit
[nunit, ntime, ntrial] = size(fnd.data{1});
r = fnd.FR({1, [100 500]});
I = nanmean(r, 2)>=1 & nanmean(r, 2)<50; % mean FR should >= 1 Hz

r_1 = mean(r(:,1:round(ntrial/2)), 2);
r_2 = mean(r(:,round(ntrial/2):end), 2);
I2 = abs((r_2 - r_1) ./ r_1) < 0.5;

fnd = fnd.set_unit_criteria('custom', I&I2);

% select trial
fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));

% classifier
coh = fnd.getp('morph_level')*100;
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');
task_set = fnd.getp('task_set');

switch classifier
    case 'coh-choice'
        trc = trial_classifier('stim_group', {[0 18], [18 32], [32 60], [60 80], [80 Inf]}, 'plus_cat', 1, 'include_0coh', false);
        ID_dpca = trc.stim_choice(coh, targ_cho, targ_cor); % correct trials only
        trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});

        param_combination = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        param_name = {'Coh', 'Choice', 'Condition-independent', 'S/C Interaction'};
        target = {{'Coh', 1}, {'Choice', 1}};

    case 'pair-choice'
        ID_dpca = task_set; ID_dpca(targ_cho==2) = ID_dpca(targ_cho==2) + max(ID_dpca(:)); % classifier for dPCA
        % trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});

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
    print(fh, '-dpdf', fullfile(FigDir, sprintf('step1_regularization_%s_%s_session%d.pdf', monkey, experiment, session_id)));
    save(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data')
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
    if regularization; opt.regularization = load(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d.mat', monkey, experiment, session_id))).data; end

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_stim_choice(r, fnd.tstamp{epoch}, opt);
    if ~isempty(fh); saveas(fh, fullfile(FigDir, sprintf('step2_dPCA_%s_%s_session%d.fig', monkey, experiment, session_id)), 'fig'); end
    save(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data');
end

%% 4. compute neural data along dPCA axes
if step(3)
    clear opt;
    opt.epoch = epoch;
    opt.target = target;
    opt.coefficient_data = load(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d.mat', monkey, experiment, session_id))).data;

    trial = fnd.getp('trial'); idx = max(trial(:))/2;
    fnd_1 = fnd.extract_trial(trial<idx);
    fnd_2 = fnd.extract_trial(trial>=idx);

    switch classifier
        case 'pair-choice'
            targ_cho = fnd_1.getp('targ_cho');
            task_set = fnd_1.getp('task_set');
            ID_1 = task_set; ID_1(targ_cho==2) = ID_1(targ_cho==2) + max(ID_1(:)); % classifier for projection

            targ_cho = fnd_2.getp('targ_cho');
            task_set = fnd_2.getp('task_set');
            ID_2 = task_set; ID_2(targ_cho==2) = ID_2(targ_cho==2) + max(ID_2(:)); % classifier for projection
    end

    data = gen_popresp_dPCA_axis(fnd_1, ID_1, opt);
    data.cutoff = [find(fnd_1.tstamp{opt.epoch}==0), find(fnd_1.tstamp{opt.epoch}==600)];
    save(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_first.mat', monkey, experiment, session_id)), 'data');

    data = gen_popresp_dPCA_axis(fnd_2, ID_2, opt);
    data.cutoff = [find(fnd_2.tstamp{opt.epoch}==0), find(fnd_2.tstamp{opt.epoch}==600)];
    save(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_second.mat', monkey, experiment, session_id)), 'data');
end

%% 5. show neural data along dPCA axes
fh_proj = [];
if step(4)
    clear opt;
    opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
    switch experiment
        case 'learnTask2'
            opt.plot.color = [0 0 0; 44 145 224; 0 0 0; 44 145 224]/255;
        case 'learnTask3'
            opt.plot.color = [0 0 0; 58 191 153; 0 0 0; 58 191 153]/255;
        case 'learnTask4'
            opt.plot.color = [0 0 0; 240 169 58; 0 0 0; 240 169 58]/255;
        case 'faceColor'
            opt.plot.color = [0 0 0; 1 0 0; 0 0 0; 1 0 0];
    end

    fh_proj_1 = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_first.mat', monkey, experiment, session_id)), opt);
    format_panel(fh_proj_1, 'ylim', [-55 55])
    print(fh_proj_1, '-dpdf', fullfile(FigDir, sprintf('step3_projection_%s_%s_session%d_first.pdf', monkey, experiment, session_id)));

    fh_proj_2 = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_second.mat', monkey, experiment, session_id)), opt);
    format_panel(fh_proj_2, 'ylim', [-55 55])
    print(fh_proj_2, '-dpdf', fullfile(FigDir, sprintf('step3_projection_%s_%s_session%d_second.pdf', monkey, experiment, session_id)));
end

end


