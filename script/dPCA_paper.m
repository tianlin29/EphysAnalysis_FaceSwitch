run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Woody'; % Nick, Woody
experiment = 'noCue'; % learnTask2, learnTask3, learnTask4, faceColor, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'dPCA', 'pair-choice'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'dPCA', 'pair-choice'); mkdir(InterimDir);

%% dPCA
for n = 1
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted_preprocessed')).fnd;
    classifier = 'pair-choice';
    epoch = 3;
    t_range = [-400 200]; % [0 600; -200 400; -400 200; -200 200]
    step = [1 1 1 1]; % regularization, dPCA, projection, plotting
    run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n, epoch, t_range);
end

%% plot 3 epochs
n = 1;
t_range = [0 600; -200 400; -400 200];
t_cutoff = [-100 600; -200 400; -400 200];

[t_all, pair_all, choice_all] = deal(cell(1,3));
for i = 1:3
    data = load(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%d_%d.mat', monkey, experiment, n, t_range(i,1), t_range(i,2)))).data;
    
    for c = 1:4
        t = data.tstamp;
        x = data.dpc(:,:,c);
        I = t>=t_cutoff(i,1) & t<=t_cutoff(i,2);
        t = t(I); x = x(:,I);

        if i==3 ; x(2,:) = -x(2,:); end % flip sign
        x = nanconv(x, fspecial('average', [1, 100]), 'same');

        t_idx = round(linspace(1, length(t), 50));
        t = t(t_idx);
        x = x(:,t_idx);
        pair_all{i}(c,:) = x(1,:);
        choice_all{i}(c,:) = x(2,:);
    end
    t_all{i} = t; 
end

clear opt
opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
opt.plot.color = [0 0 0; 44 145 224; 0 0 0; 44 145 224]/255;
opt.event_label = {'Stim on', 'Stim off', 'FP off'};
opt.ylim = [-20 20];

figure('Position', [100 100 250 150]);
ax = plot_multiepoch_trace(t_all, pair_all, [], opt);
axes(ax);
xlabel('Time (ms)');
ylabel('Pair signal');
print(gcf, '-dpdf', fullfile(FigDir, sprintf('pair_signal_throughout_trial.pdf')))

opt.ylim = [-100 100];
figure('Position', [100 400 250 150]);
ax = plot_multiepoch_trace(t_all, choice_all, [], opt);
axes(ax);
xlabel('Time (ms)');
ylabel('Choice signal');
allAxes = findobj(gcf, 'type', 'axes');
set(allAxes(3), 'YColor', 'none');
print(gcf, '-dpdf', fullfile(FigDir, sprintf('choice_signal_throughout_trial.pdf')))







%% functions
function fh_proj = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, session_id, epoch, t_range)

%% 1. load data/setup parameters
% select epoch
fnd = fnd.extract_epoch(epoch);

% select unit
% r = fnd.FR({1, [100 500]});
% I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
% fnd = fnd.set_unit_criteria('custom', I);

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
        ID = ID_dpca; % classifier for projection

        param_combination = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
        param_name = {'Coh', 'Choice', 'Condition-independent', 'S/C Interaction'};
        target = {{'Coh', 1}, {'Choice', 1}};

    case 'pair-choice'
        ID_dpca = task_set; ID_dpca(targ_cho==2) = ID_dpca(targ_cho==2) + max(ID_dpca(:)); % classifier for dPCA
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
    opt.t_range = t_range;
    opt.param_combination = param_combination;
    opt.param_name = param_name;

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_optimize_regularization(r, fnd.tstamp{epoch}, opt);
    print(fh, '-dpdf', fullfile(FigDir, sprintf('step1_regularization_%s_%s_session%d_%d_%d.pdf', monkey, experiment, session_id, t_range(1), t_range(2))));
    save(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d_%d_%d.mat', monkey, experiment, session_id, t_range(1), t_range(2))), 'data')
end

%% 3. run main dPCA
if step(2)
    clear opt;
    opt.detrend = false;
    opt.tbin = tbin;
    opt.t_range = t_range;
    opt.param_combination = param_combination;
    opt.param_name = param_name;
    opt.show_figure = false;
    if regularization; opt.regularization = load(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d.mat', monkey, experiment, session_id))).data; end

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_stim_choice(r, fnd.tstamp{epoch}, opt);
    if ~isempty(fh); saveas(fh, fullfile(FigDir, sprintf('step2_dPCA_%s_%s_session%d_%d_%d.fig', monkey, experiment, session_id, t_range(1), t_range(2))), 'fig'); end
    save(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d_%d_%d.mat', monkey, experiment, session_id, t_range(1), t_range(2))), 'data');
end

%% 4. compute neural data along dPCA axes
if step(3)
    clear opt;
    opt.epoch = epoch;
    opt.target = target;
    opt.coefficient_data = load(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d_%d_%d.mat', monkey, experiment, session_id, t_range(1), t_range(2)))).data;

    data = gen_popresp_dPCA_axis(fnd, ID, opt);
    data.cutoff = [find(fnd.tstamp{opt.epoch}==t_range(1)), find(fnd.tstamp{opt.epoch}==t_range(2))];
    save(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%d_%d.mat', monkey, experiment, session_id, t_range(1), t_range(2))), 'data');
end

%% 5. show neural data along dPCA axes
fh_proj = [];
if step(4)
    clear opt;
    opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
    opt.plot = set_plot_opt_2cond('roma', 'roma', max(ID(:))/2);
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

    fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%d_%d.mat', monkey, experiment, session_id, t_range(1), t_range(2))), opt);
    format_panel(fh_proj, 'ylim', [-55 55])
    print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('step3_projection_%s_%s_session%d_%d_%d.pdf', monkey, experiment, session_id, t_range(1), t_range(2))));
end

end


