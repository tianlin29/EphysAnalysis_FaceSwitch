run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask4'; % learnTask2, learnTask3, learnTask4, faceColor, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'dPCA', 'one_category'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'dPCA', 'one_category'); mkdir(InterimDir);

%% dPCA
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    classifier = 'category_2';
    step = [1 1 1 1]; % regularization, dPCA, projection, plotting
    run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n);
end

%% plot overlapped signal
% average signal
list = {'category_1', 'category_2'};
[task_signal] = deal(nan(801, 2, n_files, 2)); % (time, 2 cond, session, first/second)
for t = 1:2
    for n = 1:n_files
        data = load(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%s.mat', monkey, experiment, n, list{t}))).data;
        dpc = data.dpc; % (dim, time, cond) 2task*2choice
        tstamp = data.tstamp; % (1, time)

        % smooth dPC score
        for d = 1:size(dpc, 1)
            for c = 1:size(dpc, 3)
                dpc(d,:,c) = nanconv(dpc(d,:,c), fspecial('average', [1,100]), 'same');
            end
        end

        % calculate task signal difference in dPC score
        dpc_mn = squeeze(dpc);
        mn = mean(dpc_mn, 1);
        if mn(1)>mn(2) % task 1 is negative
            dpc_mn = -dpc_mn;
        end
        task_signal(:,:,n,t) = dpc_mn;
    end
end

% prepare plotting options
switch experiment
    case {'learnTask2', 'passiveLong'}
        opt.plot.color = [0 0 0; 44 145 224]/255;
    case 'learnTask3'
        opt.plot.color = [0 0 0; 58 191 153]/255;
    case 'learnTask4'
        opt.plot.color = [0 0 0; 240 169 58]/255;
    case 'faceColor'
        opt.plot.color = [0 0 0; 1 0 0];
end
opt.plot.linewidth = [1 1];
opt.plot.linestyle = {'-', '--'};

% plot overlapped signal of 2task*first/second
fh = figure('Position', [50 100 900 900]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for t = 1:2
        for c = 1:2
            plot(tstamp, task_signal(:,c,n,t), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{t}, 'LineWidth', 1)
        end
    end
    plot(xlim, [0 0], ':', 'Color', 'black')
    title(sprintf('session %d', n))
end
format_panel(gcf, 'xlabel', 'Time (ms)', 'ylabel', 'Pair signal', 'ylim', [-60 60])





%% functions
function fh_proj = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, session_id)

%% 1. load data/setup parameters
% select epoch
fnd = fnd.extract_epoch(2);

% select unit
r = fnd.FR({1, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% select trial
fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));

% classifier
coh = fnd.getp('morph_level')*100;
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');
task_set = fnd.getp('task_set');

switch classifier
    case 'category_1'
        ID_dpca = task_set; ID_dpca(targ_cho~=1) = NaN; % classifier for dPCA
        % trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});
        ID = ID_dpca; % classifier for projection

        param_combination = {{1, [1 2]}, {2}};
        param_name = {'Pair', 'Condition-independent'};
        target = {{'Pair', 1}};

    case 'category_2'
        ID_dpca = task_set; ID_dpca(targ_cho~=2) = NaN; % classifier for dPCA
        % trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});
        ID = ID_dpca; % classifier for projection

        param_combination = {{1, [1 2]}, {2}};
        param_name = {'Pair', 'Condition-independent'};
        target = {{'Pair', 1}};
end

% options
regularization = false;
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
    print(fh, '-dpdf', fullfile(FigDir, sprintf('step1_regularization_%s_%s_session%d_%s.pdf', monkey, experiment, session_id, classifier)));
    save(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d_%s.mat', monkey, experiment, session_id, classifier)), 'data')
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
    if regularization; opt.regularization = load(fullfile(InterimDir, sprintf('step1_regularization_%s_%s_session%d_%s.mat', monkey, experiment, session_id, classifier))).data; end

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_stim_choice(r, fnd.tstamp{epoch}, opt);
    if ~isempty(fh); saveas(fh, fullfile(FigDir, sprintf('step2_dPCA_%s_%s_session%d_%s.fig', monkey, experiment, session_id, classifier)), 'fig'); end
    save(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d_%s.mat', monkey, experiment, session_id, classifier)), 'data');
end

%% 4. compute neural data along dPCA axes
if step(3)
    clear opt;
    opt.epoch = epoch;
    opt.target = target;
    opt.coefficient_data = load(fullfile(InterimDir, sprintf('step2_dPCA_%s_%s_session%d_%s.mat', monkey, experiment, session_id, classifier))).data;

    data = gen_popresp_dPCA_axis(fnd, ID, opt);
    data.cutoff = [find(fnd.tstamp{opt.epoch}==0), find(fnd.tstamp{opt.epoch}==600)];
    save(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%s.mat', monkey, experiment, session_id, classifier)), 'data');
end

%% 5. show neural data along dPCA axes
fh_proj = [];
if step(4)
    clear opt;
    opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
    opt.plot = set_plot_opt('roma', 2);
    switch experiment
        case 'learnTask2'
            opt.plot.color = [0 0 0; 44 145 224]/255;
        case 'learnTask3'
            opt.plot.color = [0 0 0; 58 191 153]/255;
        case 'learnTask4'
            opt.plot.color = [0 0 0; 240 169 58]/255;
        case 'faceColor'
            opt.plot.color = [0 0 0; 1 0 0];
    end

    fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d_%s.mat', monkey, experiment, session_id, classifier)), opt);
    format_panel(fh_proj, 'ylim', [-55 55])
    print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('step3_projection_%s_%s_session%d_%s.pdf', monkey, experiment, session_id, classifier)));
end

end



