run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask3'; % learnTask2, learnTask3, learnTask4, faceColor, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'dPCA', 'pair-choice'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'dPCA', 'pair-choice'); mkdir(InterimDir);

%% dPCA
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    step = [1 1 1 0]; % regularization, dPCA, projection, plotting
    run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n);
end

%% plot all sessions
fh_proj = cell(n_files, 1);
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    step = [0 0 0 1]; % regularization, dPCA, projection, plotting
    fh_proj{n} = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, n);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [5 5], [500 500]*1.3, [-40 40], [-100 100]);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_pair_%s_%s.pdf', monkey, experiment)));
print(fh_2, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s.pdf', monkey, experiment)));

%% plot example sessions in Woody learnTask3
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask3'; % learnTask2, learnTask3, learnTask4, faceColor
session_list = [1 6 9];
fh_proj = cell(length(session_list), 1);
for n = 1:length(session_list)
    session_id = session_list(n);
    fnd = load(get_file_path(monkey, experiment, session_id, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    step = [0 0 0 1]; % regularization, dPCA, projection, plotting
    fh_proj{n} = run_dPCA(fnd, InterimDir, FigDir, classifier, step, monkey, experiment, session_id);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [2 3], [300 200]*1.2);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_pair_example_session_%s_%s.pdf', monkey, experiment)));

%% plot magnitude of pair/choice signal vs. session
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

[fh_pair, fh_choice] = plot_signal_magnitude(monkey_list, experiment_list, InterimDir);
print(fh_pair, '-dpdf', fullfile(FigDir, sprintf('pair_signal.pdf')));
print(fh_choice, '-dpdf', fullfile(FigDir, sprintf('choice_signal.pdf')));

%% pair vs. choice 2d plot
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

stat = deal(cell(2, 3)); % {monkey, experiment}
for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        [~, n_files] = get_file_path(monkey, experiment);
        stat{monkey_id, exp_id} = cell(1, n_files);
        for session_id = 1:n_files
            fnd = load(get_file_path(monkey, experiment, session_id, 'FND_sorted')).fnd;
            fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
            targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
            task_set = fnd.getp('task_set'); task_set = task_set(1,:);
            stat{monkey_id, exp_id}{session_id}(:,3) = targ_cor(:);
            stat{monkey_id, exp_id}{session_id}(:,4) = task_set(:);

            data = load(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d.mat', monkey, experiment, session_id))).data;
            dpc = data.dpc_trial; % (dim, time, trial) pair+choice
            tstamp = data.tstamp; % (1, time)

            % smooth dPC score
            for d = 1:size(dpc, 1)
                for c = 1:size(dpc, 3)
                    dpc(d,:,c) = nanconv(dpc(d,:,c), fspecial('average', [1,100]), 'same');
                end
            end

            % calculate pair signal in dPC score
            I = tstamp>200 & tstamp<600;
            dpc_pair_mn = squeeze(mean(dpc(1,I,:),2));
            stat{monkey_id, exp_id}{session_id}(:,2) = dpc_pair_mn(:);

            % calculate choice signal in dPC score
            I = tstamp>200 & tstamp<600;
            dpc_choice_mn = squeeze(mean(dpc(2,I,:),2));
            stat{monkey_id, exp_id}{session_id}(:,1) = dpc_choice_mn(:);
        end
    end
end
save(fullfile(InterimDir, sprintf('pair_vs_choice_2d_plot.mat')), 'stat')

for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        figure('Position', [50 100 600 1000]);
        [~, n_files] = get_file_path(monkey, experiment);
        for session_id = 1:n_files
            subplot(5, 3, session_id); hold on

            I = stat{monkey_id, exp_id}{session_id}(:,3)==1 & stat{monkey_id, exp_id}{session_id}(:,4)==1;
            scatter(stat{monkey_id, exp_id}{session_id}(I,2), stat{monkey_id, exp_id}{session_id}(I,1), ...
                8, 'black')
            I = stat{monkey_id, exp_id}{session_id}(:,3)==1 & stat{monkey_id, exp_id}{session_id}(:,4)==2;
            scatter(stat{monkey_id, exp_id}{session_id}(I,2), stat{monkey_id, exp_id}{session_id}(I,1), ...
                8, 'black', 'filled', 'MarkerFaceAlpha', .5)

            I = stat{monkey_id, exp_id}{session_id}(:,3)==2 & stat{monkey_id, exp_id}{session_id}(:,4)==1;
            scatter(stat{monkey_id, exp_id}{session_id}(I,2), stat{monkey_id, exp_id}{session_id}(I,1), ...
                8, 'red')
            I = stat{monkey_id, exp_id}{session_id}(:,3)==2 & stat{monkey_id, exp_id}{session_id}(:,4)==2;
            scatter(stat{monkey_id, exp_id}{session_id}(I,2), stat{monkey_id, exp_id}{session_id}(I,1), ...
                8, 'red', 'filled', 'MarkerFaceAlpha', .5)
            format_panel(gca, 'xlabel', 'Pair', 'ylabel', 'Choice')
            if session_id==1
                title('black: category 1, empty: pair 1')
            else
                title(sprintf('session %d', session_id))
            end
        end
        print(gcf, '-dpdf', fullfile(FigDir, sprintf('pair_vs_choice_2d_plot_category_%s_%s.pdf', monkey, experiment)));
    end
end

%% plot pair signal of switch vs. non-switch trials (??)
% get pair signal, cond switch, session (slow)
session_list = 1:n_files;
[pair_signal, n_after_switch, session] = deal([]);
for i = 1:length(session_list)
    n = session_list(i);
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
    n_after_switch_ = fnd.getp('n_after_switch'); n_after_switch_ = n_after_switch_(1,:);
    task_set_ = fnd.getp('task_set'); task_set_ = task_set_(1,:);
    targ_cho_ = fnd.getp('targ_cho'); targ_cho_ = targ_cho_(1,:);

    data = load(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d.mat', monkey, experiment, n))).data;
    dpc_trial = data.dpc_trial;
    pair_signal_ = squeeze(mean(dpc_trial(1,:,:), 2));

    pair_signal = [pair_signal; pair_signal_(:)];
    n_after_switch = [n_after_switch; n_after_switch_(:)];
    session = [session; i*ones(size(pair_signal_(:)))];
end
save(fullfile(InterimDir, sprintf('pair_signal_vs_switch_%s_%s.mat', monkey, experiment)), 'pair_signal', 'n_after_switch', 'session')

% stat
load(fullfile(InterimDir, sprintf('pair_signal_vs_switch_%s_%s.mat', monkey, experiment)));

stat = cell(1, n_files);
for i = 1:length(session_list)
    I = ~isnan(n_after_switch) & session==i;
    [pair_signal_mn, pair_signal_se] = calcGroupMean(pair_signal(I), n_after_switch(I), unique(n_after_switch(I)));

    idx = 30;
    pair_signal_compare = [mean(abs(pair_signal(I&n_after_switch>idx))); mean(abs(pair_signal(I&n_after_switch<=idx)))]; % non-switch vs. switch

    p = run_ttest2(abs(pair_signal(I&n_after_switch>idx))', abs(pair_signal(I&n_after_switch<=idx))', '=');

    stat{i}.pair_signal_mn = pair_signal_mn;
    stat{i}.pair_signal_se = pair_signal_se;
    stat{i}.u_n_after_switch = unique(n_after_switch(I));
    stat{i}.pair_signal_compare = pair_signal_compare;
    stat{i}.p = p;
end

% plot
pair_signal_compare = cellfun(@(x) x.pair_signal_compare, stat, 'uni', 0);
pair_signal_compare = cell2mat(pair_signal_compare);
p = cellfun(@(x) x.p, stat);

figure('Position', [50 100 300 150]);
b = bar(pair_signal_compare');
b(1).FaceColor = [0 0 0];
b(2).FaceColor = [1 0 0];
for i = 1:length(p)
    if p(i)<0.05
        text(i, max(pair_signal_compare(:))*1.1, '*', 'HorizontalAlignment', 'center', 'FontSize', 6)
    else
        text(i, max(pair_signal_compare(:))*1.1, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 6)
    end
end
format_panel(gcf, 'axis', 'normal', ...
    'xlabel', '#Session', 'ylabel', 'Pair signal', 'ylim', [0 max(pair_signal_compare(:))*1.2])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('pair_signal_vs_switch_%s_%s.pdf', monkey, experiment)));

%% stat (TBU)
fprintf('\nSpearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', task_signal, 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Pair signal did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Pair signal increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Pair signal decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

% stat
fprintf('\nSpearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', choice_signal, 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Choice signal did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Choice signal increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Choice signal decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

%% [stat] orthogonal pair and choice axis
[orth_scores, orth_scores_perm_mn, p_values] = deal(nan(n_files, 1));
for n = 1:n_files
    data = load(fullfile(InterimDir, sprintf('popresp_dPCA_axis_session%d.mat', n))).data;
    coef = load(fullfile(InterimDir, sprintf('dPCA_result_session%d.mat', n))).data;
    tdim = data.tdim; % [task axis, choice axis]
    W = coef.W; % (unit, dim)

    task_axis = W(:,tdim(1)); % (unit, 1)
    choice_axis = W(:,tdim(2)); % (unit, 1)

    % get orthogonality score: abs(v1*v2 / (norm(v1)*norm(v2))), range is [0 1]
    cos_sim = task_axis' * choice_axis / (norm(task_axis) * norm(choice_axis));
    orth_scores(n) = abs(cos_sim);
    
    % statistical test for orthogonality using permutation test
    n_permutations = 10000;
    orth_scores_perm = zeros(n_permutations, 1);
    for perm_idx = 1:n_permutations
        % randomly shuffle the choice axis while keeping task axis fixed
        perm_choice_axis = choice_axis(randperm(length(choice_axis)));
        perm_cos_sim = (task_axis' * perm_choice_axis) / (norm(task_axis) * norm(perm_choice_axis));
        orth_scores_perm(perm_idx) = abs(perm_cos_sim);
    end
    orth_scores_perm_mn(n) = mean(orth_scores_perm);
    p_values(n) = mean(orth_scores_perm >= orth_scores(n));
end

% plot orthogonality scores across sessions
fh = figure('Position', [100 100 300 300]); hold on
plot(1:n_files, orth_scores, '.-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 11)
plot(1:n_files, orth_scores_perm_mn, '.-', 'Color', 'r', 'LineWidth', 1, 'MarkerSize', 11)
legend({'Data', 'Permutated data'}); legend boxoff
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Similarity score (|cos θ|)')
title('Pair-choice axis similarity')
print(fh, '-dpdf', fullfile(FigDir, sprintf('%s_%s_similarity_scores.pdf', experiment, monkey)));

% plot p-value across sessions
fh = figure('Position', [100 100 300 300]); hold on
plot(1:n_files, p_values, '.-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 11)
plot(xlim, [0.05 0.05], '--', 'Color', 'r', 'LineWidth', 1)
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'P-value')
title('Pair-choice axis similarity')
print(fh, '-dpdf', fullfile(FigDir, sprintf('%s_%s_similarity_scores_pvalue.pdf', experiment, monkey)));

% statistical summary
fprintf('\nTest pair-choice axis similarity using permutation test:\n')
fprintf('Mean similarity score (|cos θ|): %.4f ± %.4f (SD)\n', mean(orth_scores), std(orth_scores))
fprintf('Range: [%.4f, %.4f]\n', min(orth_scores), max(orth_scores))

% count sessions with significant similarity
fprintf('Sessions with non-significant similarity (i.e., orthogonal) (p>0.05): %d/%d\n', sum(p_values>0.05), n_files)








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

    data = gen_popresp_dPCA_axis(fnd, ID, opt);
    data.cutoff = [find(fnd.tstamp{opt.epoch}==0), find(fnd.tstamp{opt.epoch}==600)];
    save(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data');
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
    opt.less_timepoint = true;

    fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('step3_projection_%s_%s_session%d.mat', monkey, experiment, session_id)), opt);
    format_panel(fh_proj, 'ylim', [-55 55])
    print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('step3_projection_%s_%s_session%d.pdf', monkey, experiment, session_id)));
end

end


