run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Nick';
experiment = 'learnTask2';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

task_color = [99 97 172; 242 128 128]/255;

%% dPCA
fh_proj = cell(n_files, 1);
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fh_proj{n} = run_dPCA(fnd, InterimDir, FigDir, n);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [4 4], [400 400]*1.5);
print(fh_1, '-dpdf', fullfile(FigDir, 'dPCA_task.pdf'));
print(fh_2, '-dpdf', fullfile(FigDir, 'dPCA_choice.pdf'));

%% [stat] dPCA pair/choice signal
[task_signal, choice_signal] = deal(nan(n_files, 1));
for n = 1:n_files
    data = load(fullfile(InterimDir, sprintf('popresp_dPCA_axis_session%d.mat', n))).data;
    dpc = data.dpc; % (dim, time, cond) 2task*2choice
    tstamp = data.tstamp; % (1, time)

    % smooth dPC score
    for d = 1:size(dpc, 1)
        for c = 1:size(dpc, 3)
            dpc(d,:,c) = nanconv(dpc(d,:,c), fspecial('average', [1,100]), 'same');
        end
    end
    
    % calculate task signal difference in dPC score
    dpc_mn = [mean(dpc(1,:,[1 3]), 3)', mean(dpc(1,:,[2 4]), 3)'];
    I = tstamp>200 & tstamp<600;
    task_signal(n) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));

    % calculate choice signal difference in dPC score
    dpc_mn = [mean(dpc(2,:,[1 2]), 3)', mean(dpc(2,:,[3 4]), 3)'];
    I = tstamp>200 & tstamp<600;
    choice_signal(n) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));
end

% plot pair signal
fh = figure('Position', [100 100 300 300]); hold on
plot(1:n_files, task_signal, '.-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 11)
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Pair signal')
print(fh, '-dpdf', fullfile(FigDir, 'pair_signal.pdf'));

% stat
fprintf('\nSpearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', task_signal, 'Type', 'Spearman');
fprintf('Pair signal decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))

% plot choice signal
fh = figure('Position', [100 100 300 300]); hold on
plot(1:n_files, choice_signal, '.-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 11)
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Choice signal')
print(fh, '-dpdf', fullfile(FigDir, 'choice_signal.pdf'));

% stat
fprintf('\nSpearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', choice_signal, 'Type', 'Spearman');
fprintf('Choice signal did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))

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
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Orthogonality score (|cos θ|)')
title('Pair-choice axis orthogonality')
print(fh, '-dpdf', fullfile(FigDir, 'orthogonality_scores.pdf'));

% plot p-value across sessions
fh = figure('Position', [100 100 300 300]); hold on
plot(1:n_files, p_values, '.-', 'Color', 'k', 'LineWidth', 1, 'MarkerSize', 11)
plot(xlim, [0.05 0.05], '--', 'Color', 'r', 'LineWidth', 1)
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'P-value')
title('Pair-choice axis orthogonality')
print(fh, '-dpdf', fullfile(FigDir, 'orthogonality_scores_pvalue.pdf'));

% statistical summary
fprintf('\nTest pair-choice axis orthogonality using permutation test:\n')
fprintf('Mean orthogonality score (|cos θ|): %.4f ± %.4f (SD)\n', mean(orth_scores), std(orth_scores))
fprintf('Range: [%.4f, %.4f]\n', min(orth_scores), max(orth_scores))

% count sessions with significant orthogonality (p < 0.05)
fprintf('Sessions with significant orthogonality (p < 0.05): %d/%d\n', sum(p_values < 0.05), n_files)



%% functions
function fh_proj = run_dPCA(fnd, InterimDir, FigDir, session_id)

%% 1. load data/setup parameters
% select unit
r = fnd.FR({2, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% classifier
coh = fnd.getp('morph_level')*100;
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');
task_set = fnd.getp('task_set');

ID_dpca = task_set; ID_dpca(targ_cho==2) = ID_dpca(targ_cho==2) + max(ID_dpca(:)); % classifier for dPCA
% trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});

ID = task_set; ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:)); % classifier for projection
% trial_classifier_result(ID, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});

% options
param_combination = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
param_name = {'Task', 'Choice', 'Condition-independent', 'S/C Interaction'};
tbin = 20; % ms
regularize = true;
epoch = 2;

%% 2. determine regularization parameter for dPCA
if regularize
    clear opt;
    opt.tbin = tbin; % ms
    opt.t_range = [0 600];
    opt.param_combination = param_combination;
    opt.param_name = param_name;

    r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
    [fh, data] = dPCA_optimize_regularization(r, fnd.tstamp{epoch}, opt);
    print(fh, '-dpdf', fullfile(FigDir, sprintf('optimal_lambda_session%d.pdf', session_id)));
    save(fullfile(InterimDir, sprintf('optimal_lambda_session%d.mat', session_id)), 'data');
end

%% 3. run main dPCA
clear opt;
opt.detrend = true;
opt.tbin = tbin;
opt.t_range = [200 600];
opt.param_combination = param_combination;
opt.param_name = param_name;
opt.show_figure = false;
if regularize; opt.regularization = load(fullfile(InterimDir, sprintf('optimal_lambda_session%d.mat', session_id))).data; end

r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
[fh, data] = dPCA_stim_choice(r, fnd.tstamp{epoch}, opt);
if ~isnan(fh); saveas(fh, fullfile(FigDir, sprintf('dPCA_session%d.fig', session_id)), 'fig'); end
save(fullfile(InterimDir, sprintf('dPCA_result_session%d.mat', session_id)), 'data');

%% 4. compute neural data along dPCA axes
clear opt;
opt.epoch = 2;
opt.target = {{'Task', 1}, {'Choice', 1}}; % target is stimulus dimension 1 and choice dimension 1
opt.coefficient_data = load(fullfile(InterimDir, sprintf('dPCA_result_session%d.mat', session_id))).data;

data = gen_popresp_dPCA_axis(fnd, ID, opt);
save(fullfile(InterimDir, sprintf('popresp_dPCA_axis_session%d.mat', session_id)), 'data');

%% 5. show neural data along dPCA axes
clear opt;
opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
opt.plot = set_plot_opt_2cond('roma', 'roma', 2);

fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('popresp_dPCA_axis_session%d.mat', session_id)), opt);
format_panel(fh_proj, 'ylim', [-55 55])
print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('popresp_dPCA_axis_session%d.pdf', session_id)));

end


