run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'dpca_Kobak2016')));
monkey = 'Nick';
experiment = 'faceColor';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'dPCA_pair-choice'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'dPCA_pair-choice'); mkdir(InterimDir);

%% dPCA
fh_proj = cell(n_files, 1);
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND')).fnd;
    classifier = 'pair-choice';
    fh_proj{n} = run_dPCA(fnd, InterimDir, FigDir, classifier, monkey, experiment, n);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [5 5], [500 500]*1.5, [-240 240], [-120 120]);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_task_%s_%s.pdf', monkey, experiment)));
print(fh_2, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s.pdf', monkey, experiment)));

%% change axis to plot
fh_proj = cell(n_files, 1);
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    fh_proj{n} = show_dPCA_only(fnd, InterimDir, FigDir, classifier, monkey, experiment, n);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [4 4], [400 400]*1.5);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_task_%s_%s.pdf', monkey, experiment)));
print(fh_2, '-dpdf', fullfile(FigDir, sprintf('dPCA_choice_%s_%s.pdf', monkey, experiment)));

%% [poster] plot example sessions
session_list = [1 6 9];
fh_proj = cell(length(session_list), 1);
for n = 1:length(session_list)
    session_id = session_list(n);
    fnd = load(get_file_path(monkey, experiment, session_id, 'FND_sorted')).fnd;
    classifier = 'pair-choice';
    fh_proj{n} = show_dPCA_only(fnd, InterimDir, FigDir, classifier, monkey, experiment, session_id);
end

[fh_1, fh_2] = plot_in_one_fig_dPCA(fh_proj, [2 3], [300 200]*1.2);
print(fh_1, '-dpdf', fullfile(FigDir, sprintf('dPCA_task_example_session_%s_%s.pdf', monkey, experiment)));

%% [poster] plot pair/choice signal vs. session
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

[task_signal, choice_signal] = deal(cell(2, 3)); % {monkey, experiment}
for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        [~, n_files] = get_file_path(monkey, experiment);
        for session_id = 1:n_files
            data = load(fullfile(InterimDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.mat', monkey, experiment, session_id))).data;
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
            task_signal{monkey_id, exp_id}(session_id) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));

            % calculate choice signal difference in dPC score
            dpc_mn = [mean(dpc(2,:,[1 2]), 3)', mean(dpc(2,:,[3 4]), 3)'];
            I = tstamp>200 & tstamp<600;
            choice_signal{monkey_id, exp_id}(session_id) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));
        end
    end
end

color_list = [44 145 224;
    58 191 153;
    240 169 58]/255;

% plot pair signal version 1 (plot 2 monkeys seperately)
% fh = figure('Position', [50 100 300*1.2 200*1.2]);
% count = 0;
% for monkey_id = 1:length(monkey_list)
%     for exp_id = 1:length(experiment_list)
%         count = count+1;
%         if count==4; continue; end
%         subplot(2,3,count); hold on
%         nses = length(task_signal{monkey_id, exp_id});
%         plot(task_signal{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7)
%         format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Pair signal', ...
%             'xlim', [0.5 nses+0.5], ...
%             'xtick', 1:2:nses)
%         xtickangle(0)
%     end
% end
% axPosition = gca().Position; % 0.7320    0.1307    0.1730    0.3204
% print(fh, '-dpdf', fullfile(FigDir, sprintf('pair_signal.pdf')));

% plot pair signal version 2 (plot 2 monkeys together)
fh = figure('Position', [50 100 300*1.2 200*1.2]);
for exp_id = 1:length(experiment_list)
    subplot(2,3,exp_id); hold on
    nses = max([length(task_signal{1, exp_id}), length(task_signal{2, exp_id})]);
    for monkey_id = 1:length(monkey_list)
        if exp_id==1 && monkey_id==2; continue; end
        
        if monkey_id==1
            plot(task_signal{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7)
        else
            plot(task_signal{monkey_id, exp_id}, '-', 'Color', color_list(exp_id,:), 'LineWidth', 1);
            scatter(1:length(task_signal{monkey_id, exp_id}), task_signal{monkey_id, exp_id}, 7, 'filled', 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', color_list(exp_id,:))
        end
    end
    format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Pair signal', ...
        'xlim', [0.5 nses+0.5], ...
        'xtick', 1:2:nses, ...
        'ylim', [10 60])
    xtickangle(0)
    title(sprintf('Pair 1 vs. %d', exp_id+1))
end
axPosition = gca().Position; % 0.7320    0.1307    0.1730    0.3204
print(fh, '-dpdf', fullfile(FigDir, sprintf('pair_signal.pdf')));

% plot choice signal
fh = figure('Position', [50 100 300*1.2 200*1.2]);
ax = subplot(2,3,6); hold on
for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        nses = length(choice_signal{monkey_id, exp_id});
        if monkey_id==1
            plot(choice_signal{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7);
        else
            plot(choice_signal{monkey_id, exp_id}, '-', 'Color', color_list(exp_id,:), 'LineWidth', 1);
            scatter(1:nses, choice_signal{monkey_id, exp_id}, 7, 'filled', 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', color_list(exp_id,:))
        end
    end
end
format_panel(fh, 'xlabel', '#Session', 'ylabel', 'Choice signal', ...
    'xlim', [0.5 15+0.5], ...
    'xtick', 1:2:15, ...
    'ylim', [20 120])
xtickangle(0)
ax.Position = axPosition;
print(fh, '-dpdf', fullfile(FigDir, sprintf('choice_signal.pdf')));

%% [poster TBU] compare pair signal between early and late sessions
for monkey_id = 1:length(monkey_list)
    monkey = monkey_list{monkey_id};
    figure('Position', [120 300 300 120]);
    for exp_id = 1:length(experiment_list)
        task_signal_ = task_signal{exp_id, monkey_id};
        nsession = length(task_signal_);
        half_ses = round(nsession/2);

        category = [ones(1,half_ses), 2*ones(1,nsession-half_ses)];
        jitter_amount = 0.2;
        jitter = (rand(size(category)) - 0.5) * jitter_amount;
        category = category + jitter;

        subplot(1,3,exp_id); hold on
        scatter(category, task_signal_, 7, color_list(exp_id,:))
        if monkey_id==1; format_panel(gca, 'ylim', [10 60], 'axis', 'normal'); end
        format_panel(gca, 'axis', 'normal', 'xtick', [1 2], 'xticklabel', {'Early', 'Late'}, 'ylabel', 'Pair signal', ...
            'xlim', [1-0.5 2+0.5])
        title(sprintf('Pair 1 vs. %d', exp_id+1))
        xtickangle(0)
        print(gcf, '-dpdf', fullfile(FigDir, sprintf('scatter_pair_signal_%s.pdf', monkey)));
    end
end

%% plot pair/choice signal vs. session in the faceColor task
monkey = 'Nick';
experiment = 'faceColor';

[~, n_files] = get_file_path(monkey, experiment);
[task_signal, choice_signal] = deal(nan(n_files, 1));
for session_id = 1:n_files
    data = load(fullfile(InterimDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.mat', monkey, experiment, session_id))).data;
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
    task_signal(session_id) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));

    % calculate choice signal difference in dPC score
    dpc_mn = [mean(dpc(2,:,[1 2]), 3)', mean(dpc(2,:,[3 4]), 3)'];
    I = tstamp>200 & tstamp<600;
    choice_signal(session_id) = abs(mean(dpc_mn(I,1)) - mean(dpc_mn(I,2)));
end

% plot pair signal
fh = figure('Position', [50 100 250 150]);
plot(task_signal, '.-', 'Color', 'black', 'LineWidth', 0.5, 'MarkerSize', 7)
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Pair signal', ...
    'xlim', [0.5 n_files+0.5], ...
    'xtick', 1:2:n_files)
xtickangle(0)
print(fh, '-dpdf', fullfile(FigDir, sprintf('pair_signal_%s_%s.pdf', monkey, experiment)));

% plot choice signal
fh = figure('Position', [50 100 250 150]);
plot(choice_signal, '.-', 'Color', 'black', 'LineWidth', 0.5, 'MarkerSize', 7)
format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Choice signal', ...
    'xlim', [0.5 n_files+0.5], ...
    'xtick', 1:2:n_files)
xtickangle(0)
print(fh, '-dpdf', fullfile(FigDir, sprintf('choice_signal_%s_%s.pdf', monkey, experiment)));

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
function fh_proj = run_dPCA(fnd, InterimDir, FigDir, classifier, monkey, experiment, session_id)

%% 1. load data/setup parameters
% select unit
r = fnd.FR({2, [100 500]});
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

        ID = ID_dpca;

    case 'pair-choice'
        ID_dpca = task_set; ID_dpca(targ_cho==2) = ID_dpca(targ_cho==2) + max(ID_dpca(:)); % classifier for dPCA
        % trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});

        ID = task_set; ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:)); % classifier for projection
        % trial_classifier_result(ID, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});
end

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
    print(fh, '-dpdf', fullfile(FigDir, sprintf('optimal_lambda_%s_%s_session%d.pdf', monkey, experiment, session_id)));
    save(fullfile(InterimDir, sprintf('optimal_lambda_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data')
end

%% 3. run main dPCA
clear opt;
opt.detrend = false;
opt.tbin = tbin;
opt.t_range = [0 600];
opt.param_combination = param_combination;
opt.param_name = param_name;
opt.show_figure = false;
if regularize; opt.regularization = load(fullfile(InterimDir, sprintf('optimal_lambda_%s_%s_session%d.mat', monkey, experiment, session_id))).data; end

r = fnd.raster_cond(ID_dpca, epoch, 'array'); r = r{1}; % (unit, trial, time, condition)
[fh, data] = dPCA_stim_choice(r, fnd.tstamp{epoch}, opt);
if ~isnan(fh); saveas(fh, fullfile(FigDir, sprintf('dPCA_%s_%s_session%d.pdf', monkey, experiment, session_id)), 'fig'); end
save(fullfile(InterimDir, sprintf('dPCA_result_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data');

%% 4. compute neural data along dPCA axes
clear opt;
opt.epoch = 2;
opt.target = {{'Task', 1}, {'Choice', 1}}; % target is stimulus dimension 1 and choice dimension 1
opt.coefficient_data = load(fullfile(InterimDir, sprintf('dPCA_result_%s_%s_session%d.mat', monkey, experiment, session_id))).data;

data = gen_popresp_dPCA_axis(fnd, ID, opt);
data.cutoff = [find(fnd.tstamp{opt.epoch}==0), find(fnd.tstamp{opt.epoch}==600)];
save(fullfile(InterimDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data');

%% 5. show neural data along dPCA axes
clear opt;
opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
opt.plot = set_plot_opt_2cond('roma', 'roma', max(ID(:))/2);

fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.mat', monkey, experiment, session_id)), opt);
format_panel(fh_proj, 'ylim', [-55 55])
print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.pdf', monkey, experiment, session_id)));

end



function fh_proj = show_dPCA_only(fnd, InterimDir, FigDir, classifier, monkey, experiment, session_id)

%% 1. load data/setup parameters
% select unit
r = fnd.FR({2, [100 500]});
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

        ID = ID_dpca;

    case 'pair-choice'
        ID_dpca = task_set; ID_dpca(targ_cho==2) = ID_dpca(targ_cho==2) + max(ID_dpca(:)); % classifier for dPCA
        % trial_classifier_result(ID_dpca, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});

        ID = task_set; ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:)); % classifier for projection
        % trial_classifier_result(ID, {'coh', 'task_set', 'targ_cho'}, {coh, task_set, targ_cho});
end

%% 4. compute neural data along dPCA axes
clear opt;
opt.epoch = 2;
opt.target = {{'Task', 1}, {'Choice', 1}}; % target is stimulus dimension 1 and choice dimension 1
opt.coefficient_data = load(fullfile(InterimDir, sprintf('dPCA_result_%s_%s_session%d.mat', monkey, experiment, session_id))).data;

data = gen_popresp_dPCA_axis(fnd, ID, opt);
data.cutoff = [find(fnd.tstamp{opt.epoch}==0), find(fnd.tstamp{opt.epoch}==600)];
save(fullfile(InterimDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.mat', monkey, experiment, session_id)), 'data');

%% 5. show neural data along dPCA axes
clear opt;
opt.conv = fspecial('average', [1, 100]); % 100 ms boxcar
opt.plot = set_plot_opt_2cond('roma', 'roma', max(ID(:))/2);
new_color_set = [0 0 0; 44 145 224; 0 0 0; 44 145 224]/255;
opt.plot.color = new_color_set;
opt.less_timepoint = 1;

fh_proj = show_popresp_dPCA_axis(fullfile(InterimDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.mat', monkey, experiment, session_id)), opt);
format_panel(fh_proj, 'ylim', [-55 55])
print(fh_proj, '-dpdf', fullfile(FigDir, sprintf('popresp_dPCA_axis_%s_%s_session%d.pdf', monkey, experiment, session_id)));

end