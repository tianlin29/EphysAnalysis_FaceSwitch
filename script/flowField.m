run('../Initialize.m');
monkey = 'Woody';
experiment = 'learnTask3';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'flowField'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'flowField'); mkdir(InterimDir);

%% check raw PCA
TASK_ID = 2;
VERBOSE = false;
PREPROCESS_FND = true;

figure('Position', [500 100 850 900]);
for n = 1:n_files
    n
    % load file
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    if PREPROCESS_FND
        % remove low FR units
        FR = fnd.FR({1, [100 500]}); % (unit, trial)
        I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);
        fprintf('%d low FR units are removed.\n', sum(~I))

        % remove units changed over time
        % I = ~isnan(fnd.misc.SNR);
        % fnd = fnd.set_unit_criteria('custom', I);
        % fprintf('%d units changed over time.\n', sum(~I))

        % select trial
        fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
        fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID);

        if TASK_ID==2
            if strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[1 2 3 4 5])
                fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=60);
            elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && n==6
                fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=18);
            elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[7 8 9 10])
                fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=12);
            elseif strcmp(monkey, 'Nick') && strcmp(experiment, 'learnTask2') && any(n==[11 12 13 14 15])
                fnd = fnd.extract_trial(abs(fnd.getp('morph_level'))*100>=6);
            end
        end
    end

    % get ID and PSTH
    task_set = fnd.getp('task_set');
    morph_level = fnd.getp('morph_level')*100;
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');

    trc = trial_classifier('stim_group', {[0 18], [18 32], [32 60], [60 80], [80 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    ID(task_set~=TASK_ID) = NaN; % one task only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % get dPC axes
    psth = fnd.PSTH(ID, {'gaussian', 20}); psth = psth{1}; % (unit, time, cond)
    [nunit, ntime, ntrial] = size(fnd.data{1});
    ncond = max(ID(:));

    target = {{'Choice', 1}, {'Choice', 2}}; % target is stimulus dimension 1 and choice dimension 1
    coef = load(fullfile(MainInterimDir, 'dPCA_pair-choice', sprintf('dPCA_result_%s_%s_session%d.mat', monkey, experiment, n))).data;

    psth = bsxfun(@rdivide, bsxfun(@minus, psth, coef.MU), coef.SIGMA);

    ntarg = length(target);
    tdim = nan(ntarg, 1); % [choice axis 1, choice axis 2, task axis 1]
    for t = 1:ntarg
        I = find(strcmp(target{t}{1}, coef.margNames));
        if isempty(I)
            error('No such dimension: %s', target{t}{1});
        end
        tmpdim = find(coef.whichMarg == I, target{t}{2});
        tdim(t) = tmpdim(end);
    end

    W = coef.W; % (unit, dim)
    coef_decision_2D = W(:,tdim([1 2])); % (unit, 2pc)

    % project to dPC
    score = coef_decision_2D' * psth(:,:);
    score = reshape(score, 2, ntime, []);

    opt.plot = set_plot_opt('vik', ncond);
    subplot(4,4,n); hold on
    for c = 1:ncond
        plot(score(1,:,c), score(2,:,c), 'Color', opt.plot.color(c,:));
    end
    format_panel(gcf, 'xlabel', 'dPC1', 'ylabel', 'dPC2');
    title(sprintf('session %d', n))

    % get PC axes
    % psth = fnd.PSTH(ID, {'gaussian', 20}, [], true); psth = psth{1}; % (unit, time, cond)
    % [nunit, ntime, ntrial] = size(fnd.data{1});
    % ncond = max(ID(:));
    % 
    % [coeff, score, latent] = pca(psth(:,:)'); % psth(:,:) ..(observation, unit)
    % score = reshape(score', nunit, ntime, ncond);
    % 
    % opt.plot = set_plot_opt('vik', ncond);
    % subplot(4,4,n); hold on
    % for c = 1:ncond
    %     plot(score(1,:,c), score(2,:,c), 'Color', opt.plot.color(c,:));
    % end
    % format_panel(gcf, 'xlabel', 'PC1', 'ylabel', 'PC2');
    % title(sprintf('session %d', n))
end


%% flow field (task 1 or 2, dPCA) 我觉得噪声更大了
TASK_ID = 1;
VERBOSE = false;
PREPROCESS_FND = true;

fh_idv = cell(n_files, 1);
for n = 1:n_files
    % load file
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    if PREPROCESS_FND
        % remove low FR units
        FR = fnd.FR({1, [100 500]}); % (unit, trial)
        I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);
        fprintf('%d low FR units are removed.\n', sum(~I))

        % remove units changed over time
        % I = ~isnan(fnd.misc.SNR);
        % fnd = fnd.set_unit_criteria('custom', I);
        % fprintf('%d units changed over time.\n', sum(~I))

        % select trial
        fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
        fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID);
    end

    % get dPC axis to plot
    coef = load(fullfile(MainInterimDir, 'dPCA_pair-choice', sprintf('dPCA_result_%s_%s_session%d.mat', monkey, experiment, n))).data;
    target = {{'Choice', 1}, {'Choice', 2}}; % actually coh 1 and choice 1

    ntarg = length(target);
    tdim = nan(ntarg, 1); % [choice axis 1, choice axis 2, task axis 1]
    for t = 1:ntarg
        I = find(strcmp(target{t}{1}, coef.margNames));
        if isempty(I)
            error('No such dimension: %s', target{t}{1});
        end
        tmpdim = find(coef.whichMarg == I, target{t}{2});
        tdim(t) = tmpdim(end);
    end

    W = coef.W; % (unit, dim)
    coef_decision_2D = W(:,tdim([1 2])); % (unit, 2pc)
    % coef_task_1D = W(:,tdim(3)); % (unit, 1pc)

    % detPSTH = fnd.PSTH([], {'boxcar', 100}, [], [], 2); detPSTH = detPSTH{2};

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    ID = targ_cho; ID(targ_cho~=targ_cor) = NaN; % correct trials only
    ID(task_set~=TASK_ID) = NaN; % one task only
    if VERBOSE; trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set}); end

    % calculate flow field
    clear opt;
    opt.PC_coef = coef_decision_2D;
    opt.detPSTH = coef.MU;
    opt.dim = [1 2];
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
    opt.bin_num = 20;

    data = calc_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, opt);

    % show flow field (PC 1-2)
    clear pltopt;
    pltopt.min_sample = size(fnd.data{1}, 3);
    pltopt.trj_color = [
        0.00  0.45  0.74;
        0.85  0.33  0.10];
    pltopt.legend = {'Choice 1', 'Choice 2'};
    pltopt.xlabel = 'dPC 1';
    pltopt.ylabel = 'dPC 2';
    pltopt.title = sprintf('Task %d', TASK_ID);

    fh_idv{n} = show_flow_field(data, pltopt);
    print(fh_idv{n}, '-dpdf', fullfile(FigDir, sprintf('flow_field_dPC%d-%d_pair%d_%s_%s_session%d.pdf', opt.dim(1), opt.dim(2), TASK_ID, monkey, experiment, n)));
end

fh_all = plot_in_one_fig_flow_field(fh_idv, [5 3], [300 500*1.2]*2);
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('flow_field_dPC%d-%d_pair%d_%s_%s_summary.pdf', opt.dim(1), opt.dim(2), TASK_ID, monkey, experiment)));

%% flow field (task 1 or 2, PCA)
TASK_ID = 1;
VERBOSE = false;
PREPROCESS_FND = true;

fh_idv = cell(n_files, 1);
for n = 1:n_files
    % load file
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    if PREPROCESS_FND
        % remove low FR units
        FR = fnd.FR({1, [100 500]}); % (unit, trial)
        I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);
        fprintf('%d low FR units are removed.\n', sum(~I))

        % remove units changed over time
        % I = ~isnan(fnd.misc.SNR);
        % fnd = fnd.set_unit_criteria('custom', I);
        % fprintf('%d units changed over time.\n', sum(~I))

        % select trial
        fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
        fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_ID);
    end

    % get PC space
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    ID = targ_cho; ID(targ_cho~=targ_cor) = NaN; % correct trials only
    ID(task_set~=TASK_ID) = NaN; % one task only
    if VERBOSE; trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set}); end

    data = fnd.PSTH(ID, {'gaussian', 20}); % the kernel is actually {'gaussian', [1 20*6], 20}
    
    clear opt;
    opt.epoch = 1;
    opt.PC_range = [250 600];
    [coef, detPSTH] = getPC_coef(fnd.tstamp, data, opt); % coeff ..(feature, pc); detPSTH ..(unit, time)

    % calculate flow field
    clear opt;
    opt.PC_coef = coef;
    opt.detPSTH = detPSTH;
    opt.dim = [1 2];
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
    opt.bin_num = 20;

    data = calc_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, opt);

    % show flow field (PC 1-2)
    clear pltopt;
    pltopt.min_sample = size(fnd.data{1}, 3);
    pltopt.trj_color = [
        0.00  0.45  0.74;
        0.85  0.33  0.10];
    pltopt.legend = {'Choice 1', 'Choice 2'};
    pltopt.xlabel = 'PC 1';
    pltopt.ylabel = 'PC 2';
    pltopt.title = sprintf('Pair %d', TASK_ID);
    
    fh_idv{n} = show_flow_field(data, pltopt);
    xlim([-100 100]); ylim([-80 80])
    print(fh_idv{n}, '-dpdf', fullfile(FigDir, sprintf('flow_field_PC%d-%d_pair%d_%s_%s_session%d.pdf', opt.dim(1), opt.dim(2), TASK_ID, monkey, experiment, n)));
end

fh_all = plot_in_one_fig_flow_field(fh_idv, [5 3], [300 500*1.2]*2);
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('flow_field_PC%d-%d_pair%d_%s_%s_summary.pdf', opt.dim(1), opt.dim(2), TASK_ID, monkey, experiment)));

%% flow field (input, task 1 or 2, PCA) 需要严格的格外检查
TASK_ID = 1;
VERBOSE = false;
PREPROCESS_FND = true;

for n = 1:n_files
    % load file
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    if PREPROCESS_FND
        % remove low FR units
        FR = fnd.FR({1, [100 500]}); % (unit, trial)
        I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);
        fprintf('%d low FR units are removed.\n', sum(~I))

        % remove units changed over time
        % I = ~isnan(fnd.misc.SNR);
        % fnd = fnd.set_unit_criteria('custom', I);
        % fprintf('%d units changed over time.\n', sum(~I))

        % select trial
        fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
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
    opt.PC_range = [0 700];
    opt.conv_kernel = fspecial('average', [1 100]);
    opt.epoch = 1;

    data = fnd.PSTH(ID, {'gaussian', 20});
    [coef, detPSTH] = getPC_coef(fnd.tstamp, data, opt);

    % define external input
    clear opt;
    opt.latency = 100; % ms
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);

    sdur = fnd.getp('stim_dur') * 1e3;
    morph = fnd.getp('morph_level');
    % task = fnd.getp('task_set');
    % task(task(:) == 2) = -1; % 1 vs -1
    task = zeros(size(morph));
    input = gen_input_matrix(fnd.tstamp{1}, sdur(1,:), {morph(1,:), task(1,:)}, opt); % don't use task as input

    % calculate flow field
    clear opt;
    opt.PC_coef = coef;
    opt.detPSTH = detPSTH;
    opt.dim = [1 2];
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
    opt.bin_num = 20;
    data = calc_input_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, input, opt);

    % show flow field
    clear pltopt;
    pltopt.min_sample = 1000;
    pltopt.trj_color = [
        0.00  0.45  0.74;
        0.85  0.33  0.10];
    pltopt.legend = {'Choice 1', 'Choice 2'};
    pltopt.input_col = [0 0 0;.8 .8 .8];
    pltopt.xlabel = 'PC 1';
    pltopt.ylabel = 'PC 2';
    fh = show_flow_field(data, pltopt);
    title('K arrow: morph');
    print(fh, '-dpdf', fullfile(FigDir, sprintf('input_flow_field_PC%d-%d_pair%d_%s_%s_session%d.pdf', opt.dim(1), opt.dim(2), TASK_ID, monkey, experiment, n)));
end

%% Fano rate (mean FR / std)
fano_mn = nan(4, n_files);
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

    % select unit
    r = fnd.FR({2, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));

    % get Fano factor
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    ID = targ_cho; ID(targ_cho~=targ_cor) = NaN; % correct trials only
    ID(task_set==2) = ID(task_set==2) + max(ID(:));
    trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set});

    r = fnd.raster_cond(ID, 2, 'cell'); r = r{1}; % when 'format' is 'cell': {1, condition}, (unit, trial, time)

    I = fnd.tstamp{2}>100 & fnd.tstamp{2}<500;
    r = cellfun(@(x) mean(x(:,:,I),3), r, 'uni', 0);

    r_mn = cellfun(@(x) mean(x,2), r, 'uni', 0);
    r_std = cellfun(@(x) std(x,[],2), r, 'uni', 0);
    fano = cellfun(@(x) mean(x,2) ./ std(x,[],2), r, 'uni', 0);

    fano_mn(:,n) = cellfun(@(x) mean(x), fano)';
end
save(fullfile(InterimDir, 'fano_mn.mat'), 'fano_mn');

% plot
figure;
subplot(1,2,1); hold on
plot(fano_mn(1,:), '.-', 'MarkerSize', 7); % chocie 1
plot(fano_mn(2,:), '.-', 'MarkerSize', 7); % chocie 2
legend({'Chocie 1', 'Choice 2'}); title('Task 1')

subplot(1,2,2); hold on
plot(fano_mn(3,:), '.-', 'MarkerSize', 7); % chocie 1
plot(fano_mn(4,:), '.-', 'MarkerSize', 7); % chocie 2
legend({'Chocie 1', 'Choice 2'}); title('New task')

format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Fano factor')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('%s_%s_fano_factor.pdf', experiment, monkey)));

%% [stat] Fano factor 无显著趋势
fano_mn = load(fullfile(InterimDir, 'fano_mn.mat')).fano_mn;
fano_mn_task1 = mean(fano_mn([1 2],:), 1);
fano_mn_task2 = mean(fano_mn([3 4],:), 1);

fprintf('\n%s %s:\n', monkey, experiment)
fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', fano_mn_task1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Fano factor of pair 1       did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Fano factor of pair 1       increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Fano factor of pair 1       decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end
[rho, p_val] = corr((1:n_files)', fano_mn_task2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Fano factor of the new pair did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Fano factor of the new pair increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Fano factor of the new pair decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

%% raw flow field
for n = 1:n_files

% get axes
target = {{'Choice', 1}, {'Choice', 2}}; % target is stimulus dimension 1 and choice dimension 1
coef = load(fullfile(MainInterimDir, 'dPCA_pair-choice', sprintf('dPCA_result_%s_%s_session%d.mat', monkey, experiment, n))).data;

ntarg = length(target);
tdim = nan(ntarg, 1); % [choice axis 1, choice axis 2, task axis 1]
for t = 1:ntarg
    I = find(strcmp(target{t}{1}, coef.margNames));
    if isempty(I)
        error('No such dimension: %s', target{t}{1});
    end
    tmpdim = find(coef.whichMarg == I, target{t}{2});
    tdim(t) = tmpdim(end);
end

W = coef.W; % (unit, dim)
coef_decision_2D = W(:,tdim([2 1])); % (unit, 2pc)
% coef_task_1D = W(:,tdim(3)); % (unit, 1pc)

% load fnd
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;

% select unit
r = fnd.FR({2, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% select trial
fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));
[nunit, ntime, ntrial] = size(fnd.data{2});

% project to pair axis
FR = fnd.FR({2, [100 500]}); % (unit, trial)
FR_proj_pair = coef_task_1D' * FR; % (1, trial)

task_set = fnd.getp('task_set'); task_set = task_set(1,:);
FR_proj_pair_1 = mean(FR_proj_pair(task_set==1));
FR_proj_pair_2 = mean(FR_proj_pair(task_set==2));
if FR_proj_pair_1>=FR_proj_pair_2
    fprintf('Flip pair axis. Pair 1 should have more negative pair signal.\n')
    coef_task_1D = -coef_task_1D;
    FR_proj_pair = coef_task_1D' * FR; % (1, trial)
end

nbin = 8;
pair_signal_list = prctile(FR_proj_pair, linspace(0,100,nbin+1));

% reshape
figure('Position', [130 500 1400 110]);
for b = 1:nbin
    trial_idx = (FR_proj_pair>pair_signal_list(b)) & (FR_proj_pair<=pair_signal_list(b+1));

    % select trial
    fnd_new = fnd.extract_trial(repmat(trial_idx, [nunit, 1]));

    % PSTH
    task_set = fnd_new.getp('task_set');
    morph_level = fnd_new.getp('morph_level')*100;
    targ_cho = fnd_new.getp('targ_cho');
    targ_cor = fnd_new.getp('targ_cor');

    % coh ID
    trc = trial_classifier('stim_group', {[0 12], [12 24], [24 40], [40 60], [60 Inf]}, 'plus_cat', 1, 'include_0coh', false);
    [ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
    trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho', 'task_set'}, {morph_level, targ_cor, targ_cho, task_set});

    % DIY ID
    % [~, ~, ID] = unique(morph_level(1,:));
    % ntr = histcounts(ID, 1:max(ID)+1);
    % for c = 1:max(ID)
    %     if ntr(c)<3
    %         ID(ID==c) = NaN;
    %     end
    % end
    % ID = repmat(ID', [nunit,1]);
    % trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {morph_level, targ_cor, targ_cho});

    % project to decision space
    % psth = fnd_new.PSTH(ID, {'boxcar', 100}, [], [], 2); psth = psth{2};
    % psth_ = psth(:,:);
    % psth_proj_dec = coef_decision_2D' * psth_;
    % psth_proj_dec = reshape(psth_proj_dec, 2, ntime, []);
    % 
    % opt.plot = set_plot_opt('vik', max(ID(:)));
    % subplot(1,nbin,b); hold on
    % for c = 1:max(ID(:))
    %     plot(psth_proj_dec(1,:,c), psth_proj_dec(2,:,c), 'Color', opt.plot.color(c,:))
    % end

    % PSTH project to PC
    psth = fnd_new.PSTH(ID, {'boxcar', 100}, [], true, 2); psth = psth{2};
    [coeff, score, latent] = pca(psth(:,:)'); % psth(:,:) ..(observation, unit); coeff ..(feature, pc)
    score = reshape(score', nunit, ntime, []);

    opt.plot = set_plot_opt('vik', max(ID(:)));
    subplot(1,nbin,b); hold on
    for c = 1:max(ID(:))
        plot(score(1,:,c), score(2,:,c), 'Color', opt.plot.color(c,:))
    end
    if b==1
        title('Pair 1')
    elseif b==nbin
        title('New Pair')
    end

    % % raster project to PC
    % r = fnd_new.raster(2); r = r{1}*1000; % (unit, time, trial)
    % r_mn = mean(r,3);
    % r_detrend = r - r_mn;
    % 
    % score = coeff * r(:,:);
    % score = reshape(score', nunit, ntime, []);
    % 
    % for c = 1:size(score,3)
    %     plot(score(1,:,c), score(2,:,c))
    % end
    % if b==1
    %     title('Pair 1')
    % elseif b==nbin
    %     title('New Pair')
    % end
end
format_panel(gcf, 'xlim', [-60 60], 'ylim', [-60 60], 'axis', 'normal')

end


