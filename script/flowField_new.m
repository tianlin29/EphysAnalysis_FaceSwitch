run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, faceColor_passiveLong, threeExemplar
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'flowField'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'flowField'); mkdir(InterimDir);

%% raw flow field
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);
fnd = fnd.extract_trial(fnd.getp('task_set')==2);

% get smooth detrended raster
r = fnd.raster(1); r = r{1}*1000; % (unit, time, trial)
[nunit, ntime, ntrial] = size(r);
r_mn = mean(r, 3);
r_detrend = r - r_mn;

kernel = fspecial('average', [1, 100]);
r_detrend = nanconv(r_detrend, kernel, 'same');

% PCA
targ_cho = fnd.getp('targ_cho');
ID = targ_cho;
ncond = max(ID(:));
psth = nan(nunit, ntime, ncond);
for c = 1:ncond
    psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
end

[coeff, score, latent] = pca(psth(:,:)');
npc = 2;
r_pc = reshape(coeff(:,1:npc)'*r_detrend(:,:), npc, ntime, ntrial);
psth_pc = reshape(score(:,1:npc)', npc, ntime, ncond);

% creat mesh
ngrid = 20;
x_edge = linspace(min(r_pc(:)), max(r_pc(:)), ngrid);
y_edge = linspace(min(r_pc(:)), max(r_pc(:)), ngrid);
[X, Y] = meshgrid(x_edge, y_edge);

% calculate velocity
[U_data, V_data] = deal(cell(ngrid, ngrid));
for i = 1:ngrid
    for j = 1:ngrid
        U_data{i,j} = [];
        V_data{i,j} = [];
    end
end

for t = 1:ntime-1
    t
    for tr = 1:ntrial
        x_t0 = r_pc(1,t,tr); x_t1 = r_pc(1,t+1,tr);
        y_t0 = r_pc(2,t,tr); y_t1 = r_pc(2,t+1,tr);
        col_idx = find(x_t0 >= x_edge(1:end-1) & x_t0 <= x_edge(2:end), 1);
        row_idx = find(y_t0 >= y_edge(1:end-1) & y_t0 <= y_edge(2:end), 1);

        u_t0 = x_t1 - x_t0;
        v_t0 = y_t1 - y_t0;
        U_data{row_idx, col_idx} = cat(1, U_data{row_idx, col_idx}, u_t0);
        V_data{row_idx, col_idx} = cat(1, V_data{row_idx, col_idx}, v_t0);
    end
end

[U, V] = deal(nan(ngrid, ngrid));
for i = 1:ngrid
    for j = 1:ngrid
        nsample = length(U_data{i,j});
        if nsample>=1000
            U(i,j) = mean(U_data{i,j});
            V(i,j) = mean(V_data{i,j});
        end
    end
end

% plot
opt.plot = set_plot_opt('vik', ncond);
figure; hold on
quiver(X, Y, U, V)
for c = 1:ncond
    plot(psth_pc(1,:,c), psth_pc(2,:,c), 'Color', opt.plot.color(c,:));
end


%% ###### Load data (necessary for all steps below) #######
TASK_SET_ID = 2;
[u, v, x, y] = deal(cell(n_files, 1));
for n = 1:n_files
    n
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % fnd = fnd.extract_trial(fnd.getp('task_set')==TASK_SET_ID);

    % get PC space to plot
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor'));

    ID = targ_cho;
    ID(correct~=1) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cho', 'correct'}, {task_set, targ_cho, correct})

    psth = fnd.PSTH(ID, {'gaussian', 20});

    clear opt;
    opt.PC_range = [250 600];
    opt.epoch = 1;
    [coef, detPSTH] = getPC_coef(fnd.tstamp, psth, opt);

    % get dPC axis
    data = load(fullfile(MainInterimDir, 'dPCA\cat-choice\', sprintf('step2_dPCA_%s_%s_session%d.mat', monkey, experiment, n))).data;

    target = {{'Choice', 1}, {'Category', 1}};
    ntarg = length(target);
    tdim = nan(ntarg,1);
    for ta = 1:ntarg
        I = find(strcmp(target{ta}{1}, data.margNames));
        if isempty(I)
            error('No such dimension: %s', target{ta}{1});
        end
        tmpdim = find(data.whichMarg == I, target{ta}{2});
        tdim(ta) = tmpdim(end);
    end

    coef = data.W(:, tdim);

    % calculate flow field
    clear opt;
    opt.PC_coef = coef;
    opt.detPSTH = detPSTH;
    opt.dim = [1 2];
    opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
    opt.bin_num = 20;
    psth = calc_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, opt);

    % show flow field
    clear pltopt;
    pltopt.min_sample = 1000;
    pltopt.trj_color = [
        0.00  0.45  0.74;
        0.85  0.33  0.10];
    pltopt.legend = {'Choice 1', 'Choice 2'};
    pltopt.xlabel = 'PC 1';
    pltopt.ylabel = 'PC 2';
    [fh, stat] = show_flow_field(psth, pltopt);
    title(sprintf('Pair %d', TASK_SET_ID))
    print(fh, '-dpdf', fullfile(FigDir, 'PCA', sprintf('flow_field_pair%d_%s_%s_session%d.pdf', TASK_SET_ID, monkey, experiment, n)))

    % smoothness
    u{n} = stat.U; v{n} = stat.V;
    x{n} = stat.Xgrid; y{n} = stat.Ygrid;
end
save(fullfile(InterimDir, sprintf('stat_pair%d_%s_%s.mat', TASK_SET_ID, monkey, experiment)), 'u', 'v', 'x', 'y')

%% check smoothness
roughness_score = nan(n_files, 2);
for TASK_SET_ID = 1:2
    load(fullfile(InterimDir, sprintf('stat_pair%d_%s_%s.mat', TASK_SET_ID, monkey, experiment)))
    roughness_score(:, TASK_SET_ID) = get_roughness(u,v,x,y);
end

opt.plot = set_plot_opt('roma', 2);
figure('Position', [100 100 220 150]); hold on
for t = 1:2
    plot(roughness_score(:,t), '.-', 'Color', opt.plot.color(t,:))
end
legend({'Pair 1', 'Pair 2'})
xlabel('#Session'); ylabel('Roughness')

%% check dPCA axis
for n = 1:n_files
    % get dPCA axis
    data = load(fullfile(MainInterimDir, 'dPCA\cat-choice\', sprintf('step2_dPCA_%s_%s_session%d.mat', monkey, experiment, n))).data;

    target = {{'Choice', 1}, {'Category', 1}};
    ntarg = length(target);
    tdim = nan(ntarg,1);
    for ta = 1:ntarg
        I = find(strcmp(target{ta}{1}, data.margNames));
        if isempty(I)
            error('No such dimension: %s', target{ta}{1});
        end
        tmpdim = find(data.whichMarg == I, target{ta}{2});
        tdim(ta) = tmpdim(end);
    end

    coef = data.W(:, tdim);

    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % check activity
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    ID = targ_cor;
    ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:));
    trial_classifier_result(ID, {'targ_cor', 'targ_cho'}, {targ_cor, targ_cho})

    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));
    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    psth_choice = coef(:,1)' * psth(:,:); psth_choice = reshape(psth_choice, ntime, ncond);
    psth_cat = coef(:,2)' * psth(:,:); psth_cat = reshape(psth_cat, ntime, ncond);

    tstamp = fnd.tstamp{1};
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
    figure('Position', [100 100 450 150]);
    subplot(1,2,1); hold on
    for c = 1:ncond
        plot(tstamp, psth_choice(:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
    end
    title('Choice axis'); xlabel('Time (ms)')
    subplot(1,2,2); hold on
    for c = 1:ncond
        plot(tstamp, psth_cat(:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
    end
    title('Category axis'); xlabel('Time (ms)')

    figure('Position', [100 400 300 250]); hold on
    for c = 1:ncond
        plot(psth_choice(:,c), psth_cat(:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
    end
    xlabel('Choice axis'); ylabel('Category axis'); xlim([-0.05 0.05]); ylim([-0.03 0.03])

end


%% check PCA axis
for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % check activity
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    ID = targ_cor;
    ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:));
    trial_classifier_result(ID, {'targ_cor', 'targ_cho'}, {targ_cor, targ_cho})

    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));
    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % psth_pca = nan(nunit, ntime, 2);
    % for c = 1:2
    %     psth_pca(:,:,c) = mean(r_detrend(:,:,targ_cho(1,:)==c), 3);
    % end
    psth_pca = psth;
    [coeff, score, latent] = pca(psth_pca(:,:)');
    coef = coeff(:,1:2);

    psth_choice = coef(:,1)' * psth(:,:); psth_choice = reshape(psth_choice, ntime, ncond);
    psth_cat = coef(:,2)' * psth(:,:); psth_cat = reshape(psth_cat, ntime, ncond);

    tstamp = fnd.tstamp{1};
    opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
    figure('Position', [100 100 450 150]);
    subplot(1,2,1); hold on
    for c = 1:ncond
        plot(tstamp, psth_choice(:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
    end
    title('PC 1'); xlabel('Time (ms)')
    subplot(1,2,2); hold on
    for c = 1:ncond
        plot(tstamp, psth_cat(:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
    end
    title('PC 2'); xlabel('Time (ms)')

    figure('Position', [100 400 300 250]); hold on
    for c = 1:ncond
        plot(psth_choice(:,c), -psth_cat(:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
    end
    xlabel('PC 1'); ylabel('PC 2'); xlim([-0.05 0.05]); ylim([-0.03 0.03])

end

%% PC trajectories in all sessions
[psth_choice, psth_cat] = deal(cell(n_files, 2));
for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % check activity
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    task_set = fnd.getp('task_set');
    for t = 1:2
        ID = targ_cho;
        ID(targ_cho~=targ_cor) = NaN; ID(task_set~=t) = NaN;
        trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set})

        psth = fnd.PSTH(ID, {'boxcar', 100}, [], true); psth = psth{1}; % (unit, time, condition)

        psth_pca = psth;
        [coeff, score, latent] = pca(psth_pca(:,:)');
        coef = coeff(:,1:2);

        psth_choice{n,t} = coef(:,1)' * psth(:,:); psth_choice{n,t} = reshape(psth_choice{n,t}, ntime, ncond);
        psth_cat{n,t} = coef(:,2)' * psth(:,:); psth_cat{n,t} = reshape(psth_cat{n,t}, ntime, ncond);
    end
end

%%
for n = 1:n_files
    for t = 1:2
        if mean(psth_choice{n,t}(1:400,2))<0
            psth_choice{n,t} = -psth_choice{n,t};
        end
        if mean(psth_cat{n,t}(1:400,2))<0
            psth_cat{n,t} = -psth_cat{n,t};
        end
    end
end

tstamp = fnd.tstamp{1};
opt.plot = set_plot_opt('roma', 2);
figure('Position', [100 100 1000 1000]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for t = 1:2
        c = 1;
        h(t) = plot(psth_choice{n,t}(:,c)/1000, psth_cat{n,t}(:,c)/1000, 'Color', opt.plot.color(t,:), 'LineStyle', '-');
        c = 2;
        plot(psth_choice{n,t}(:,c)/1000, psth_cat{n,t}(:,c)/1000, 'Color', opt.plot.color(t,:), 'LineStyle', '--')
    end
    if n==1; legend(h, {'Pair 1', 'Pair 2'}); end
    title(sprintf('Session %d', n))
end
% xlabel('PC 1'); ylabel('PC 2'); xlim([-0.05 0.05]); ylim([-0.03 0.03])
format_panel(gcf, 'xlabel', 'PC 1', 'ylabel', 'PC 2', 'xlim', [-0.06 0.06], 'ylim', [-0.06 0.06])






