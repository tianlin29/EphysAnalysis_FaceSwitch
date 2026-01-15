run('../Initialize.m');
monkey = 'Woody'; % Nick, Woody
experiment = 'learnTask4'; % learnTask2, learnTask3, learnTask4, faceColor, faceColor_passiveLong, threeExemplar
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'distance'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'distance'); mkdir(InterimDir);

%% compare dist, proj, angle between category axes (average)
[angle, dist] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % perform PCA on psth
    tstamp = fnd.tstamp{1};
    tmp = psth(:, tstamp>=200 & tstamp<=500, :);
    npc = nunit;
    if npc<nunit
        [coeff, score, latent] = pca(psth(:,:)');
        score = reshape(coeff' * psth(:,:), nunit, ntime, []);
        psth = score(1:npc, :, :);
    end

    % get angle
    vector_pair1_cat = psth(:,:,1) - psth(:,:,3); % pointing to cat 1
    vector_pair2_cat = psth(:,:,2) - psth(:,:,4); % pointing to cat 1
    angle{n} = nan(ntime, 1);
    for t = 1:ntime
        cos_theta = dot(vector_pair1_cat(:,t), vector_pair2_cat(:,t)) / norm(vector_pair1_cat(:,t)) / norm(vector_pair2_cat(:,t));
        cos_theta = max(min(cos_theta, 1), -1);
        angle{n}(t) = rad2deg(acos(cos_theta));
    end

    % get dist
    mid_point_pair1_cat = (psth(:,:,1) + psth(:,:,3)) / 2;
    mid_point_pair2_cat = (psth(:,:,2) + psth(:,:,4)) / 2;
    dist{n} = nan(ntime, 1);
    for t = 1:ntime
        dist{n}(t) = norm(mid_point_pair1_cat(:,t) - mid_point_pair2_cat(:,t));
    end
end


angle_all = cell2mat(angle'); angle_mn = mean(angle_all(tstamp>=200,:), 1);
dist_all = cell2mat(dist'); dist_mn = mean(dist_all(tstamp>=200,:), 1);
save(fullfile(InterimDir, sprintf('angle_and_dist_%s_%s.mat', monkey, experiment)), ...
    'angle', 'dist', 'angle_mn', 'dist_mn')

load(fullfile(InterimDir, sprintf('angle_and_dist_%s_%s.mat', monkey, experiment)))

% angle
opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 570 200]);
subplot(1,2,1); hold on
for n = 1:n_files
    plot(-100:700, angle{n}, 'Color', opt.plot.color(n,:))
end
format_panel(gca, 'xlabel', 'Time (ms)', 'xlim', [-100 700], 'ylabel', 'Angle (deg)')
subplot(1,2,2); hold on
imagesc(-100:700, 1:n_files, cell2mat(angle')', [0 40]); colorbar
set(gca, 'YDir', 'reverse')
format_panel(gca, 'xlabel', 'Time (ms)', 'xlim', [-100 700], 'ylabel', '#Session', 'ylim', [1 n_files])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_%dPC_%s_%s.pdf', npc, monkey, experiment)))

% dist
opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 570 200]);
subplot(1,2,1); hold on
for n = 1:n_files
    plot(-100:700, dist{n}, 'Color', opt.plot.color(n,:))
end
format_panel(gca, 'xlabel', 'Time (ms)', 'xlim', [-100 700], 'ylabel', 'Distance')
subplot(1,2,2); hold on
imagesc(-100:700, 1:n_files, cell2mat(dist')')
set(gca, 'YDir', 'reverse')
format_panel(gca, 'xlabel', 'Time (ms)', 'xlim', [-100 700], 'ylabel', '#Session', 'ylim', [1 n_files])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('distance_%dPC_%s_%s.pdf', npc, monkey, experiment)))

%% paper format
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

[angle_mn, dist_mn] = deal(cell(2, 3)); % {monkey, experiment}
for monkey_id = 1:length(monkey_list)
    for exp_id = 1:length(experiment_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        data = load(fullfile(InterimDir, sprintf('angle_and_dist_%s_%s.mat', monkey, experiment)));
        angle_mn{monkey_id, exp_id} = data.angle_mn;
        dist_mn{monkey_id, exp_id} = data.dist_mn;
    end
end

% plot
color_list = [44 145 224;
    58 191 153;
    240 169 58]/255;

figure('Position', [50 100 300*1.2 200*1.2]);
for exp_id = 1:length(experiment_list)
    subplot(2,3,exp_id); hold on
    nses = max([length(angle_mn{1, exp_id}), length(angle_mn{2, exp_id})]);
    for monkey_id = 1:length(monkey_list)
        if exp_id==1 && monkey_id==2; continue; end
        
        if monkey_id==1
            plot(angle_mn{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7)
        else
            plot(angle_mn{monkey_id, exp_id}, '-', 'Color', color_list(exp_id,:), 'LineWidth', 1);
            scatter(1:length(angle_mn{monkey_id, exp_id}), angle_mn{monkey_id, exp_id}, 7, 'filled', 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', color_list(exp_id,:))
        end
    end
    format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Angle (deg)', ...
        'xlim', [0.5 nses+0.5], ...
        'xtick', 1:2:nses, ...
        'ylim', [20 60])
    xtickangle(0)
    title(sprintf('Pair 1 vs. %d', exp_id+1))
end
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle.pdf')))


figure('Position', [50 100 300*1.2 200*1.2]);
for exp_id = 1:length(experiment_list)
    subplot(2,3,exp_id); hold on
    nses = max([length(dist_mn{1, exp_id}), length(dist_mn{2, exp_id})]);
    for monkey_id = 1:length(monkey_list)
        if exp_id==1 && monkey_id==2; continue; end
        
        if monkey_id==1
            plot(dist_mn{monkey_id, exp_id}, '.-', 'Color', color_list(exp_id,:), 'LineWidth', 0.5, 'MarkerSize', 7)
        else
            plot(dist_mn{monkey_id, exp_id}, '-', 'Color', color_list(exp_id,:), 'LineWidth', 1);
            scatter(1:length(dist_mn{monkey_id, exp_id}), dist_mn{monkey_id, exp_id}, 7, 'filled', 'o', ...
                'MarkerEdgeColor', 'black', 'MarkerFaceColor', color_list(exp_id,:))
        end
    end
    format_panel(gca, 'xlabel', '#Session', 'ylabel', 'Distance', ...
        'xlim', [0.5 nses+0.5], ...
        'xtick', 1:2:nses, ...
        'ylim', [0.015 0.07])
    xtickangle(0)
    title(sprintf('Pair 1 vs. %d', exp_id+1))
end
print(gcf, '-dpdf', fullfile(FigDir, sprintf('distance.pdf')))


%% compare dist, proj, angle between category axes (single-trial)
[dist, proj, angle, targ_cho_pair2_cat1] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get dist
    psth_pair1_cat1 = psth(:,:,1);
    r_pair2_cat1 = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1); ntrial = size(r_pair2_cat1, 3);
    targ_cho_pair2_cat1{n} = targ_cho(1,task_set(1,:)==2 & targ_cor(1,:)==1);
    dist{n} = nan(ntime, ntrial);
    for t = 1:ntime
        dist{n}(t,:) = vecnorm(psth_pair1_cat1(:,t) - r_pair2_cat1(:,t,:), 2, 1);
    end

    % get proj
    vector_pair1_cat = psth(:,:,1) - psth(:,:,3); % pointing to cat 1
    r_pair2_cat1 = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1); ntrial = size(r_pair2_cat1, 3);
    % targ_cho_pair2_cat1{n} = targ_cho(1,task_set(1,:)==2 & targ_cor(1,:)==1);
    proj{n} = nan(ntime, ntrial);
    for t = 1:ntime
        proj{n}(t,:) = vector_pair1_cat(:,t)' * squeeze(r_pair2_cat1(:,t,:));
    end

    % get angle
    vector_pair1_cat = psth(:,:,3) - psth(:,:,1); % (unit, time)
    r_pair2_cat1 = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1); ntrial = size(r_pair2_cat1, 3);
    vector_pair2_cat = psth(:,:,4) - r_pair2_cat1; % (unit, time, trial)
    % targ_cho_pair2_cat1{n} = targ_cho(1,task_set(1,:)==2 & targ_cor(1,:)==1);
    angle{n} = nan(ntime, ntrial);
    for t = 1:ntime
        for tr = 1:ntrial
            cos_theta = dot(vector_pair1_cat(:,t), vector_pair2_cat(:,t,tr)) / norm(vector_pair1_cat(:,t)) / norm(vector_pair2_cat(:,t,tr));
            cos_theta = max(min(cos_theta, 1), -1);
            angle{n}(t,tr) = rad2deg(acos(cos_theta));
        end
    end
end

save(fullfile(InterimDir, sprintf('compare_%s_%s.mat', monkey, experiment)), ...
    'dist', 'proj', 'angle', 'targ_cho_pair2_cat1')


% dist
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    imagesc([dist{n}(:,targ_cho_pair2_cat1{n}==1), dist{n}(:,targ_cho_pair2_cat1{n}==2)]', [0.1 0.4])
    plot(xlim, [sum(targ_cho_pair2_cat1{n}==1), sum(targ_cho_pair2_cat1{n}==1)], 'r', 'LineWidth', 1)
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,0,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Distance')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('dist_single_trial_%s_%s.pdf', monkey, experiment)))


figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(1:801, ...
        mean(dist{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(dist{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    plot_trace(1:801, ...
        mean(dist{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(dist{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);

    % plot(mean(dist{n}(:,targ_cho_pair2_cat1{n}==1), 2))
    % plot(mean(dist{n}(:,targ_cho_pair2_cat1{n}==2), 2))
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Distance', 'xlim', [1 801])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('dist_correct_vs_wrong_%s_%s.pdf', monkey, experiment)))

opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 450 250]);
subplot(1,2,1); hold on; title({'Correct', 'Red: session 1  Blue: last session'})
for n = 1:n_files
    plot_trace(1:801, ...
        mean(dist{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(dist{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        opt.plot.color(n,:));
end
subplot(1,2,2); hold on; title('Error')
for n = 1:n_files
    plot_trace(1:801, ...
        mean(dist{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(dist{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        opt.plot.color(n,:));
end
format_panel(gcf, 'lim_match', {0,1,0}, 'xlabel', 'Time', 'ylabel', 'Distance', 'xlim', [1 801])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('dist_session_%s_%s.pdf', monkey, experiment)))



% proj
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    imagesc([proj{n}(:,targ_cho_pair2_cat1{n}==1), proj{n}(:,targ_cho_pair2_cat1{n}==2)]', [-0.016 0.014])
    plot(xlim, [sum(targ_cho_pair2_cat1{n}==1), sum(targ_cho_pair2_cat1{n}==1)], 'r', 'LineWidth', 1)
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,0,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Projection')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('proj_single_trial_%s_%s.pdf', monkey, experiment)))

figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(1:801, ...
        mean(proj{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(proj{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    plot_trace(1:801, ...
        mean(proj{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(proj{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);

    % plot(mean(proj{n}(:,targ_cho_pair2_cat1{n}==1), 2))
    % plot(mean(proj{n}(:,targ_cho_pair2_cat1{n}==2), 2))
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Projection', 'xlim', [1 801])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('proj_correct_vs_wrong_%s_%s.pdf', monkey, experiment)))

opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 450 250]);
subplot(1,2,1); hold on; title({'Correct', 'Red: session 1  Blue: last session'})
for n = 1:n_files
    plot_trace(1:801, ...
        mean(proj{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(proj{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        opt.plot.color(n,:));
end
subplot(1,2,2); hold on; title('Error')
for n = 1:n_files
    plot_trace(1:801, ...
        mean(proj{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(proj{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        opt.plot.color(n,:));
end
format_panel(gcf, 'lim_match', {0,1,0}, 'xlabel', 'Time', 'ylabel', 'Projection', 'xlim', [1 801])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('proj_session_%s_%s.pdf', monkey, experiment)))





% angle
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    imagesc([angle{n}(:,targ_cho_pair2_cat1{n}==1), angle{n}(:,targ_cho_pair2_cat1{n}==2)]', [50 100])
    plot(xlim, [sum(targ_cho_pair2_cat1{n}==1), sum(targ_cho_pair2_cat1{n}==1)], 'r', 'LineWidth', 1)
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,0,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Angle')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_single_trial_%s_%s.pdf', monkey, experiment)))

figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(1:801, ...
        mean(angle{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(angle{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    plot_trace(1:801, ...
        mean(angle{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(angle{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);

    % plot(mean(angle{n}(:,targ_cho_pair2_cat1{n}==1), 2))
    % plot(mean(angle{n}(:,targ_cho_pair2_cat1{n}==2), 2))
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Angle', 'xlim', [1 801])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_correct_vs_wrong_%s_%s.pdf', monkey, experiment)))

opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 450 250]);
subplot(1,2,1); hold on; title({'Correct', 'Red: session 1  Blue: last session'})
for n = 1:n_files
    plot_trace(1:801, ...
        mean(angle{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(angle{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        opt.plot.color(n,:));
end
subplot(1,2,2); hold on; title('Error')
for n = 1:n_files
    plot_trace(1:801, ...
        mean(angle{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(angle{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        opt.plot.color(n,:));
end
format_panel(gcf, 'lim_match', {0,1,0}, 'xlabel', 'Time', 'ylabel', 'Angle', 'xlim', [1 801])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_session_%s_%s.pdf', monkey, experiment)))


%% regression of difficulty (待删)
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);

% get smooth detrended raster
r = fnd.raster(1); r = r{1}; % (unit, time, trial)
r_mn = mean(r, 3);
r_detrend = r - r_mn;

kernel = fspecial('average', [1, 100]);
r_detrend = nanconv(r_detrend, kernel, 'same');

% get ID
morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:)'; morph_level = abs(morph_level);
task_set = fnd.getp('task_set'); task_set = task_set(1,:)';
targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:)';
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:)';
trial = (1:size(morph_level,1))'; mod_trial = mod(trial, 2);

% regression
training = r_detrend(:,:,task_set==1 & correct==1);
train_value = morph_level(task_set==1 & correct==1);

tstamp = fnd.tstamp{1};
t_list = [linspace(min(tstamp)+50, max(tstamp)-50, 20)' - 50, linspace(min(tstamp)+50, max(tstamp)-50, 20)' + 50];
ntime = size(t_list,1);

ntrial = length(train_value);
test_value = nan(ntrial, ntime);
for t = 1:ntime
    t
    I = tstamp>=t_list(t,1) & tstamp<=t_list(t,2);
    trial_id = (1:length(train_value))';
    for tr = 1:ntrial
        W = regress(train_value(trial_id~=tr), [squeeze(mean(training(:,I,trial_id~=tr),2))', ones(size(training,3)-1,1)]); % (trial, unit+1)
        weight = W(1:end-1);
        offset = W(end);

        testing = squeeze(mean(training(:,I,trial_id==tr),2))';
        test_value(tr,t) = testing * weight + offset;
    end
end

ucoh = unique(train_value);
ntrial = arrayfun(@(x) sum(train_value==x), ucoh);
ucoh_pair1 = ucoh(ntrial>=10);
ncond = length(ucoh_pair1);

test_value_mn = nan(ncond, ntime);
for c = 1:ncond
    test_value_mn(c,:) = mean(test_value(train_value==ucoh_pair1(c),:), 1);
end

opt.plot = set_plot_opt('vik', ncond);
figure; hold on
for c = 1:ncond
    plot([1 20], [ucoh_pair1(c) ucoh_pair1(c)], '--', 'Color', opt.plot.color(c,:))
    plot(test_value_mn(c,:), 'Color', opt.plot.color(c,:), 'LineWidth', 1.5)
end


%% difficulty axis of pair 1 itself (待删)
[angle] = deal(cell(n_files, 1));
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:); morph_level = abs(morph_level);
    task_set = fnd.getp('task_set'); task_set = task_set(1,:);
    targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:);
    targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
    correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:);
    trial = 1:size(morph_level,2);

    ucoh = unique(morph_level(task_set==1));
    ntrial = arrayfun(@(x) sum(morph_level(task_set==1)==x), ucoh);
    ucoh_pair1 = ucoh(ntrial>=10);

    ID = nan(size(morph_level));
    ID(morph_level<median(ucoh_pair1) & task_set==1 & mod(trial,2)==0) = 1;
    ID(morph_level>median(ucoh_pair1) & task_set==1 & mod(trial,2)==0) = 2;
    ID(morph_level<median(ucoh_pair1) & task_set==1 & mod(trial,2)==1) = 3;
    ID(morph_level>median(ucoh_pair1) & task_set==1 & mod(trial,2)==1) = 4;
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'morph_level', 'correct'}, {task_set, targ_cor, morph_level, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID==c), 3);
    end

    % angle between difficulty axes
    tstamp = fnd.tstamp{1};
    t_list = [linspace(min(tstamp)+50, max(tstamp)-50, 20)' - 50, linspace(min(tstamp)+50, max(tstamp)-50, 20)' + 50];
    ntime = size(t_list,1);

    vector_pair1 = psth(:,:,2) - psth(:,:,1); % (unit, time)
    vector_pair2 = psth(:,:,4) - psth(:,:,3);
    angle{n} = nan(ntime, 1);
    for t = 1:ntime
        I = tstamp>=t_list(t,1) & tstamp<=t_list(t,2);
        cos_theta = dot(mean(vector_pair1(:,I),2), mean(vector_pair2(:,I),2)) / norm(mean(vector_pair1(:,I),2)) / norm(mean(vector_pair2(:,I),2));
        cos_theta = max(min(cos_theta, 1), -1);
        angle{n}(t) = rad2deg(acos(cos_theta));
    end
end



figure;
imagesc(cell2mat(angle')'); colorbar


opt.plot = set_plot_opt('vik', n_files);
figure; hold on
for n = 1:n_files
    plot(angle{n}, 'Color', opt.plot.color(n,:))
end


%% difficulty axis
[angle_correct, angle_wrong] = deal(cell(n_files, 1));
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:); morph_level = abs(morph_level);
    task_set = fnd.getp('task_set'); task_set = task_set(1,:);
    targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:);
    targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
    correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:);

    ucoh = unique(morph_level(task_set==1));
    ntrial = arrayfun(@(x) sum(morph_level(task_set==1)==x), ucoh);
    ucoh_pair1 = ucoh(ntrial>=0);

    ucoh = unique(morph_level(task_set==2 & correct==1));
    ntrial = arrayfun(@(x) sum(morph_level(task_set==2 & correct==1)==x), ucoh);
    ucoh_pair2_correct = ucoh(ntrial>=0);

    ucoh = unique(morph_level(task_set==2 & correct==0));
    ntrial = arrayfun(@(x) sum(morph_level(task_set==2 & correct==0)==x), ucoh);
    ucoh_pair2_wrong = ucoh(ntrial>=0);

    ID = nan(size(morph_level));
    ID(morph_level<median(ucoh_pair1) & task_set==1 & correct==1) = 1;
    ID(morph_level>median(ucoh_pair1) & task_set==1 & correct==1) = 2;
    ID(morph_level<median(ucoh_pair2_correct) & task_set==2 & correct==1) = 3;
    ID(morph_level>median(ucoh_pair2_correct) & task_set==2 & correct==1) = 4;
    ID(morph_level<median(ucoh_pair2_wrong) & task_set==2 & correct==0) = 5;
    ID(morph_level>median(ucoh_pair2_wrong) & task_set==2 & correct==0) = 6;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'morph_level', 'correct'}, {task_set, targ_cor, morph_level, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID==c), 3);
    end

    % angle between difficulty axes
    vector_pair1 = psth(:,:,1) - psth(:,:,2); % (unit, time)
    vector_pair2_correct = psth(:,:,3) - psth(:,:,4);
    vector_pair2_wrong = psth(:,:,5) - psth(:,:,6);

    [angle_correct{n}, angle_wrong{n}] = deal(nan(ntime, 1));
    for t = 1:ntime
        cos_theta = dot(vector_pair1(:,t), vector_pair2_correct(:,t)) / norm(vector_pair1(:,t)) / norm(vector_pair2_correct(:,t));
        cos_theta = max(min(cos_theta, 1), -1);
        angle_correct{n}(t) = rad2deg(acos(cos_theta));

        cos_theta = dot(vector_pair1(:,t), vector_pair2_wrong(:,t)) / norm(vector_pair1(:,t)) / norm(vector_pair2_wrong(:,t));
        cos_theta = max(min(cos_theta, 1), -1);
        angle_wrong{n}(t) = rad2deg(acos(cos_theta));
    end
end


save(fullfile(InterimDir, sprintf('difficulty_axis_%s_%s.mat', monkey, experiment)), ...
    'angle_correct', 'angle_wrong')


load(fullfile(InterimDir, sprintf('difficulty_axis_%s_%s.mat', monkey, experiment)))

figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot(-100:700, angle_correct{n}, 'Color', 'black')
    plot(-100:700, angle_wrong{n}, 'Color', 'red')
    plot([-100 700], [90 90], ':', 'Color', 'black')
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Angle', 'xlim', [-100 700], 'ylim', [40 110])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_difficulty_axis_%s_%s.pdf', monkey, experiment)))



opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 300 250]); hold on
for n = 1:n_files
    plot(-100:700, angle_correct{n}, 'Color', opt.plot.color(n,:))
end
plot([-100 700], [90 90], ':', 'Color', 'black')
format_panel(gcf, 'xlabel', 'Time', 'ylabel', 'Angle', 'xlim', [-100 700], 'ylim', [40 110])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_difficulty_axis_session_%s_%s.pdf', monkey, experiment)))


figure('Position', [100 100 300 250]); hold on
imagesc(cell2mat(angle_correct')', [40 110]); colorbar
format_panel(gcf, 'xlabel', 'Time', 'ylabel', 'Session')
set(gca, 'YDir', 'reverse');
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_difficulty_axis_heatmap_%s_%s.pdf', monkey, experiment)))











%% Procrustes transformation (300-400 ms + single-trial)
[theta, scale, abs_translation, dist, targ_cho_pair2_cat1] = deal(cell(n_files, 1));
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    % trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end
    
    % get PCA
    npc = 2;
    if npc==nunit
        psth_pc = psth; % (unit, time, condition)
        r_detrend_pc = r_detrend;
    else
        [coeff, score, latent] = pca(psth(:,:)');
        psth_pc = reshape(score', nunit, ntime, ncond);
        psth_pc = psth_pc(1:npc,:,:); % (unit/pc, time, condition)

        r_detrend_pc = coeff' * r_detrend(:,:);
        r_detrend_pc = reshape(r_detrend_pc, nunit, ntime, []);
        r_detrend_pc = r_detrend_pc(1:npc,:,:); % (unit/pc, time, condition)
    end
    
    % get single-trial transformation
    tstamp = fnd.tstamp{1};
    t_list = [linspace(min(tstamp)+50, max(tstamp)-50, 20)' - 50, linspace(min(tstamp)+50, max(tstamp)-50, 20)' + 50];
    ntime = size(t_list,1);
    for t = 1:ntime
        I = tstamp>=t_list(t,1) & tstamp<=t_list(t,2);
        points_pair1_cat1 = mean(psth_pc(:,I,1), 2)'; % (condition, unit)
        points_pair1_cat2 = mean(psth_pc(:,I,3), 2)';
        points_pair2_cat1 = permute(mean(r_detrend_pc(:,I,task_set(1,:)==2 & targ_cor(1,:)==1), 2), [3 1 2]); % (trial, unit)
        targ_cho_pair2_cat1{n} = targ_cho(1,task_set(1,:)==2 & targ_cor(1,:)==1);
        points_pair2_cat2 = mean(psth_pc(:,I,4), 2)';

        ntrial = size(points_pair2_cat1, 1);
        for tr = 1:ntrial
            [~, Z, transform] = procrustes([points_pair1_cat1(:,:); points_pair1_cat2(:,:)], [points_pair2_cat1(tr,:); points_pair2_cat2(:,:)], ...
                'scaling', false, 'reflection', false); % (condition, unit) -> (point, dim)
            if npc==2; theta{n}(t,tr) = atan2d(transform.T(2,1),transform.T(1,1)); end
            scale{n}(t,tr) = transform.b;
            abs_translation{n}(t,tr) = norm(transform.c(1,:));
            dist{n}(t,tr) = norm(mean([points_pair1_cat1(:,:); points_pair1_cat2(:,:)]) - mean([points_pair2_cat1(tr,:); points_pair2_cat2(:,:)]));
        end
    end
end

save(fullfile(InterimDir, sprintf('Procrustes_%s_%s.mat', monkey, experiment)), ...
    'theta', 'scale', 'abs_translation', 'dist', 'targ_cho_pair2_cat1')

load(fullfile(InterimDir, sprintf('Procrustes_%s_%s.mat', monkey, experiment)))

% plot dots
n = 15;
figure;
the = abs_theta{n}(4,targ_cho_pair2_cat1{n}==1);
rrr = abs_translation{n}(10,targ_cho_pair2_cat1{n}==1);
polarplot(the, rrr, 'o', 'Color', 'black', 'LineWidth', 2)

figure;
the = abs_theta{n}(4,targ_cho_pair2_cat1{n}==2);
rrr = abs_translation{n}(10,targ_cho_pair2_cat1{n}==2);
polarplot(the, rrr, 'o', 'Color', 'red', 'LineWidth', 2)



% translation
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(1:ntime, ...
        mean(abs_translation{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(abs_translation{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    plot_trace(1:ntime, ...
        mean(abs_translation{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(abs_translation{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);
    title(sprintf('session %d', n))

    p = nan(20, 1);
    for t = 1:20
        p(t) = run_ttest2(abs_translation{n}(t,targ_cho_pair2_cat1{n}==1), abs_translation{n}(t,targ_cho_pair2_cat1{n}==2), '=', 'off');
    end
    p_plot = find(p<0.05);
    scatter(p_plot, 0.08*ones(size(p_plot)), 7, 'black', '*')
end
format_panel(gcf, 'xlable', 'Time', 'ylabel', 'Translation', 'ylim', [0 0.1])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_translation_%s_%s.pdf', monkey, experiment)))



% difference in distance and translation
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot(mean(dist{n}(:,targ_cho_pair2_cat1{n}==2), 2) - mean(dist{n}(:,targ_cho_pair2_cat1{n}==1), 2))
    plot(mean(abs_translation{n}(:,targ_cho_pair2_cat1{n}==2), 2) - mean(abs_translation{n}(:,targ_cho_pair2_cat1{n}==1), 2))
    if n==1; legend({'Distance', 'Translation'}); end
    plot(xlim, [0 0], ':', 'Color', [0 0 0])
end
format_panel(gcf, 'xlabel', 'Time', 'ylabel', {'Difference between', 'correct and error'})


% theta
abs_theta = cellfun(@(x) abs(x), theta, 'uni', 0);
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(1:ntime, ...
        mean(abs_theta{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(abs_theta{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    plot_trace(1:ntime, ...
        mean(abs_theta{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(abs_theta{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);
    plot(xlim, [90 90], ':', 'Color', [0 0 0])
    title(sprintf('session %d', n))

    p = nan(20, 1);
    for t = 1:20
        p(t) = run_ttest2(abs_theta{n}(t,targ_cho_pair2_cat1{n}==1), abs_theta{n}(t,targ_cho_pair2_cat1{n}==2), '=', 'off');
    end
    p_plot = find(p<0.05);
    scatter(p_plot, 120*ones(size(p_plot)), 7, 'black', '*')
end
format_panel(gcf, 'xlable', 'Time', 'ylabel', 'Theta')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_theta_%s_%s.pdf', monkey, experiment)))

% dist
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(1:ntime, ...
        mean(dist{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(dist{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    plot_trace(1:ntime, ...
        mean(dist{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(dist{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);
    title(sprintf('session %d', n))

    p = nan(20, 1);
    for t = 1:20
        p(t) = run_ttest2(dist{n}(t,targ_cho_pair2_cat1{n}==1), dist{n}(t,targ_cho_pair2_cat1{n}==2), '=', 'off');
    end
    p_plot = find(p<0.05);
    scatter(p_plot, 0.08*ones(size(p_plot)), 7, 'black', '*')
end
format_panel(gcf, 'xlable', 'Time', 'ylabel', 'Distance', 'ylim', [0 0.1])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_direct_distance_%s_%s.pdf', monkey, experiment)))


%% Procrustes transformation (300-400 ms)
% [scale, theta, abs_translation, targ_cho_pair2_cat1, psth_register, psth_new] = deal(cell(n_files, 1));
[scale, theta, abs_translation] = deal(nan(n_files, 1));
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    % trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end
    
    % get PCA
    npc = 2;
    if npc==nunit
        psth_pc = psth; % (unit, time, condition)
        r_detrend_pc = r_detrend;
    else
        [coeff, score, latent] = pca(psth(:,:)');
        psth_pc = reshape(score', nunit, ntime, ncond);
        psth_pc = psth_pc(1:npc,:,:); % (unit/pc, time, condition)
        psth_register{n} = psth_pc;

        r_detrend_pc = coeff' * r_detrend(:,:);
        r_detrend_pc = reshape(r_detrend_pc, nunit, ntime, []);
        r_detrend_pc = r_detrend_pc(1:npc,:,:); % (unit/pc, time, condition)
    end
    
    % get single-trial transformation
    tstamp = fnd.tstamp{1};
    I = tstamp>=300 & tstamp<=400;
    points_pair1_cat1 = mean(psth_pc(:,I,1), 2)'; % (condition, unit)
    points_pair1_cat2 = mean(psth_pc(:,I,3), 2)';
    points_pair2_cat1 = mean(psth_pc(:,I,2), 2)';
    points_pair2_cat2 = mean(psth_pc(:,I,4), 2)';

    [~, Z, transform] = procrustes([points_pair1_cat1(:,:); points_pair1_cat2(:,:)], [points_pair2_cat1(:,:); points_pair2_cat2(:,:)], ...
        'scaling', false, 'reflection', false); % (condition, unit) -> (point, dim)
    if npc==2; theta(n) = atan2d(transform.T(2,1),transform.T(1,1)); end
    scale(n) = transform.b;
    abs_translation(n) = norm(transform.c(1,:));
end


figure('Position', [100 100 740 180]);
subplot(1,3,1); hold on; title('Theta')
plot(abs(theta), '.-')
subplot(1,3,2); hold on; title('Scale')
plot(scale, '.-')
subplot(1,3,3); hold on; title('Translation')
plot(abs_translation, '.-')
format_panel(gcf, 'lim_match', {1,0,0}, 'xlabel', 'Time')










%% [more serious] angle between pair1 cat axis and cat1 pair axis
[angle] = deal(cell(n_files, 1));
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);

    tstamp = fnd.tstamp{1};
    t_list = [linspace(min(tstamp)+50, max(tstamp)-50, 20)' - 50, linspace(min(tstamp)+50, max(tstamp)-50, 20)' + 50];
    ntime_range = size(t_list,1);
    angle{n} = nan(ntime_range, 1);
    for t = 1:ntime_range
        I = tstamp>=t_list(t,1) & tstamp<=t_list(t,2);
        % pair
        ID = task_set;
        ID(~correct) = NaN;
        % trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

        [nunit, ~, ~] = size(r);
        ntime = sum(I);
        ncond = max(ID(:));
        psth = nan(nunit, ntime, ncond);
        for c = 1:ncond
            psth(:,:,c) = mean(r_detrend(:,I,ID(1,:)==c), 3);
        end

        [coeff, score, latent] = pca(psth(:,:)'); % coeff ..(unit, pc)
        pair_axis = coeff(:,1);

        % category 1
        ID = targ_cor;
        ID(task_set~=1) = NaN;
        ID(~correct) = NaN;
        % trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

        [nunit, ~, ~] = size(r);
        ntime = sum(I);
        ncond = max(ID(:));
        psth = nan(nunit, ntime, ncond);
        for c = 1:ncond
            psth(:,:,c) = mean(r_detrend(:,I,ID(1,:)==c), 3);
        end

        [coeff, score, latent] = pca(psth(:,:)'); % coeff ..(unit, pc)
        cat1_axis = coeff(:,1);

        % angle
        cos_theta = dot(pair_axis, cat1_axis) / norm(pair_axis) / norm(cat1_axis);
        cos_theta = max(min(cos_theta, 1), -1);
        angle{n}(t) = rad2deg(acos(cos_theta));
    end
end


% permutation
ndim = 50;
npermute = 100000;
v_min = -1; v_max = 1;
angle_permute = nan(npermute, 1);
for p = 1:npermute
    v1 = v_min + (v_max-v_min)*rand(ndim, 1);
    v2 = v_min + (v_max-v_min)*rand(ndim, 1);
    cos_theta = dot(v1, v2) / norm(v1) / norm(v2);
    cos_theta = max(min(cos_theta, 1), -1);
    angle_permute(p) = rad2deg(acos(cos_theta));
end
angle_boundary = prctile(angle_permute, 0.025*100);

figure('Position', [100 100 200 200]); hold on
histogram(angle_permute, 40, 'Normalization', 'probability', 'EdgeColor', 'None')
plot([angle_boundary angle_boundary], ylim, 'Color', 'red', 'LineWidth', 1)
format_panel(gcf, 'xlabel', 'Angle', 'ylabel', 'Probability', 'xtick', [80 90 100])


% plot
figure('Position', [100 100 200 200]); hold on
opt.plot = set_plot_opt('vik', n_files);
t_middle = mean(t_list, 2);
for n = 1:n_files
    plot(t_middle, angle{n}, '.-', 'Color', opt.plot.color(n,:))
end
plot(xlim, [83 83], 'Color', 'red', 'LineWidth', 1)
format_panel(gcf, 'xlabel', 'Time (ms)', 'ylabel', 'Angle between pair and cat1 axis', 'axis', 'normal', 'xlim', [t_middle(1) t_middle(end)])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('anlge_pair_cat1_%s_%s.pdf', monkey, experiment)))


%% angle between pair1 cat axis and cat1 pair axis
[angle, tstamp] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    tstamp{n} = fnd.tstamp{1};

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get angle
    vector_pair1_cat = psth(:,:,3) - psth(:,:,1); % (unit, time)
    vector_cat1_pair = psth(:,:,2) - psth(:,:,1); % (unit, time)
    angle{n} = nan(ntime, 1);
    for t = 1:ntime
        cos_theta = dot(vector_pair1_cat(:,t), vector_cat1_pair(:,t)) / norm(vector_pair1_cat(:,t)) / norm(vector_cat1_pair(:,t));
        cos_theta = max(min(cos_theta, 1), -1);
        angle{n}(t) = rad2deg(acos(cos_theta));
    end
end

save(fullfile(InterimDir, sprintf('angle_cat_pair_%s_%s.mat', monkey, experiment)), ...
    'angle', 'tstamp')

load(fullfile(InterimDir, sprintf('angle_cat_pair_%s_%s.mat', monkey, experiment)))

angle_all = cell2mat(angle');

figure('Position', [100 100 400 130]);
imagesc(tstamp{1}, 1:n_files, angle_all', [50 90])
colorbar
format_panel(gcf, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', '#Session')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_cat_pair_%s_%s.pdf', monkey, experiment)))



%% projection: 把pair2 cat1分别投射到pair2或者1的cat axis上，一个更能预测cat，一个更能预测choice
[proj_to_pair2, proj_to_pair1, targ_cho_pair2_cat1, tstamp] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    tstamp{n} = fnd.tstamp{1};

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % pair 2 cat 1 projects to pair 2 cat axis -> cat/choice?
    vector_pair2_cat = mean(r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1), 3) - ...
        mean(r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==2), 3); % pointing to cat 1
    % vector_pair2_cat = psth(:,:,2) - psth(:,:,4); % pointing to cat 1
    r_pair2_cat1 = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1); ntrial = size(r_pair2_cat1, 3);
    targ_cho_pair2_cat1{n} = targ_cho(1,task_set(1,:)==2 & targ_cor(1,:)==1);
    proj_to_pair2{n} = nan(ntime, ntrial);
    for t = 1:ntime
        proj_to_pair2{n}(t,:) = vector_pair2_cat(:,t)' * squeeze(r_pair2_cat1(:,t,:));
    end

    % pair 2 cat 1 projects to pair 1 cat axis -> cat/choice?
    vector_pair1_cat = mean(r_detrend(:,:,task_set(1,:)==1 & targ_cor(1,:)==1), 3) - ...
        mean(r_detrend(:,:,task_set(1,:)==1 & targ_cor(1,:)==2), 3); % pointing to cat 1
    % vector_pair1_cat = psth(:,:,1) - psth(:,:,3); % pointing to cat 1
    r_pair2_cat1 = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1); ntrial = size(r_pair2_cat1, 3);
    % targ_cho_pair2_cat1 = targ_cho(1,task_set(1,:)==2 & targ_cor(1,:)==1);
    proj_to_pair1{n} = nan(ntime, ntrial);
    for t = 1:ntime
        proj_to_pair1{n}(t,:) = vector_pair1_cat(:,t)' * squeeze(r_pair2_cat1(:,t,:));
    end
end


save(fullfile(InterimDir, sprintf('proj_to_pair1_or_2_%s_%s.mat', monkey, experiment)), ...
    'proj_to_pair2', 'proj_to_pair1', 'targ_cho_pair2_cat1', 'tstamp')

load(fullfile(InterimDir, sprintf('proj_to_pair1_or_2_%s_%s.mat', monkey, experiment)))

% merge
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    plot_trace(tstamp{n}, ...
        mean(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [44, 145, 224]/255);
    plot_trace(tstamp{n}, ...
        mean(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [150, 200, 240]/255);
    plot([tstamp{n}(1) tstamp{n}(end)], [0 0], ':', 'Color', 'black')
    title(sprintf('session %d', n))

    h1 = plot_trace(tstamp{n}, ...
        mean(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    h2 = plot_trace(tstamp{n}, ...
        mean(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [.5 .5 .5]);
    if n==1; legend([h1 h2], {'Correct', 'Error'}); end
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Projection to pair 1/2 cat axis', 'ylim', [-0.004 0.004])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('proj_to_pair1_or_2_%s_%s.pdf', monkey, experiment)))



% proj to pair 2 cat axis
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    h1 = plot_trace(tstamp{n}, ...
        mean(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    h2 = plot_trace(tstamp{n}, ...
        mean(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(proj_to_pair2{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);
    plot([tstamp{n}(1) tstamp{n}(end)], [0 0], ':', 'Color', 'black')
    if n==1; legend([h1 h2], {'Correct', 'Error'}); end
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Projection to pair 2 cat axis', 'ylim', [-0.004 0.004])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('proj_to_pair2_%s_%s.pdf', monkey, experiment)))


% proj to pair 1 cat axis
figure('Position', [100 100 800 950]);
for n = 1:n_files
    subplot(5,5,n); hold on
    h1 = plot_trace(tstamp{n}, ...
        mean(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==1), 2), ...
        std(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==1), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==1)), ...
        [0 0 0]);
    h2 = plot_trace(tstamp{n}, ...
        mean(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==2), 2), ...
        std(proj_to_pair1{n}(:,targ_cho_pair2_cat1{n}==2), [], 2)./sqrt(sum(targ_cho_pair2_cat1{n}==2)), ...
        [1 0 0]);
    plot([tstamp{n}(1) tstamp{n}(end)], [0 0], ':', 'Color', 'black')
    if n==1; legend([h1 h2], {'Correct', 'Error'}); end
    title(sprintf('session %d', n))
end
format_panel(gcf, 'lim_match', {0,1,0}, 'axis', 'normal', 'xlabel', 'Time', 'ylabel', 'Projection to pair 1 cat axis', 'ylim', [-0.004 0.004])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('proj_to_pair1_%s_%s.pdf', monkey, experiment)))







%% 如果把pair2的试次投射到pair1的category axis上，能否预测该trial的选择？
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);
tstamp = fnd.tstamp{1};

% get smooth detrended raster
r = fnd.raster(1); r = r{1}; r = r(:,tstamp>=0,:); % (unit, time, trial)
r_mn = mean(r, 3);
r_detrend = r - r_mn;

kernel = fspecial('average', [1, 100]);
r_detrend = nanconv(r_detrend, kernel, 'same');

% get ID
task_set = fnd.getp('task_set'); task_set = task_set(1,:)';
targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:)';
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:)';

% linear classification
training = r_detrend(:,:,task_set==2);
train_label = targ_cho(task_set==2); % choice of pair 2
testing = r_detrend(:,:,task_set==2);
test_label = targ_cho(task_set==2); % choice of pair 2

ntime = size(r_detrend,2);
cvmdl = cell(ntime, 1);
accuracy_test = nan(ntime, ntime); % (time of pair 1, time of pair 2)
parfor t = 1:ntime
    cvmdl{t} = fitclinear(squeeze(training(:,t,:))', train_label); % (observation, unit)
end

for t = 1:ntime
    t
    cvmdl_tmp = cvmdl{t};
    y_hat_test = predict(cvmdl_tmp, testing(:,:)'); % (observation, unit)
    ntrial = size(testing, 3);
    y_hat_test = reshape(y_hat_test', ntime, ntrial);
    for tt = 1:ntime
        accuracy_test(t,tt) = mean(y_hat_test(tt,:)'==test_label);
    end
end

% plot
figure('Position', [100 100 300 220]); hold on
imagesc(accuracy_test', [0.35 0.85]); colorbar
plot([1 701], [1 701], 'r', 'LineWidth', 1)
title(sprintf('session %d', n))
xlabel('Train on pair 1'); ylabel('Test on pair 2')


%% visualize single-trial location
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    tstamp = fnd.tstamp{1};

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get PC axis
    psth_time = psth(:,tstamp>=150 & tstamp<=400,:);
    [coeff, score, latent] = pca(psth_time(:,:)'); % coeff ..(unit, pc)

    % project raster to PC axis
    r_time = squeeze(mean(r_detrend(:,tstamp>=150 & tstamp<=400,:), 2));
    r_pc = coeff(:,[1 2])' * r_time;

    cor = correct(1,:);
    [C, ia, ic] = unique(morph_level(1,:));
    opt.plot = set_plot_opt('vik', length(C));

    figure; hold on
    for tr = 1:size(r_pc,2)
        if cor(tr)
            scatter(r_pc(1,tr), r_pc(2,tr), 10, opt.plot.color(ic(tr),:), 'filled')
        else
            scatter(r_pc(1,tr), r_pc(2,tr), 10, opt.plot.color(ic(tr),:), 'filled', 'LineWidth', 1, 'MarkerEdgeColor', 'red')
        end
    end
end





%% why no discovery in single-trial pair signal
a = [0 0];
n_session = 15;
b = [zeros(n_session, 1), linspace(0.054, 0.03, n_session)'];
center = a;
num_points = 1000;
radius = 0.054 * (1:7);

P = cell(1, length(radius));
ddd = nan(num_points, size(b,1), length(radius));
for ra = 1:length(radius)
    for n = 1:size(b,1)
        % sample points
        theta = 2 * pi * rand(num_points, 1);
        r = radius(ra) * sqrt(rand(num_points, 1));
        x = center(1) + r .* cos(theta);
        y = center(2) + r .* sin(theta);
        points = [x, y]; % (point, 2)
        if n==1; P{ra} = points; end

        % get distance
        vector = points - b(n,:);
        ddd(:,n,ra) = vecnorm(vector, 2, 2);
    end
end

figure; histogram(ddd(:,1,1), 'BinEdges', linspace(0, 0.12, 20))

figure('Position', [100 100 280 1000]);
ax = nan(length(radius), 2);
for ra = 1:length(radius)
    ax(ra,1) = subplot(length(radius),2,ra*2-1); hold on
    scatter(P{ra}(:,1), P{ra}(:,2), 1, [.5 .5 .5], 'filled')
    scatter(a(1), a(2), 2, 'red', 'filled')
    scatter(b(:,1), b(:,2), 2, 'blue', 'filled')
    if ra==1; legend({'Sample around pair 1', 'Pair 1', 'Pair 2 moving toward pair 1'}); end

    mn = mean(ddd(:,:,ra), 1);
    ax(ra,2) = subplot(length(radius),2,ra*2); hold on
    plot(mn, '.-')
end
format_panel(ax(:,1), 'lim_match', {1,1,0}, 'axis', 'equal')
format_panel(ax(:,2), 'lim_match', {0,0,0}, 'xlabel', '#Session', 'ylabel', 'Distance to blue points')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('simulation_single_trial_distance.pdf')))


%% decoding using fitclinear
% pair 2, decode category, train on correct, test on correct
figure;
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    tstamp = fnd.tstamp{1};

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; r = r(:,tstamp>=0,:); % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:)';
    task_set = fnd.getp('task_set'); task_set = task_set(1,:)';
    targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:)';
    targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
    correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:)';

    % linear classification
    training = r_detrend(:,:,task_set==2 & correct==1);
    train_label = targ_cor(task_set==2 & correct==1);
    train_morph = abs(morph_level(task_set==2 & correct==1));
    testing = r_detrend(:,:,task_set==2 & correct==0);
    test_label = targ_cor(task_set==2 & correct==0);
    test_morph = abs(morph_level(task_set==2 & correct==0));

    t = 350;

    % correct -> correct
    cvmdl = fitclinear(squeeze(training(:,t,:))', train_label, 'KFold', 5); % (observation, unit)
    y_hat = cvmdl.kfoldPredict;
    cor_CC = y_hat==train_label;

    % correct -> error
    cvmdl = fitclinear(squeeze(training(:,t,:))', train_label); % (observation, unit)
    y_hat = predict(cvmdl, squeeze(testing(:,t,:))'); % (observation, unit)
    cor_CE = y_hat==test_label;

    % error -> error
    cvmdl = fitclinear(squeeze(testing(:,t,:))', test_label, 'KFold', 5); % (observation, unit)
    y_hat = cvmdl.kfoldPredict;
    cor_EE = y_hat==test_label;

    % plot
    [p_CC, pse_CC] = calcGroupMean(cor_CC, train_morph, unique(train_morph), 'binary');
    [p_CE, pse_CE] = calcGroupMean(cor_CE, test_morph, unique(test_morph), 'binary');
    [p_EE, pse_EE] = calcGroupMean(cor_EE, test_morph, unique(test_morph), 'binary');

    opt.color = [0 0 0; 1 0 0; 0 0 1];
    subplot(4,4,n); hold on
    plot(unique(train_morph), p_CC, '.-', 'markers', 7, 'Color', opt.color(1,:));
    % cerrorbar(unique(train_morph), p_CC, pse_CC, 'Color', opt.color(1,:));
    plot(unique(test_morph), p_CE, '.-', 'markers', 7, 'Color', opt.color(2,:));
    % cerrorbar(unique(test_morph), p_CE, pse_CE, 'Color', opt.color(2,:));
    plot(unique(test_morph), p_EE, '.-', 'markers', 7, 'Color', opt.color(3,:));
    % cerrorbar(unique(test_morph), p_EE, pse_EE, 'Color', opt.color(3,:));
    plot(xlim, [.5 .5], ':', 'Color', [0 0 0])
end

format_panel(gcf, 'xlabel', 'Morph', 'ylabel', 'Accuracy')


%% decoding using fitclinear (cross-pair decoding of category)
n = n_files;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);
tstamp = fnd.tstamp{1};

% get smooth detrended raster
r = fnd.raster(1); r = r{1}; r = r(:,tstamp>=0,:); % (unit, time, trial)
r_mn = mean(r, 3);
r_detrend = r - r_mn;

kernel = fspecial('average', [1, 100]);
r_detrend = nanconv(r_detrend, kernel, 'same');

% get ID
task_set = fnd.getp('task_set'); task_set = task_set(1,:)';
targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:)';
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:)';

% linear classification
training = r_detrend(:,:,task_set==1);
train_label = targ_cor(task_set==1);
testing = r_detrend(:,:,task_set==2);
test_label = targ_cor(task_set==2);

ntime = size(r_detrend,2);
cvmdl = cell(ntime, 1);
accuracy_test = nan(ntime, ntime); % (time of pair 1, time of pair 2)
parfor t = 1:ntime
    cvmdl{t} = fitclinear(squeeze(training(:,t,:))', train_label); % (observation, unit)
end

for t = 1:ntime
    t
    cvmdl_tmp = cvmdl{t};
    y_hat_test = predict(cvmdl_tmp, testing(:,:)'); % (observation, unit)
    ntrial = size(testing, 3);
    y_hat_test = reshape(y_hat_test', ntime, ntrial);
    for tt = 1:ntime
        accuracy_test(t,tt) = mean(y_hat_test(tt,:)'==test_label);
    end
end

save(fullfile(InterimDir, sprintf('accuracy_test_%s_%s_session%d.mat', monkey, experiment, n)), 'accuracy_test')

accuracy_test = load(fullfile(InterimDir, sprintf('accuracy_test_%s_%s_session%d.mat', monkey, experiment, n))).accuracy_test;

% plot
figure('Position', [100 100 300 220]); hold on
imagesc(accuracy_test', [0.35 0.85]); colorbar
plot([1 701], [1 701], 'r', 'LineWidth', 1)
title(sprintf('session %d', n))
xlabel('Train on pair 1'); ylabel('Test on pair 2')



%% single-trial trajectory
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    tstamp = fnd.tstamp{1};

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    % kernel = fspecial('average', [1, 100]);
    kernel = fspecial('gaussian', [1, 20*6], 20);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = targ_cor;
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get PC, then project raster onto it
    [coeff, score, latent] = pca(psth(:,:)');
    % coeff ..(unit, pc) for projection
    score = reshape(score', nunit, ntime, ncond); % (pc, time, cond)

    % get single trial
    r_target = r_detrend(:, :, correct(1,:)==1);
    targ_cor_target = targ_cor(:, correct(1,:)==1);
    [nunit, ntime, ntrial] = size(r_target);
    r_proj = r_target(:,:)' * coeff;
    r_score = reshape(r_proj', nunit, ntime, ntrial); % (pc, time, trial)

    % plot
    opt.plot.color = [1 0 0; 0 0 1];
    I = tstamp>=0 & tstamp<=500;
    figure; hold on
    for c = 1:ncond
        % plot(score(1,I,c), 'Color', opt.plot.color(c,:))
        % plot3(tstamp(I), score(1,I,c), score(2,I,c), 'Color', opt.plot.color(c,:))
        plot3(score(1,I,c), score(2,I,c), score(3,I,c), 'Color', opt.plot.color(c,:), 'LineWidth', 1)
    end
    for tr = 1:ntrial
        c = targ_cor_target(1,tr);
        % plot(r_score(1,I,tr), 'Color', opt.plot.color(c,:) * .5 + [.5 .5 .5])
        % plot3(tstamp(I), r_score(1,I,tr), r_score(2,I,tr), 'Color', opt.plot.color(c,:) * .5 + [.5 .5 .5])
        plot3(r_score(1,I,tr), r_score(2,I,tr), r_score(3,I,tr), 'Color', opt.plot.color(c,:) * .5 + [.5 .5 .5])
    end
    grid on
end

%% dissociate distance to angle-related and -unrelated (single trial)
[angle, half_length_pair2, dist_pair_signal, base_length, direct_length, correct_target] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get single trial for test (both correct and wrong trials)
    r_target = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1);
    correct_target{n} = correct(:,task_set(1,:)==2 & targ_cor(1,:)==1);

    % get distance
    dist_cat1 = r_target - psth(:,:,1);
    
    vector_pair1 = psth(:,:,3) - psth(:,:,1);
    vector_pair2 = psth(:,:,4) - r_target;

    ntrial = size(r_target, 3);
    [angle{n}, half_length_pair2{n}, dist_pair_signal{n}, base_length{n}, direct_length{n}] = deal(nan(ntime, ntrial));
    for t = 1:ntime
        for tr = 1:ntrial
            cos_theta = dot(vector_pair1(:,t), vector_pair2(:,t,tr)) / norm(vector_pair1(:,t)) / norm(vector_pair2(:,t,tr));
            cos_theta = max(min(cos_theta, 1), -1);
            angle{n}(t,tr) = rad2deg(acos(cos_theta));

            half_length_pair2{n}(t,tr) = 0.5 * norm(vector_pair2(:,t,tr));
            dist_pair_signal{n}(t,tr) = norm(dist_cat1(:,t,tr));
            base_length{n}(t,tr) = sqrt(2 * half_length_pair2{n}(t,tr)^2 * (1 - cos(deg2rad(angle{n}(t,tr)))));
            direct_length{n}(t,tr) = sqrt(dist_pair_signal{n}(t,tr)^2 - base_length{n}(t,tr)^2);
        end
    end
end

save(fullfile(InterimDir, sprintf('summary_cat1_%s_%s.mat', monkey, experiment)), ...
    'angle', 'half_length_pair2', 'dist_pair_signal', 'base_length', 'direct_length', 'correct_target')

load(fullfile(InterimDir, sprintf('summary_cat1_%s_%s.mat', monkey, experiment)))

% plot average
angle_correct = cellfun(@(x,y) reshape(x(300:400,y(1,:)==1), [], 1), angle, correct_target, 'uni', 0);
angle_wrong = cellfun(@(x,y) reshape(x(300:400,y(1,:)==0), [], 1), angle, correct_target, 'uni', 0);
mn_correct = cellfun(@(x) mean(x), angle_correct);
mn_wrong = cellfun(@(x) mean(x), angle_wrong);
figure; hold on
plot(mn_correct, '.-')
plot(mn_wrong, '.-')
legend({'Correct', 'Wong'})
format_panel(gcf, 'xlabel', '#Session', 'ylabel', 'Angle')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_change_%s_%s.pdf', monkey, experiment)))

% plot distance
color_range = [0.1 0.4];
for n = 1:n_files
    figure('Position', [100 100 950 285]);
    ax = nan(2,2);

    ax(1,1) = subplot(2,5,1); hold on
    imagesc([angle{n}(:,correct_target{n}(1,:)==1), angle{n}(:,correct_target{n}(1,:)==0)]', [50 100])
    n_correct = sum(correct_target{n}(1,:)==1);
    plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
    title('Angle between cat axis of pair 1&2')

    ax(2,1) = subplot(2,5,1+5); hold on
    dist1 = angle{n}(300:400, correct_target{n}(1,:)==1);
    dist2 = angle{n}(300:400, correct_target{n}(1,:)==0);
    histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(50, 120, 40))
    histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(50, 120, 40))
    plot([90 90], ylim, 'r', 'LineWidth', 2)
    xlabel('Data'); ylabel('Probability')
    
    % ax(1,2) = subplot(2,5,2); hold on
    % imagesc([half_length_pair2{n}(:,correct_target{n}(1,:)==1), half_length_pair2{n}(:,correct_target{n}(1,:)==0)]', color_range)
    % plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
    % title('Dist between origin and pair2 cat1')
    % 
    % ax(2,2) = subplot(2,5,2+5); hold on
    % dist1 = half_length_pair2{n}(300:400, correct_target{n}(1,:)==1);
    % dist2 = half_length_pair2{n}(300:400, correct_target{n}(1,:)==0);
    % histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    % histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    
    ax(1,2) = subplot(2,5,2); hold on
    imagesc([dist_pair_signal{n}(:,correct_target{n}(1,:)==1), dist_pair_signal{n}(:,correct_target{n}(1,:)==0)]', color_range)
    plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
    title('Cat1 pair signal')

    ax(2,2) = subplot(2,5,2+5); hold on
    dist1 = dist_pair_signal{n}(300:400, correct_target{n}(1,:)==1);
    dist2 = dist_pair_signal{n}(300:400, correct_target{n}(1,:)==0);
    histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    title('Blue: Correct  Orange: Wrong')

    % ax(1,4) = subplot(2,5,4); hold on
    % imagesc([base_length{n}(:,correct_target{n}(1,:)==1), base_length{n}(:,correct_target{n}(1,:)==0)]', color_range)
    % plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
    % title('Angle-related pair signal')
    % 
    % ax(2,4) = subplot(2,5,4+5); hold on
    % dist1 = base_length{n}(300:400, correct_target{n}(1,:)==1);
    % dist2 = base_length{n}(300:400, correct_target{n}(1,:)==0);
    % histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    % histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    % 
    % ax(1,5) = subplot(2,5,5); hold on
    % imagesc([direct_length{n}(:,correct_target{n}(1,:)==1), direct_length{n}(:,correct_target{n}(1,:)==0)]', color_range)
    % plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
    % title('Angle-unrelated pair signal')
    % 
    % ax(2,5) = subplot(2,5,5+5); hold on
    % dist1 = direct_length{n}(300:400, correct_target{n}(1,:)==1);
    % dist2 = direct_length{n}(300:400, correct_target{n}(1,:)==0);
    % histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))
    % histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(0.1, 0.4, 40))

    format_panel(ax(1,:), 'xlim', [0 801], 'xlabel', 'Time (ms)', 'ylabel', '#Trial')
    format_panel(ax(2,:), 'lim_match', {0,0,0})
    print(gcf, '-dpdf', fullfile(FigDir, sprintf('summary_cat1_%s_%s_session%d.pdf', monkey, experiment, n)))
end




%% dissociate distance to angle-related and -unrelated
[angle, half_length_pair2, dist_pair_signal, base_length, direct_length] = deal(nan(801, n_files));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get distance
    dist_cat1 = psth(:,:,2) - psth(:,:,1);
    
    vector_pair1 = psth(:,:,3) - psth(:,:,1);
    vector_pair2 = psth(:,:,4) - psth(:,:,2);

    for t = 1:ntime
        cos_theta = dot(vector_pair1(:,t), vector_pair2(:,t)) / norm(vector_pair1(:,t)) / norm(vector_pair2(:,t));
        cos_theta = max(min(cos_theta, 1), -1);
        angle(t,n) = rad2deg(acos(cos_theta));

        half_length_pair2(t,n) = 0.5 * norm(vector_pair2(:,t));
        dist_pair_signal(t,n) = norm(dist_cat1(:,t));
        base_length(t,n) = sqrt(2 * half_length_pair2(t,n)^2 * (1 - cos(deg2rad(angle(t,n)))));
        direct_length(t,n) = sqrt(dist_pair_signal(t,n)^2 - base_length(t,n)^2);
    end
end

opt.plot = set_plot_opt('vik', n_files);
figure('Position', [50 100 1170 150]);
subplot(1,5,1); hold on
for n = 1:n_files; plot(angle(:,n), 'Color', opt.plot.color(n,:)); end

subplot(1,5,2); hold on
for n = 1:n_files; plot(half_length_pair2(:,n), 'Color', opt.plot.color(n,:)); end; title('Dist from origin to cat 1'); ylim([0.01 0.1])

subplot(1,5,3); hold on
for n = 1:n_files; plot(dist_pair_signal(:,n), 'Color', opt.plot.color(n,:)); end; title('Pair signal'); ylim([0.01 0.1])

subplot(1,5,4); hold on
for n = 1:n_files; plot(base_length(:,n), 'Color', opt.plot.color(n,:)); end; title('Angle related'); ylim([0.01 0.1])

subplot(1,5,5); hold on
for n = 1:n_files; plot(direct_length(:,n), 'Color', opt.plot.color(n,:)); end; title('Angle unrelated'); ylim([0.01 0.1])

%% check pair distance
% 单试次的pair distance可能不适合用来计算
[p1_c1_raster, p2_c1_raster, p1_c1_psth, p2_c1_psth, dist1, dist2, dist3] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get data
    p1_c1_raster{n} =  r_detrend(:,:,ID(1,:)==1);
    p2_c1_raster{n} =  r_detrend(:,:,ID(1,:)==2);
    p1_c1_psth{n} = mean(r_detrend(:,:,ID(1,:)==1), 3);
    p2_c1_psth{n} = mean(r_detrend(:,:,ID(1,:)==2), 3);

    % get distance
    p2_psth_to_p1_psth = p2_c1_psth{n} - p1_c1_psth{n};
    p2_raster_to_p1_psth = p2_c1_raster{n} - p1_c1_psth{n};
    p2_raster_to_p2_psth = p2_c1_raster{n} - p2_c1_psth{n};

    ntime = size(p2_raster_to_p1_psth, 2);
    ntrial = size(p2_raster_to_p1_psth, 3);
    dist1{n} = nan(ntime, 1); dist2{n} = nan(ntime, ntrial); dist3{n} = nan(ntime, ntrial);
    for t = 1:ntime
        dist1{n}(t,1) = norm(p2_psth_to_p1_psth(:,t));
        for tr = 1:ntrial
            dist2{n}(t,tr) = norm(p2_raster_to_p1_psth(:,t,tr));
            dist3{n}(t,tr) = norm(p2_raster_to_p2_psth(:,t,tr));
        end
    end
end


figure('Position', [100 100 840 788]);
for n = 1:n_files
    subplot(4,4,n); hold on
    histogram(p1_c1_raster{n}(:), 'Normalization', 'probability', 'BinEdges', linspace(-0.0981, 0.1853, 100))
    histogram(p2_c1_raster{n}(:), 'Normalization', 'probability', 'BinEdges', linspace(-0.0981, 0.1853, 100))
    histogram(p2_c1_psth{n}(:), 'Normalization', 'probability', 'BinEdges', linspace(-0.0981, 0.1853, 100))
    if n==1
        legend({'pair1 cat1 raster', 'pair2 cat2 raster', 'pair2 cat1 PSTH'});
        xlabel('Data'); ylabel('Probability')
    end
    title(sprintf('session %d', n))
end
print(gcf, '-dpdf', fullfile(FigDir, sprintf('check_raster_or_psth_data_%s_%s.pdf', monkey, experiment)))

figure('Position', [100 100 840 788]);
for n = 1:n_files
    subplot(4,4,n); hold on
    histogram(dist1{n}, 'Normalization', 'probability', 'BinEdges', linspace(0.0085, 0.4400, 100))
    histogram(dist2{n}(:), 'Normalization', 'probability', 'BinEdges', linspace(0.0085, 0.4400, 100))
    histogram(dist3{n}(:), 'Normalization', 'probability', 'BinEdges', linspace(0.0085, 0.4400, 100))
    if n==1
        legend({'p2_psth_to_p1_psth', 'p2_raster_to_p1_psth', 'p2_raster_to_p2_psth'}, 'Interpreter', 'none');
        xlabel('Data'); ylabel('Probability')
    end
    title(sprintf('session %d', n))
end
print(gcf, '-dpdf', fullfile(FigDir, sprintf('check_pair_distance_%s_%s.pdf', monkey, experiment)))

% whether distance changes with learning
avg_dist1 = cellfun(@(x) mean(x), dist1);
avg_dist2 = cellfun(@(x) mean(x(:)), dist2);
avg_dist3 = cellfun(@(x) mean(x(:)), dist3);
figure('Position', [100 100 560 420]);
subplot(1,3,1); plot(avg_dist1, '.-'); title('p2_psth_to_p1_psth', 'Interpreter', 'none');
subplot(1,3,2); plot(avg_dist2, '.-'); title('p2_raster_to_p1_psth', 'Interpreter', 'none')
subplot(1,3,3); plot(avg_dist3, '.-'); title('p2_raster_to_p2_psth', 'Interpreter', 'none')
format_panel(gcf, 'lim_match', {0,0,0}, 'xlabel', '#Session', 'ylabel', 'Distance')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('check_change_of_pair_distance_%s_%s.pdf', monkey, experiment)))



%% distance for estimating pair signal (single trial)
[dist_cat1, correct_target] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get single trial for test (both correct and wrong trials)
    r_target = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1);
    correct_target{n} = correct(:,task_set(1,:)==2 & targ_cor(1,:)==1);

    % get angle
    vector = r_target - psth(:,:,1);

    ntrial = size(vector, 3);
    dist_cat1{n} = nan(ntime, ntrial);
    for t = 1:ntime
        for tr = 1:ntrial
            dist_cat1{n}(t,tr) = norm(vector(:,t,tr));
        end
    end
end

save(fullfile(InterimDir, sprintf('dist_cat1_single_trial_%s_%s.mat', monkey, experiment)), ...
    'dist_cat1', 'correct_target')


% plot trial-averaged distance
figure; hold on
opt.plot = set_plot_opt('vik', n_files);
for n = 1:n_files
    avg = mean(dist_cat1{n}, 2);
    plot(avg, 'Color', opt.plot.color(n,:))
end

% plot distance
figure('Position', [100 100 870 770]);
for n = 1:n_files
    subplot(5,4,n); hold on
    imagesc([dist_cat1{n}(:,correct_target{n}(1,:)==1), dist_cat1{n}(:,correct_target{n}(1,:)==0)]', [0.1 0.4])
    colorbar

    n_correct = sum(correct_target{n}(1,:)==1);
    plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
end

% plot distribution at certain time range
figure('Position', [100 100 870 770]);
for n = 1:n_files
    subplot(5,4,n); hold on
    dist1 = dist_cat1{n}(300:400, correct_target{n}(1,:)==1);
    dist2 = dist_cat1{n}(300:400, correct_target{n}(1,:)==0);

    histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(0.15, 0.36, 40))
    histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(0.15, 0.36, 40))
end

%% distance for estimating pair signal
[dist_cat1, dist_cat2, dist_pair_signal] = deal(nan(801, n_files));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get angle
    vector_cat1 = psth(:,:,2) - psth(:,:,1);
    vector_cat2 = psth(:,:,4) - psth(:,:,3);

    for t = 1:ntime
        dist_cat1(t,n) = norm(vector_cat1(:,t));
        dist_cat2(t,n) = norm(vector_cat2(:,t));
        dist_pair_signal(t,n) = mean([dist_cat1(t,n) dist_cat2(t,n)]);
    end
end

opt.plot = set_plot_opt('vik', size(dist_cat1,2));
figure; hold on
for n = 1:n_files
    plot(dist_cat1(:,n), 'Color', opt.plot.color(n,:))
end

figure;
imagesc(dist_cat1')

%% angle of single trials
[angle, correct_target, morph_target, fano_correct, fano_wrong, fano_all] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    correct = double(targ_cho==targ_cor);
    ID = task_set;
    ID(targ_cor==2) = ID(targ_cor==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'correct'}, {task_set, targ_cor, correct})

    % get trial-averaged baseline (correct trials only)
    [nunit, ntime, ~] = size(r);
    ncond = max(ID(:));

    psth = nan(nunit, ntime, ncond);
    for c = 1:ncond
        psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
    end

    % get single trial for test (both correct and wrong trials)
    r_target = r_detrend(:,:,task_set(1,:)==2 & targ_cor(1,:)==1);
    correct_target{n} = correct(:,task_set(1,:)==2 & targ_cor(1,:)==1);
    morph_target{n} = morph_level(:,task_set(1,:)==2 & targ_cor(1,:)==1);
    r_fano = r(:,:,task_set(1,:)==2 & targ_cor(1,:)==1);

    % get angle
    vector_pair1 = psth(:,:,3) - psth(:,:,1);
    vector_pair2 = psth(:,:,4) - r_target;

    ntrial = size(vector_pair2, 3);
    angle{n} = nan(ntime, ntrial);
    for t = 1:ntime
        for tr = 1:ntrial
            cos_theta = dot(vector_pair1(:,t), vector_pair2(:,t,tr)) / norm(vector_pair1(:,t)) / norm(vector_pair2(:,t,tr));
            cos_theta = max(min(cos_theta, 1), -1);
            angle{n}(t,tr) = rad2deg(acos(cos_theta));
        end
    end

    % get Fano factor
    r_correct = sum(r_fano(:,:,correct_target{n}(1,:)==1), 2); mn_correct = mean(r_correct, 3); var_correct = var(r_correct, [], 3);
    fano_correct{n} = var_correct ./ mn_correct;
    r_wrong = sum(r_fano(:,:,correct_target{n}(1,:)==0), 2); mn_wrong = mean(r_wrong, 3); var_wrong= var(r_wrong, [], 3);
    fano_wrong{n} = var_wrong ./ mn_wrong;

    r_all = sum(r_fano, 2); mn_correct = mean(r_correct, 3); var_correct = var(r_correct, [], 3);
    fano_all{n} = var_correct ./ mn_correct;
end

save(fullfile(InterimDir, sprintf('angle_single_trial_%s_%s.mat', monkey, experiment)), ...
    'angle', 'correct_target', 'morph_target', 'fano_correct', 'fano_wrong', 'fano_all')

% relation between angle and coh
[p_error, pse_error, coh_error, p_correct, pse_correct, coh_correct] = deal(cell(n_files, 1));
for n = 1:n_files
    I = correct_target{n}(1,:)==0;
    coh_error{n} = unique(morph_target{n}(1,I)); ncoh = length(coh_error{n}); ntime = size(angle{n},1);
    [p_error{n}, pse_error{n}] = deal(nan(ntime, ncoh));
    for t = 1:ntime
        [p_error{n}(t,:), pse_error{n}(t,:)] = calcGroupMean(mean(angle{n}(t,I), 1), morph_target{n}(1,I), unique(morph_target{n}(1,I)), 'continuous');
    end

    I = correct_target{n}(1,:)==1;
    coh_correct{n} = unique(morph_target{n}(1,I)); ncoh = length(coh_correct{n}); ntime = size(angle{n},1);
    [p_correct{n}, pse_correct{n}] = deal(nan(ntime, ncoh));
    for t = 1:ntime
        [p_correct{n}(t,:), pse_correct{n}(t,:)] = calcGroupMean(mean(angle{n}(t,I), 1), morph_target{n}(1,I), unique(morph_target{n}(1,I)), 'continuous');
    end
end

% plot time sequence
coh_list = coh_correct{end};
opt.plot = set_plot_opt('bamako', length(coh_list));
figure('Position', [100 100 800 650]);
for n = 1:n_files
    subplot(4,4,n); hold on
    for c = 1:length(coh_list)
        if any(coh_correct{n}==coh_list(c))
            plot(-100:700, p_correct{n}(:, coh_correct{n}==coh_list(c)), 'Color', opt.plot.color(c,:))
        end
    end
end
format_panel(gcf, 'ylim', [70 95], 'xlim', [-100 700], 'xlabel', 'Time (ms)', 'ylabel', 'Angle')

figure('Position', [100 100 800 650]);
for n = 1:n_files
    subplot(4,4,n); hold on
    for c = 1:length(coh_list)
        if any(coh_error{n}==coh_list(c))
            plot(-100:700, p_error{n}(:, coh_error{n}==coh_list(c)), 'Color', opt.plot.color(c,:))
        end
    end
end
format_panel(gcf, 'ylim', [70 95], 'xlim', [-100 700], 'xlabel', 'Time (ms)', 'ylabel', 'Angle')


% plot average across time
opt.color = [1 0 0; 0 0 0];
figure('Position', [100 100 800 650]);
for n = 1:n_files
    subplot(4,4,n); hold on
    plot(coh_error{n}, mean(p_error{n}(300:400,:), 1), '.-', 'markers', 7, 'Color', opt.color(1,:));
    cerrorbar(coh_error{n}, mean(p_error{n}(300:400,:), 1), mean(pse_error{n}(300:400,:), 1), 'Color', opt.color(1,:));
    plot(coh_correct{n}, mean(p_correct{n}(300:400,:), 1), '.-', 'markers', 7, 'Color', opt.color(2,:));
    cerrorbar(coh_correct{n}, mean(p_correct{n}(300:400,:), 1), mean(pse_correct{n}(300:400,:), 1), 'Color', opt.color(2,:));
    plot([0 1], [90 90], ':', 'Color', 'black')
end
format_panel(gcf, 'ylim', [70 95], 'xlabel', 'Morph', 'ylabel', 'Angle')


% plot Fano factor of both correct and wrong trials
mn = cellfun(@(x) nanmean(x), fano_all);
figure('Position', [100 100 250 200]); hold on
plot(mn, '.-')
xlabel('#Session'); ylabel('Fano Factor'); title(sprintf('%s %s', monkey, experiment))

% plot mean of Fano factor
mn_correct = cellfun(@(x) nanmean(x), fano_correct);
mn_wrong = cellfun(@(x) nanmean(x), fano_wrong);
figure('Position', [100 100 250 200]); hold on
plot(mn_correct, '.-')
plot(mn_wrong, '.-')
legend({'Correct', 'Wrong'}, 'Location', 'best')
xlabel('#Session'); ylabel('Fano Factor'); title(sprintf('%s %s', monkey, experiment))
print(gcf, '-dpdf', fullfile(FigDir, sprintf('fano_factor_%s_%s.pdf', monkey, experiment)))

% plot distribution of Fano factor
figure('Position', [100 100 870 770]);
for n = 1:n_files
    subplot(5,4,n); hold on
    histogram(fano_correct{n}, 'Normalization', 'probability', 'BinEdges', linspace(0.22, 6.3, 20))
    histogram(fano_wrong{n}, 'Normalization', 'probability', 'BinEdges', linspace(0.22, 6.3, 20))
end

% plot angle
figure('Position', [100 100 870 770]);
for n = 1:n_files
    subplot(5,4,n); hold on
    imagesc([angle{n}(:,correct_target{n}(1,:)==1), angle{n}(:,correct_target{n}(1,:)==0)]', [50 100])
    colorbar

    n_correct = sum(correct_target{n}(1,:)==1);
    plot([0 801], [n_correct n_correct], 'r', 'LineWidth', 2)
end
% print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle_single_trial_%s_%s.pdf', monkey, experiment)))

% plot distribution at certain time range
figure('Position', [100 100 870 770]);
for n = 1:n_files
    subplot(5,4,n); hold on
    dist1 = angle{n}(300:400, correct_target{n}(1,:)==1);
    dist2 = angle{n}(300:400, correct_target{n}(1,:)==0);

    histogram(dist1(:), 'Normalization', 'probability', 'BinEdges', linspace(50, 120, 40))
    histogram(dist2(:), 'Normalization', 'probability', 'BinEdges', linspace(50, 120, 40))
    plot([90 90], ylim, 'r', 'LineWidth', 2)
end



%% vector angle
angle = nan(801, n_files);
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor')); % correct trials only

    task_set = fnd.getp('task_set');
    targ_cho = fnd.getp('targ_cho');
    ID = task_set;
    ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:));
    psth = fnd.PSTH(ID, {'boxcar', 100}, [], true); psth = psth{1}; % (unit, time, condition)

    [nunit, ntime, ncond] = size(psth);
    vector_pair1 = psth(:,:,3) - psth(:,:,1);
    vector_pair2 = psth(:,:,4) - psth(:,:,2);

    for t = 1:ntime
        cos_theta = dot(vector_pair1(:,t), vector_pair2(:,t)) / norm(vector_pair1(:,t)) / norm(vector_pair2(:,t));
        cos_theta = max(min(cos_theta, 1), -1);
        angle(t,n) = rad2deg(acos(cos_theta));
    end
end
save(fullfile(InterimDir, sprintf('angle_%s_%s.mat', monkey, experiment)), 'angle')

%% plot angle
monkey_list = {'Nick', 'Nick', 'Nick', 'Woody', 'Woody', 'Nick'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4', 'learnTask3', 'learnTask4', 'faceColor'};

figure('Position', [100 100 560 820]);
tstamp = -100:700;
[ax_time, ax_session] = deal(nan(length(monkey_list), 1));
for n = 1:length(monkey_list)
    monkey = monkey_list{n};
    experiment = experiment_list{n};
    angle = load(fullfile(InterimDir, sprintf('angle_%s_%s.mat', monkey, experiment))).angle;

    ax_time(n) = subplot(length(monkey_list), 2, 2*n-1); hold on
    opt.plot = set_plot_opt('vik', size(angle,2));
    for m = 1:size(angle,2)
        plot(tstamp, angle(:,m), 'Color', opt.plot.color(m,:))
    end
    title(sprintf('%s %s', monkey, experiment))

    ax_session(n) = subplot(length(monkey_list), 2, 2*n);
    plot(angle(tstamp==400, :), '.-')
    if n==1; title('Angle at 400 ms'); end
end
format_panel(ax_time, 'xlabel', 'Time (ms)', 'ylabel', 'Angle (deg)', 'xlim', [tstamp(1) tstamp(end)], 'ylim', [20 110], 'axis', 'normal')
format_panel(ax_session, 'xlabel', '#Session', 'ylabel', 'Angle (deg)', 'ylim', [20 60], 'axis', 'normal')

print(gcf, '-dpdf', fullfile(FigDir, sprintf('angle.pdf')))


%% Euclidean distance without choice signal (for each choice)
EPOCH = [2 3];
for n = [1 n_files]
    stat = cell(1, length(EPOCH));
    for e = 1:length(EPOCH)
        % load data
        fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
        fnd = fnd.extract_epoch(EPOCH(e));

        % select unit
        r = fnd.FR({1, [100 500]});
        I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);

        % select trial
        % fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor')); % correct trials only

        % get ID
        task_set = fnd.getp('task_set');
        targ_cho = fnd.getp('targ_cho');
        targ_cor = fnd.getp('targ_cor');

        ID = task_set;
        ID(targ_cho==2) = ID(targ_cho==2)+max(ID(:));
        trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set});

        % process PSTH
        r = fnd.raster; r = r{1}*1000; % (unit, time, trial)
        for c = 1:2 % c for choice
            mn = nanmean(r(:,:,targ_cho(1,:)==c), 3);
            % sd = nanstd(r(:,:,targ_cho(1,:)==c), [], 3); sd(sd==0) = 0.001;
            r(:,:,targ_cho(1,:)==c) = (r(:,:,targ_cho(1,:)==c) - mn);
        end

        PSTH_unsmooth = nan(size(r,1), size(r,2), max(ID(:)));
        for c = 1:max(ID(:))
            PSTH_unsmooth(:,:,c) = nanmean(r(:,:,ID(1,:)==c), 3);
        end
        kernel = fspecial('average', [1, 100]);
        psth = nanconv(PSTH_unsmooth, kernel, 'same');

        % get distance
        stat{e}.tstamp = fnd.tstamp{1};
        [nunit, ntime, ntrial] = size(fnd.data{1});

        dist_choice1 = nan(ntime, 1);
        dist_choice2 = nan(ntime, 1);
        for t = 1:ntime
            dist_choice1(t,1) = pdist([psth(:,t,1)'; psth(:,t,2)']);
            dist_choice2(t,1) = pdist([psth(:,t,3)'; psth(:,t,4)']);
        end
        stat{e}.dist_choice1 = dist_choice1;
        stat{e}.dist_choice2 = dist_choice2;
    end

    % plot
    fh = figure('Position', [50 100 450 200]);
    for e = 1:length(EPOCH)
        subplot(1,2,e); hold on
        plot(stat{e}.tstamp, stat{e}.dist_choice1)
        plot(stat{e}.tstamp, stat{e}.dist_choice2)
        if e==1; title(sprintf('Session %d', n)); xlabel('Time from stim on (ms)'); end
        if e==2; legend({'Choice 1', 'Choice 2'}, 'Location', 'best'); xlabel('Time from resp (ms)'); end
        format_panel(gca, 'ylabel', 'Distance', 'ylim', [20 140])
    end
    print(fh, '-dpdf', fullfile(FigDir, sprintf('distance_%s_%s_session%d.pdf', monkey, experiment, n)))
end

%% Euclidean distance (for each choice)
EPOCH = [2 3];
for n = [1 n_files]
    stat = cell(1, length(EPOCH));
    for e = 1:length(EPOCH)
        % load data
        fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
        fnd = fnd.extract_epoch(EPOCH(e));

        % select unit
        r = fnd.FR({1, [100 500]});
        I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
        fnd = fnd.set_unit_criteria('custom', I);

        % select trial
        fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor')); % correct trials only

        % get ID
        task_set = fnd.getp('task_set');
        targ_cho = fnd.getp('targ_cho');
        targ_cor = fnd.getp('targ_cor');

        ID = task_set;
        ID(targ_cho==2) = ID(targ_cho==2)+max(ID(:));
        trial_classifier_result(ID, {'targ_cor', 'targ_cho', 'task_set'}, {targ_cor, targ_cho, task_set});

        % get distance
        stat{e}.tstamp = fnd.tstamp{1};
        [nunit, ntime, ntrial] = size(fnd.data{1});
        psth = fnd.PSTH(ID, {'boxcar', 100}); psth = psth{1}; % (unit, time, condition)

        dist_choice1 = nan(ntime, 1);
        dist_choice2 = nan(ntime, 1);
        for t = 1:ntime
            dist_choice1(t,1) = pdist([psth(:,t,1)'; psth(:,t,2)']);
            dist_choice2(t,1) = pdist([psth(:,t,3)'; psth(:,t,4)']);
        end
        stat{e}.dist_choice1 = dist_choice1;
        stat{e}.dist_choice2 = dist_choice2;
    end

    % plot
    fh = figure('Position', [50 100 450 200]);
    for e = 1:length(EPOCH)
        subplot(1,2,e); hold on
        plot(stat{e}.tstamp, stat{e}.dist_choice1)
        plot(stat{e}.tstamp, stat{e}.dist_choice2)
        if e==1; title(sprintf('Session %d', n)); xlabel('Time from stim on (ms)'); end
        if e==2; legend({'Choice 1', 'Choice 2'}, 'Location', 'best'); xlabel('Time from resp (ms)'); end
        format_panel(gca, 'ylabel', 'Distance', 'ylim', [20 140])
    end
    print(fh, '-dpdf', fullfile(FigDir, sprintf('distance_%s_%s_session%d.pdf', monkey, experiment, n)))
end




