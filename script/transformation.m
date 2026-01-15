run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'faceColor'; % learnTask2, learnTask3, learnTask4, faceColor, faceColor_passiveLong, threeExemplar
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'distance'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'distance'); mkdir(InterimDir);

%% Procrustes transformation (trial-averaged)
[theta, scale, abs_translation, psth_register, psth_new, explained_var] = deal(cell(n_files, 1));
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

    % get task variables
    morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:);
    task_set = fnd.getp('task_set'); task_set = task_set(1,:);
    targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:);
    targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
    correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor')); correct = correct(1,:);

    % get ID
    % ucoh = unique(morph_level(task_set==1));
    % ntrial = arrayfun(@(x) sum(morph_level(task_set==1)==x), ucoh);
    % ucoh_pair1 = ucoh(ntrial>=10);
    % 
    % ucoh = unique(morph_level(task_set==2));
    % ntrial = arrayfun(@(x) sum(morph_level(task_set==2)==x), ucoh);
    % ucoh_pair2 = ucoh(ntrial>=10);
    % 
    % ucoh = intersect(ucoh_pair1, ucoh_pair2);
    % 
    % [~, ID] = ismember(morph_level, ucoh); ID(ID==0) = NaN;
    % ID(task_set==2) = ID(task_set==2) + max(ID(:));
    % ID(~correct) = NaN;
    % trial_classifier_result(ID, {'task_set', 'targ_cor', 'morph_level', 'correct'}, {task_set, targ_cor, morph_level, correct})

    % get ID
    ID = targ_cor;
    ID(task_set==2) = ID(task_set==2) + max(ID(:));
    ID(~correct) = NaN;
    trial_classifier_result(ID, {'task_set', 'targ_cor', 'morph_level', 'correct'}, {task_set, targ_cor, morph_level, correct})

    % get ID
    % ucoh = unique(morph_level(task_set==1));
    % ntrial = arrayfun(@(x) sum(morph_level(task_set==1)==x), ucoh);
    % ucoh_pair1 = ucoh(ntrial>=10);
    % 
    % ucoh = unique(morph_level(task_set==2));
    % ntrial = arrayfun(@(x) sum(morph_level(task_set==2)==x), ucoh);
    % ucoh_pair2 = ucoh(ntrial>=10);
    % 
    % ucoh = intersect(ucoh_pair1, ucoh_pair2);
    % ncoh = length(ucoh); split_points = round([0, ncoh/3, 2*ncoh/3, ncoh]);
    % segments = arrayfun(@(i) ucoh_pair2(split_points(i)+1:split_points(i+1)), 1:3, 'UniformOutput', false);
    % [seg1, seg2, seg3] = deal(segments{:});
    % 
    % ID = nan(size(morph_level));
    % for c = 1:3
    %     ID(ismember(morph_level, segments{c})) = c;
    % end
    % ID(task_set==2) = ID(task_set==2) + max(ID(:));
    % ID(~correct) = NaN;
    % trial_classifier_result(ID, {'task_set', 'targ_cor', 'morph_level', 'correct'}, {task_set, targ_cor, morph_level, correct})

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

        explained_var{n} = cumsum(latent/sum(latent));
    end
    
    % get single-trial transformation
    tstamp = fnd.tstamp{1};
    t_list = [linspace(min(tstamp)+50, max(tstamp)-50, 20)' - 50, linspace(min(tstamp)+50, max(tstamp)-50, 20)' + 50];
    ntime = size(t_list,1);
    [scale{n}, theta{n}, abs_translation{n}] = deal(nan(ntime, 1));
    psth_register{n} = nan(npc, ntime, ncond);
    psth_new{n} = nan(npc, ntime, ncond/2);
    for t = 1:ntime
        I = tstamp>=t_list(t,1) & tstamp<=t_list(t,2);
        psth_register{n}(:,t,:) = mean(psth_pc(:,I,:), 2);
        points_pair1 = permute(mean(psth_pc(:,I,1:ncond/2), 2), [3 1 2]); % (condition, unit) -> (point, dim)
        points_pair2 = permute(mean(psth_pc(:,I,ncond/2+1:ncond), 2), [3 1 2]);

        % R_est = findMinAngleRotation(points_pair1'-mean(points_pair1)', points_pair2'-mean(points_pair2)'); % (dim, point)
        % theta{n}(t) = acos((trace(R_est) - 1)/2);

        [~, Z, transform] = procrustes(points_pair1, points_pair2, 'scaling', true, 'reflection', false); % (condition, unit) -> (point, dim)
        if npc==2; theta{n}(t) = atan2d(transform.T(2,1),transform.T(1,1)); end
        if npc==3; theta{n}(t) = acosd((trace(transform.T)-1)/2); end
        if npc>3
            v_unit = ones(npc,1);
            v_trans = transform.T * v_unit;
            cos_theta = dot(v_unit, v_trans) / norm(v_unit) / norm(v_trans);
            cos_theta = max(min(cos_theta, 1), -1);
            theta{n}(t) = rad2deg(acos(cos_theta));
        end
        scale{n}(t) = transform.b;
        abs_translation{n}(t) = norm(transform.c(1,:));
        psth_new{n}(:,t,:) = Z';
    end
end

save(fullfile(InterimDir, sprintf('Procrustes_average_%s_%s.mat', monkey, experiment)), ...
    'theta', 'scale', 'abs_translation', 'psth_register', 'psth_new', 'explained_var')

% plot pc
figure('Position', [100 100 200 200]); hold on
opt.plot = set_plot_opt('vik', n_files);
for n = 1:n_files
    plot(explained_var{n}, 'Color', opt.plot.color(n,:))
end
plot(xlim, [0.9 0.9], ':', 'Color', 'black')
format_panel(gcf, 'xlabel', '#PC', 'ylabel', 'Explained variance', 'ylim', [0 1], 'xtick', [1 50 100 150 200 250])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_PC_%s_%s.pdf', monkey, experiment)))

% plot imagesc
theta_all = abs(cell2mat(theta'));
scale_all = cell2mat(scale');
abs_translation_all = cell2mat(abs_translation');

figure('Position', [100 100 740 180]);
subplot(1,3,1); hold on; title('Theta')
imagesc(theta_all', [0 90]); colorbar
set(gca, 'YDir', 'reverse');

subplot(1,3,2); hold on; title('Scale')
imagesc(scale_all'); colorbar
set(gca, 'YDir', 'reverse');

subplot(1,3,3); hold on; title('Translation')
imagesc(abs_translation_all'); colorbar
set(gca, 'YDir', 'reverse');

format_panel(gcf, 'lim_match', {0,0,0}, 'xlabel', 'Time', 'ylabel', '#Session')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_heatmap_%s_%s.pdf', monkey, experiment)))


% plot change of time
opt.plot = set_plot_opt('vik', n_files);
figure('Position', [100 100 740 180]);
subplot(1,3,1); hold on; title('Theta')
for n = 1:n_files
    plot(theta_all(:,n), 'Color', opt.plot.color(n,:))
end
subplot(1,3,2); hold on; title('Scale')
for n = 1:n_files
    plot(scale_all(:,n), 'Color', opt.plot.color(n,:))
end
subplot(1,3,3); hold on; title('Translation')
for n = 1:n_files
    plot(abs_translation_all(:,n), 'Color', opt.plot.color(n,:))
end
format_panel(gcf, 'lim_match', {1,0,0}, 'xlabel', 'Time')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_time_%s_%s.pdf', monkey, experiment)))

% plot change of session
t = 10;
figure('Position', [100 100 740 180]);
subplot(1,3,1)
plot(theta_all(t,:))
title('Theta')
subplot(1,3,2)
plot(scale_all(t,:))
title('Scale')
subplot(1,3,3)
plot(abs_translation_all(t,:))
title('Translation')
format_panel(gcf, 'lim_match', {1,0,0}, 'xlabel', '#Session')
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_session_%s_%s.pdf', monkey, experiment)))

% plot alignment (2d)
if npc==2
    figure('Position', [100 100 800 950]);
    for n = 1:n_files
        ncond = size(psth_register{n},3);
        opt.plot = set_plot_opt_2cond('vik', 'vik', ncond/2);
        subplot(5,5,n); hold on; title(sprintf('session %d', n))
        for c = 1:ncond
            plot(psth_register{n}(1,:,c), psth_register{n}(2,:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
        end
        for c = 1:ncond/2
            plot(psth_new{n}(1,:,c), psth_new{n}(2,:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{end}, 'LineWidth', 2)
        end
    end
end
format_panel(gcf, 'lim_match', {1,1,1}, 'xlabel', 'PC 1', 'ylabel', 'PC 2', 'xlim', [-0.06 0.06], 'ylim', [-0.06 0.06])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('Procrustes_alignment_%s_%s.pdf', monkey, experiment)))

% plot alignment (3d)
t = 4;
if npc==3
    figure('Position', [100 100 800 950]);
    for n = n_files
        ncond = size(psth_register{n},3);
        opt.plot = set_plot_opt_2cond('roma', 'roma', ncond/2);
        % subplot(5,5,n);
        hold on; title(sprintf('session %d', n))
        for c = 1:ncond
            plot3(psth_register{n}(1,:,c), psth_register{n}(2,:,c), psth_register{n}(3,:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
            scatter3(psth_register{n}(1,t,c), psth_register{n}(2,t,c), psth_register{n}(3,t,c), 40, opt.plot.color(c,:), 'filled')
        end
        for c = 1:ncond/2
            plot3(psth_new{n}(1,:,c), psth_new{n}(2,:,c), psth_new{n}(3,:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{end}, 'LineWidth', 2)
            scatter3(psth_new{n}(1,t,c), psth_new{n}(2,t,c), psth_new{n}(3,t,c), 40, opt.plot.color(c,:), 'filled', 'MarkerEdgeColor', 'red', 'LineWidth', 2)
        end
    end
    view([45 45]); grid on; axis equal
end
format_panel(gcf, 'lim_match', {1,1,1}, 'xlabel', 'PC 1', 'ylabel', 'PC 2', 'xlim', [-0.06 0.06], 'ylim', [-0.06 0.06])

%% functions
function R = findMinAngleRotation(A, B)
% 寻找旋转角度最小的旋转矩阵
% 输入：
%   A, B: n×m 矩阵，每组m个n维点（m < n）
% 输出：
%   R: n×n 旋转矩阵（正交，det=1），使得 ||R*A - B||_F 最小

% 确保输入正确
[n, m] = size(A);
assert(all(size(A) == size(B)), 'A和B维度必须相同');
assert(m < n, 'm必须小于n');

% 中心化点云
A_mean = mean(A, 2);
B_mean = mean(B, 2);
A_centered = A - A_mean;
B_centered = B - B_mean;

% 计算协方差矩阵
H = A_centered * B_centered';

% SVD分解
[U, ~, V] = svd(H);

% 计算初始旋转矩阵
R_init = V * U';

% 确保是旋转矩阵（det=1）
if det(R_init) < 0
    % 处理反射情况
    V(:, end) = -V(:, end);
    R_init = V * U';
end

% 现在需要找到所有可能的旋转矩阵并选择角度最小的
% 由于m<n，解空间有自由度

% 方法1：使用四元数表示法寻找最小角度解
if n == 3
    % 对于3D情况，可以直接计算
    R = findMinAngleRotation3D(A, B);
else
    % 对于n维情况，使用优化方法
    R = optimizeMinAngleRotation(A, B, R_init);
end
end

function R = findMinAngleRotation3D(A, B)
% 3D特殊情况
[U, ~, V] = svd(A * B');

% 确保是旋转矩阵
R = V * U';
if det(R) < 0
    % 处理反射
    V(:, 3) = -V(:, 3);
    R = V * U';
end

% 计算旋转角度
theta = acos((trace(R) - 1)/2);

% 如果需要寻找更小角度的解，可以尝试添加2π的倍数
% 但实际上SVD已经给出了最小Frobenius范数误差的解
end

function R_opt = optimizeMinAngleRotation(A, B, R_init)
% 使用优化寻找最小角度旋转矩阵

% 参数：使用旋转矩阵的李代数表示
n = size(A, 1);

% 定义目标函数：最小化旋转角度 + 对齐误差
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');

% 使用李代数参数化SO(n)
if n == 3
    % 3D情况：使用轴角表示
    x0 = [0; 0; 0];  % 初始角度为0
else
    % n维情况：使用反对称矩阵参数化
    x0 = zeros(n*(n-1)/2, 1);  % 李代数参数
end

% 优化
x_opt = fmincon(@(x)rotationObjective(x, A, B, R_init), x0, ...
                [], [], [], [], [], [], [], options);
            
% 从优化参数重建旋转矩阵
R_opt = reconstructRotation(x_opt, n);
end

function cost = rotationObjective(x, A, B, R_init)
% 目标函数：旋转角度 + 对齐误差权重
n = size(A, 1);
R = reconstructRotation(x, n);

% 计算旋转角度（通过矩阵指数）
if n == 3
    % 3D：直接计算角度
    theta = acos(min(max((trace(R) - 1)/2, -1), 1));
else
    % n维：计算对数得到角度
    [U, D, V] = svd(R);
    S = U * diag(sign(diag(D))) * V';
    R_corrected = U * S * V';
    logR = logm(R_corrected);
    theta = norm(logR, 'fro') / sqrt(2);
end

% 计算对齐误差
alignment_error = norm(R * A - B, 'fro');

% 总代价：主要最小化角度，加上小权重的对齐误差
cost = theta + 0.01 * alignment_error;
end

function R = reconstructRotation(x, n)
% 从参数重建旋转矩阵
if n == 3
    % 3D：轴角表示
    theta = norm(x);
    if theta < 1e-10
        R = eye(3);
    else
        axis = x / theta;
        K = [0, -axis(3), axis(2);
             axis(3), 0, -axis(1);
             -axis(2), axis(1), 0];
        R = eye(3) + sin(theta) * K + (1 - cos(theta)) * K^2;
    end
else
    % n维：反对称矩阵指数
    % 重建反对称矩阵
    idx = 1;
    Omega = zeros(n);
    for i = 1:n
        for j = i+1:n
            Omega(i, j) = x(idx);
            Omega(j, i) = -x(idx);
            idx = idx + 1;
        end
    end
    
    % 矩阵指数得到旋转矩阵
    R = expm(Omega);
end
end

function [R, theta_min] = svdRotationWithMinAngle(A, B)
% 另一种方法：使用SVD并选择最小角度解
[n, ~] = size(A);

% 中心化
A_mean = mean(A, 2);
B_mean = mean(B, 2);
A_c = A - A_mean;
B_c = B - B_mean;

% SVD
[U, ~, V] = svd(A_c * B_c');

% 基本解
R_base = V * U';

% 确保是旋转矩阵
if det(R_base) < 0
    U(:, n) = -U(:, n);
    R_base = V * U';
end

% 对于m<n的情况，存在自由度
% 我们可以添加任意旋转到零空间部分
[U_A, S_A, V_A] = svd(A_c);
[U_B, S_B, V_B] = svd(B_c);

% 非零奇异值对应的子空间
r = rank(A_c);
U_A_r = U_A(:, 1:r);
U_B_r = U_B(:, 1:r);

% 零空间
if r < n
    U_A_null = U_A(:, r+1:end);
    U_B_null = U_B(:, r+1:end);
    
    % 在零空间中寻找最小角度旋转
    % 零空间中的任意正交变换不会影响对齐
    R_null = U_B_null * U_A_null';  % 将零空间对齐
    
    % 组合旋转矩阵
    R_full = U_B_r * U_A_r' + U_B_null * U_A_null';
    
    % 确保正交性（可能需要重新正交化）
    [U_full, ~, V_full] = svd(R_full);
    R = U_full * V_full';
    
    % 计算角度
    theta = acos(min(max((trace(R) - 1)/2, -1), 1));
    
    % 选择R_base和R中角度较小的
    theta_base = acos(min(max((trace(R_base) - 1)/2, -1), 1));
    
    if theta < theta_base
        R_opt = R;
        theta_min = theta;
    else
        R_opt = R_base;
        theta_min = theta_base;
    end
else
    R_opt = R_base;
    theta_min = acos(min(max((trace(R_base) - 1)/2, -1), 1));
end

R = R_opt;
end




