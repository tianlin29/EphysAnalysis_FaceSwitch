run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'faceColor'; % learnTask2, learnTask3, learnTask4, faceColor, faceColor_passiveLong, threeExemplar
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'PCA'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'PCA'); mkdir(InterimDir);

%%
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);
fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor')); % correct trials only

task_set = fnd.getp('task_set');
targ_cho = fnd.getp('targ_cho');
ID = task_set;
ID(targ_cho==2) = ID(targ_cho==2) + max(ID(:));
psth = fnd.PSTH(ID, {'boxcar', 100}, [], true); psth = psth{1}; % (unit, time, condition)

% tstamp = fnd.tstamp{1};
% psth = psth(:, tstamp>=0 & tstamp<=400, :);

[nunit, ntime, ncond] = size(psth);

[coeff, score, latent] = pca(psth(:,:)');
% coeff ..(unit, pc) 主成分系数，每一列是原空间中的一个PC轴，可用于投影
% score ..(observation, pc) 主成分分数，每一行是一个原始数据点在PC空间中的坐标，也就是降维后的数据
% latent ..(pc, 1) 主成分方差，每一行是一个主成分的方差

psth_pc = reshape(score', nunit, ntime, ncond);

% check PCA result
opt.plot = set_plot_opt_2cond('roma', 'roma', 2);
figure; hold on
for c = 1:ncond
    plot3(psth_pc(1,:,c), psth_pc(2,:,c), psth_pc(3,:,c), 'Color', opt.plot.color(c,:), 'LineStyle', opt.plot.linestyle{c})
end
format_panel(gcf, 'xlabel', 'PC 1', 'ylabel', 'PC 2', 'zlabel', 'PC 3')
view([45 45])

% check variance explained
latent_cumsum = cumsum(latent/sum(latent));

figure;
plot(latent_cumsum, '.-')
format_panel(gcf, 'xlabel', '#PC', 'ylabel', 'Explained variance', 'ylim', [0 1])

% project pair 1 to the pc space
psth_pair1 = psth(:,:,[1 2]);
score_pair1 = psth_pair1(:,:)' * coeff;
latent_pair1 = var(score_pair1)';
explained_pair1 = latent_pair1 / sum(latent_pair1);

% project pair 2 to the pc space
psth_pair2 = psth(:,:,[3 4]);
score_pair2 = psth_pair2(:,:)' * coeff;
latent_pair2 = var(score_pair2)';
explained_pair2 = latent_pair2 / sum(latent_pair2);

figure;
bar([explained_pair1(1:6), explained_pair2(1:6)])
legend({'Pair 1', 'New pair'})
format_panel(gcf, 'xlabel', '#PC', 'ylabel', 'Explained variance', 'ylim', [0 1])

%%
psth_pair1 = psth(:,:,[1 2]);
psth_pair2 = psth(:,:,[3 4]);

[coeff, score, latent] = pca(psth_pair1(:,:)'); latent_cumsum_pair1 = cumsum(latent/sum(latent));
[~, ~, latent] = pca(psth_pair2(:,:)'); latent_cumsum_pair2 = cumsum(latent/sum(latent));

% coeff ..(unit, pc) 主成分系数，每一列是原空间中的一个PC轴，可用于投影
% score ..(observation, pc) 主成分分数，每一行是一个原始数据点在PC空间中的坐标，也就是降维后的数据
% latent ..(pc, 1) 主成分方差，每一行是一个主成分的方差

psth_pair2 = psth_pair2 - mean(psth_pair1(:,:)')';
score_pair2 = psth_pair2(:,:)' * coeff;
latent_pair2 = var(score_pair2)';
latent_cumsum_cross = cumsum(latent_pair2/sum(latent_pair2));

figure; hold on
plot(latent_cumsum_pair1(1:6), '.-')
plot(latent_cumsum_pair2(1:6), '.-')
plot(latent_cumsum_cross(1:6), '.-')
format_panel(gcf, 'xlabel', '#PC', 'ylabel', 'Explained variance', 'ylim', [0 1])




