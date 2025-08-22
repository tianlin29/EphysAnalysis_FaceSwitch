clear;


%% PCA：主成分系数、分数和方差
load hald
size(ingredients); % (13 observation, 4 feature) = (cond*time, unit)

[coeff, score, latent] = pca(ingredients);
% coeff ..(4 feature, 4 pc) 主成分系数，每一列是原空间中的一个PC轴，可用于投影
% score ..(13 observation, 4 pc) 主成分分数，每一行是一个原始数据点在PC空间中的坐标，也就是降维后的数据
% latent ..(4 pc, 1) 主成分方差，每一行是一个主成分的方差

k = 3; % 选择前三个主成分
score_reduced = score(:, 1:k);
coeff_reduced = coeff(:, 1:k);

X_denoised = score_reduced * coeff_reduced'; % 降噪后的数据投射回原始空间
noise = X - X_denoised; % 噪声





