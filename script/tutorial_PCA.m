clear;


%% PCA�����ɷ�ϵ���������ͷ���
load hald
size(ingredients); % (13 observation, 4 feature) = (cond*time, unit)

[coeff, score, latent] = pca(ingredients);
% coeff ..(4 feature, 4 pc) ���ɷ�ϵ����ÿһ����ԭ�ռ��е�һ��PC�ᣬ������ͶӰ
% score ..(13 observation, 4 pc) ���ɷַ�����ÿһ����һ��ԭʼ���ݵ���PC�ռ��е����꣬Ҳ���ǽ�ά�������
% latent ..(4 pc, 1) ���ɷַ��ÿһ����һ�����ɷֵķ���

k = 3; % ѡ��ǰ�������ɷ�
score_reduced = score(:, 1:k);
coeff_reduced = coeff(:, 1:k);

X_denoised = score_reduced * coeff_reduced'; % ����������Ͷ���ԭʼ�ռ�
noise = X - X_denoised; % ����





