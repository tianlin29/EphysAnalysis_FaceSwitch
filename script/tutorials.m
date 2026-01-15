clear;

%% fitglm
load patients
tbl = table(Age, Weight, Gender, Smoker, VariableNames=["Age", "Weight", "Gender", "Smoker"]);
modelspec = "Smoker ~ Age*Weight*Gender - Age:Weight:Gender"; % 减去三阶交互项，只剩主效应和二阶交互项
mdl = fitglm(tbl, modelspec, Distribution="binomial");

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

%% Find Procrustes distance and plot superimposed shape
X = [40 88; 51 88; 35 78; 36 75; 39 72; 44 71; 48 71; 52 74; 55 77]; % (point, dim)
Y = [36 43; 48 42; 31 26; 33 28; 37 30; 40 31; 45 30; 48 28; 51 24]; % (point, dim)
figure; hold on
plot(X(:,1),X(:,2),"x")
plot(Y(:,1),Y(:,2),"o")
xlim([0 100])
ylim([0 100])
legend("Target shape (X)","Comparison shape (Y)")

[d,Z] = procrustes(X,Y);
% d ..squared Procrustes distance
% Z ..Procrustes transformation on Y

plot(Z(:,1),Z(:,2),"s")
legend("Target shape (X)","Comparison shape (Y)", ...
    "Transformed shape (Z)")
hold off

%% Analyze Prorustes transformation including rotation
rng("default")
n = 10;  
Y = normrnd(0,1,[n 2]); % (point, dim)

S = [cos(pi/3) -sin(pi/3); sin(pi/3) cos(pi/3)];
X = normrnd(0.5*Y*S+2,0.05,n,2);

figure; hold on
plot(X(:,1),X(:,2),"x")
plot(Y(:,1),Y(:,2),"o")
legend("Target shape (X)","Comparison shape (Y)")

[~,Z,transform] = procrustes(X,Y);
% transform.T ..(dim, dim) rotation or reflection
% transform.b ..(1) scalar
% transform.c ..(point, dim) translation

% move Y to X:
% Y' = b*Y*T + c

plot(Z(:,1),Z(:,2),"s")
legend("Target shape (X)","Comparison shape (Y)", ...
    "Transformed shape (Z)")
hold off




