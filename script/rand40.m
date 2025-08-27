run('../Initialize.m');
monkey = 'Nick';
experiment = 'rand40';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

%% [preprocess] coordinates of rand40 in the individual human/monkey face space
% 一共制作了41张人类面孔和41张猴子面孔，从中随机挑选除了20张人类面孔和20张猴子面孔
fname = load('C:\Engine\rome_external\PassiveFixationMask\randomFace40\process.mat').fname; % {40,1}
info_human = load('..\data_face_switch\info\rand40\Passive_fixation_faces_Human.mat').struc; % pcscore, name
info_monkey = load('..\data_face_switch\info\rand40\Passive_fixation_faces_Monkey.mat').struc;

[pcscore_human, pcscore_monkey] = deal(nan(20, 50)); % (20 stim, 50 dimension)
for n = 1:20
    I = strcmp(info_human.name, fname{n});
    pcscore_human(n,:) = info_human.pcscore{I};
    I = strcmp(info_monkey.name, fname{n+20});
    pcscore_monkey(n,:) = info_monkey.pcscore{I};
end
save(fullfile(MainInterimDir, experiment, 'pcscore_human.mat'), 'pcscore_human');
save(fullfile(MainInterimDir, experiment, 'pcscore_monkey.mat'), 'pcscore_monkey');

%% [preprocess] coordinates of rand40 in the shared human-monkey face space
% 人猴共用空间中共有189张脸（109张人脸和80张猴脸），找到rand40中的40张脸
fname = load('C:\Engine\rome_external\PassiveFixationMask\randomFace40\process.mat').fname;
faces_org = load('D:\Engine_D2\RNN_modeling\FaceSwitch\data\preproc\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.faces_org; % (189 stim, 50 dimension)

[pcscore_human, pcscore_monkey] = deal(nan(20, 50)); % (20 stim, 50 dimension)
for n = 1:20
    idx = str2double(fname{n}([6 7])); % e.g., 'human03.bmp' -> 3
    pcscore_human(n,:) = faces_org(idx,:);
    idx = str2double(fname{n+20}([7 8]));  % e.g., 'monkey01.bmp' -> 1
    pcscore_monkey(n,:) = faces_org(idx+109,:);
end
save(fullfile(MainInterimDir, experiment, 'pcscore_human_share_rand40.mat'), 'pcscore_human');
save(fullfile(MainInterimDir, experiment, 'pcscore_monkey_share_rand40.mat'), 'pcscore_monkey');

%% [preprocess] coordinates of pair 1 to 5 in the shared human-monkey face space
faces_org = load('D:\Engine_D2\RNN_modeling\FaceSwitch\data\preproc\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.faces_org; % (189 stim, 50 dimension)
idx_human = [109, 76, 89, 71, 2]; % pair 1 to 5
idx_monkey = [80, 4, 38, 32, 21];

[pcscore_human, pcscore_monkey] = deal(nan(5, 50)); % (5 stim, 50 dimension)
for n = 1:5
    pcscore_human(n,:) = faces_org(idx_human(n),:);
    pcscore_monkey(n,:) = faces_org(idx_monkey(n)+109,:);
end
save(fullfile(MainInterimDir, experiment, 'pcscore_human_share_pari1to5.mat'), 'pcscore_human');
save(fullfile(MainInterimDir, experiment, 'pcscore_monkey_share_pari1to5.mat'), 'pcscore_monkey');

%% training parameters
% Trial structure is similar to the main task. Only stimulus set was
% changed to randomFace40. All stimuli were prototypes.

%        stim_cat targ_cor stim_id
% human  2        1        1-20
% monkey 1        2        21-40

fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;

stim_id = fnd.getp('stim_id'); stim_id = stim_id(1,:)'; unique(stim_id); % 40 stim id
stim_cat = fnd.getp('stim_cat'); stim_cat = stim_cat(1,:)'; unique(stim_cat); % 1 or 2
tmp = stim_id(stim_cat==1); unique(tmp); % for stim_cat=1, stim_id is 21-40
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)'; unique(targ_cor); % 1 or 2, targ_cor is flipped compared to stim_cat
set_id = fnd.getp('set_id'); set_id = set_id(1,:)'; unique(set_id); % always 2, this variable could be meaningless
check = table(stim_id, stim_cat, targ_cor); % human is category 2, target 1

%% ** plot human faces in human face space
% get PC coordinates
pcscore_human = load(fullfile(MainInterimDir, experiment, 'pcscore_human.mat')).pcscore_human;
main_pc = pcscore_human(:, [1 26]);

% get iamges
fname = load('C:\Engine\rome_external\PassiveFixationMask\randomFace40\process.mat').fname; % (40 stim, 1)
ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';

% plot
n_image = size(main_pc,1);
fh = figure('Position', [50 100 150 150]); hold on;
colormap(gray);
for i = 1:n_image
    % load image
    img = imread([ImgPath, fname{i}]);

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.003; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    image([main_pc(i,1)-width/2, main_pc(i,1)+width/2], [main_pc(i,2)-height/2, main_pc(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
end
format_panel(fh, 'xlabel', 'Feature PC 1', 'ylabel', 'Appearance PC 1', 'axis', 'equal') % make axis equal so that faces are not stretched
print(fh, '-dpdf', fullfile(FigDir, 'human_in_face_space.pdf'));

%% ** plot monkey faces in monkey face space
% get PC coordinates
pcscore_monkey = load(fullfile(MainInterimDir, experiment, 'pcscore_monkey.mat')).pcscore_monkey;
main_pc = pcscore_monkey(:, [1 26]);

% get iamges
fname = load('C:\Engine\rome_external\PassiveFixationMask\randomFace40\process.mat').fname; % (40 stim, 1)
ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';

% plot
n_image = size(main_pc,1);
fh = figure('Position', [50 100 150 150]); hold on;
colormap(gray);
for i = 1:n_image
    % load image
    img = imread([ImgPath, fname{i + 20}]); % load monkey faces

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.003; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    image([main_pc(i,1)-width/2, main_pc(i,1)+width/2], [main_pc(i,2)-height/2, main_pc(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
end
format_panel(fh, 'xlabel', 'Feature PC 1', 'ylabel', 'Appearance PC 1', 'axis', 'equal')
print(fh, '-dpdf', fullfile(FigDir, 'monkey_in_face_space.pdf'));

%% plot rand40 faces and pair 1-5 faces in the shared face space
pcscore_human = load(fullfile(MainInterimDir, experiment, 'pcscore_human_share_rand40.mat')).pcscore_human;
pcscore_monkey = load(fullfile(MainInterimDir, experiment, 'pcscore_monkey_share_rand40.mat')).pcscore_monkey;
pcscore_human_ = load(fullfile(MainInterimDir, experiment, 'pcscore_human_share_pari1to5.mat')).pcscore_human;
pcscore_monkey_ = load(fullfile(MainInterimDir, experiment, 'pcscore_monkey_share_pari1to5.mat')).pcscore_monkey;

dim = [1 26]; % feature PC 1, appearance PC 1
pcscore_human = pcscore_human(:, dim); pcscore_human_ = pcscore_human_(:, dim);
pcscore_monkey = pcscore_monkey(:, dim); pcscore_monkey_ = pcscore_monkey_(:, dim);

h = nan(4,1);
fh = figure('Position', [50 100 270 150]); hold on;
h(1) = scatter(pcscore_human(:,1), pcscore_human(:,2), 5, 'blue', 'o');
h(2) = scatter(pcscore_human_(:,1), pcscore_human_(:,2), 5, 'blue', 'o', 'filled', 'MarkerEdgeColor', 'black');
h(3) = scatter(pcscore_monkey(:,1), pcscore_monkey(:,2), 5, 'red', 'o');
h(4) = scatter(pcscore_monkey_(:,1), pcscore_monkey_(:,2), 5, 'red', 'o', 'filled', 'MarkerEdgeColor', 'black');
plot([pcscore_human_(:,1)'; pcscore_monkey_(:,1)'], [pcscore_human_(:,2)'; pcscore_monkey_(:,2)'], 'Color', 'black', 'LineWidth', 0.5)
legend(h, {'Human (rand40)', 'Human (pair1-5)', 'Monkey (rand40)', 'Monkey (pair1-5)'}, 'Location', 'eastoutside');
if all(dim==[1 2])
    format_panel(fh, 'xlabel', 'Feature PC 1', 'ylabel', 'Feature PC 2', 'axis', 'normal')
    print(fh, '-dpdf', fullfile(FigDir, 'faces_rand40_pair1to5_view1.pdf'));
elseif all(dim==[1 26])
    format_panel(fh, 'xlabel', 'Feature PC 1', 'ylabel', 'Appearance PC 1', 'axis', 'normal')
    print(fh, '-dpdf', fullfile(FigDir, 'faces_rand40_pair1to5_view2.pdf'));
elseif all(dim==[2 27])
    format_panel(fh, 'xlabel', 'Feature PC 2', 'ylabel', 'Appearance PC 2', 'axis', 'normal')
end

%% PSTH (population)
fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;
ID = fnd.getp('stim_id');
psth_data = fnd.PSTH(ID, {'boxcar', 100});

clear opt
opt.cutoff = fnd.cutoff();
opt.plot = set_plot_opt('vik', max(ID(:)));
opt.event_label = {fnd.alignto.event};
fh = showPopPSTH(fnd.tstamp, psth_data, opt);
print(fh, '-dpdf', fullfile(FigDir, 'PopPSTH.pdf'));

%% ** PCA (stim id)
fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;
ID = fnd.getp('stim_id');

clear opt;
opt.plot = set_plot_opt('vik', max(ID(:)));
opt.plot.markersize = 2*ones(max(ID(:)),1); opt.plot.linewidth = 0.5*ones(max(ID(:)),1); % report format
opt.epoch = 1;
opt.PC_range = [250 600];
opt.PSTH_conv = {'boxcar', 100};
opt.PC_kernel = {'bartlett', 200};
fh = showPCA_trajectory(fnd, ID, opt);

% format
view([0 90])
format_panel(gcf, 'fig_size', [150 150], 'xlim', [-100 100], 'ylim', [-100 100])
set(gca, 'LineWidth', 0.5)
xtickangle(0)
print(fh, '-dpdf', fullfile(FigDir, 'PCA_stimID.pdf'));

%% predict FR from face pc
%        stim_cat targ_cor stim_id
% human  2        1        1-20
% monkey 1        2        21-40

% face pc
pcscore_human = load(fullfile(MainInterimDir, experiment, 'pcscore_human.mat')).pcscore_human; % (20 faces, 50 face pc)
nfeature = 2;
pcscore_human = pcscore_human(:, [1 26]); % (20 faces, face pc)

% FR
fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;
fnd = fnd.extract_trial(fnd.getp('stim_cat')==2);
r = fnd.FR({1, [50 400]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

[nunit, ntime, ntrial] = size(fnd.data{1});
stim_id = fnd.getp('stim_id'); stim_id = stim_id(1,:)'; % (trial, 1)
FR = fnd.FR({1, [50 400]})'; % (trial, unit) = (observation, feature)

% prepare data for regression
X = zeros(ntrial, nfeature);
for tr = 1:ntrial
    X(tr,:) = pcscore_human(stim_id(tr),:); % (trial, face pc)
end

Y = FR; % (trial, unit)

%
rng(1);
nfold = 5;
cv = cvpartition(ntrial, 'KFold', nfold);

[modelCorr, modelMSE] = deal(zeros(nunit, 1));
betaCoefs = zeros(nfeature, nunit);
for u = 1:nunit
    u

    y = Y(:, u);
    YPredict = zeros(size(y));
    for f = 1:cv.NumTestSets
        trainIdx = cv.training(f);
        testIdx = cv.test(f);

        XTrain = X(trainIdx, :); % (part of trial, face pc)
        yTrain = y(trainIdx); % (part of trial, 1)

        XTest = X(testIdx, :);

        % regression
        mdl = fitrlinear(XTrain, yTrain, 'Learner', 'leastsquares', ... 
                         'Regularization', 'ridge', ...
                         'Lambda', 'auto', ...
                         'Solver', 'sgd');

        % predict
        y_pred = predict(mdl, XTest);
        YPredict(testIdx) = y_pred;
    end

    % register
    modelCorr(u) = corr(y, YPredict);
    modelMSE(u) = mean((y - YPredict).^2);

end

% stat
meanCorr = mean(modelCorr);
meanMSE = mean(modelMSE);

fprintf('Average correlation between predicted and actual firing rate: %.4f\n', meanCorr);
fprintf('Average mean squared error: %.4f\n', meanMSE);

% stat 2
significantUnits = sum(modelCorr > 0.1);
fprintf('%d out of %d units (%.2f%%) were significantly encoded by the face space.\n', ...
        significantUnits, nunit, (significantUnits/nunit)*100);

% stat 3
figure;
subplot(1,2,1);
histogram(modelCorr, 20);
xlabel('Prediction Correlation');
ylabel('Number of Units');
title('Distribution of Encoding Model Performance');
xline(0.1, 'r--', 'LineWidth', 2); % 添加参考线

subplot(1,2,2);
boxplot(modelCorr);
ylabel('Prediction Correlation');
title('Performance Across Units');

% stat 4
% 选一个表现好（或不好）的神经元来可视化
[~, exampleUnit] = max(modelCorr); 

figure;
plot(y, YPredict, 'k.');
xlabel('Actual Firing Rate');
ylabel('Predicted Firing Rate');
title(sprintf('Unit %d: Corr = %.3f', exampleUnit, modelCorr(exampleUnit)));
lsline; % 添加拟合线


B = nan(nfeature+1, nunit); % (face pc, unit)
for u = 1:nunit
    [b, ~, ~, ~, stats] = regress(r(:,u), [ones(20,1), main_pc]); % (20 faces, u) = (20 faces, 19 face pc + 1) * (19 face pc + 1, u)
    B(:,u) = b(1:end);
end

r_pred = [ones(20,1), main_pc] * B;



