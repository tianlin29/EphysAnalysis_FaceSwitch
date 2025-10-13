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

%% [RNN] randomly generalized RNN coordinates
sigma = 1;
a = normrnd(0, sigma, [1 10]);

cat1 = nan(10,4);
for n = 1:4
    cat1(:,n) = normrnd(a, sigma);
end

cat2 = nan(10,4);
for n = 1:4
    cat2(:,n) = normrnd(-a, sigma);
end

cat1 = cat1*4;
cat2 = cat2*4;

% mean
figure('Position', [50 100 200 200]); hold on
plot(a, 'LineWidth', 1, 'Color', 'red')
plot(-a, 'LineWidth', 1, 'Color', 'black')
format_panel(gcf, 'xlabel', 'Dimension', 'ylabel', 'Score')

% 4 pairs
figure('Position', [250 100 200 200]); hold on
plot(cat1, 'LineWidth', 1, 'Color', 'red')
plot(cat2, 'LineWidth', 1, 'Color', 'black')
format_panel(gcf, 'xlabel', 'Dimension', 'ylabel', 'Score')

pc = [cat1, cat2]'; % (stim, dimension)
save('D:\Engine_D2\RNN_modeling\FaceSwitch\data\preproc\generated_pc\four_pairs.mat', 'pc')

%% [RNN] randomly generalized RNN coordinates
sigma = 1;
a = nan(4,10);
for n = 1:4
    a(n,:) = normrnd(0, sigma, [1 10]);
end

cat1 = nan(10,4);
for n = 1:4
    cat1(:,n) = normrnd(a(n,:), sigma);
end

cat2 = nan(10,4);
for n = 1:4
    cat2(:,n) = normrnd(-a(n,:), sigma);
end

cat1 = cat1*4;
cat2 = cat2*4;

% mean
figure('Position', [50 100 600 200]);
for n = 1:4
    subplot(1,4,n); hold on
    plot(a(n,:), 'LineWidth', 1, 'Color', 'red')
    plot(-a(n,:), 'LineWidth', 1, 'Color', 'black')
end
format_panel(gcf, 'xlabel', 'Dimension', 'ylabel', 'Score')

% 4 pairs
figure('Position', [250 100 200 200]); hold on
plot(cat1, 'LineWidth', 1, 'Color', 'red')
plot(cat2, 'LineWidth', 1, 'Color', 'black')
format_panel(gcf, 'xlabel', 'Dimension', 'ylabel', 'Score')

pc = [cat1, cat2]'; % (stim, dimension)
save('D:\Engine_D2\RNN_modeling\FaceSwitch\data\preproc\generated_pc\four_pairs_ungeneralizable.mat', 'pc')

%% [RNN] randomly generalized RNN coordinates
dim = 10;
num_vectors = 4;
cat1 = create_orthogonal_vectors(dim, num_vectors) * 3;
cat2 = -cat1;

% plot mean
figure('Position', [50 100 600 200]);
for n = 1:4
    subplot(1,4,n); hold on
    plot(cat1(:,n), 'LineWidth', 1, 'Color', 'red')
    plot(cat2(:,n), 'LineWidth', 1, 'Color', 'black')
end
format_panel(gcf, 'xlabel', 'Dimension', 'ylabel', 'Score')

% plot 4 pairs
figure('Position', [250 100 200 200]); hold on
plot(cat1, 'LineWidth', 1, 'Color', 'red')
plot(cat2, 'LineWidth', 1, 'Color', 'black')
format_panel(gcf, 'xlabel', 'Dimension', 'ylabel', 'Score')

% check orthogonality
dot_product = dot(cat1(:,1), cat1(:,4));
fprintf('dot product: %.10f\n', dot_product);

% save
pc = [cat1, cat2]'; % (stim, dimension)
save('D:\Engine_D2\RNN_modeling\FaceSwitch\data\preproc\generated_pc\four_pairs_ungeneralizable_orth.mat', 'pc')


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

%% predict FR from face pc (TBU)
%        stim_cat targ_cor stim_id
% human  2        1        1-20
% monkey 1        2        21-40

% face pc
pcscore_human = load(fullfile(MainInterimDir, experiment, 'pcscore_human.mat')).pcscore_human; % (20 faces, 50 face pc)
pcscore_human = pcscore_human(:, [1 2 3 26 27 28]); % (20 faces, face pc)

% FR
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_trial(fnd.getp('stim_cat')==2);
r = fnd.FR({1, [50 400]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

[nunit, ntime, ntrial] = size(fnd.data{1});
condID = fnd.getp('stim_id');
psth = fnd.PSTH(condID, {'boxcar', 100}); psth = psth{1};

r = fnd.FR({1, [100 400]}, condID); % (unit, condition)
FR = fnd.FR({1, [100 400]})'; % (trial, unit) = (observation, feature)

% PCA
[coeff, score, latent] = pca(psth(:,:)');
score


coeff = pca(r'); % coeff ..(unit, pc)
dim = 1:6;
coeff = coeff(:,dim);

% regression
% coeff (unit, pc) = w * pcscore_human (stim, face pc)

X = pcscore_human(:,1);
for n = 1:6 % pc
    y = coeff(:,n);
    b = regress(y', X');
end


b = nan(nunit,7);
stats = nan(nunit,4);
X = pcscore_human(stim_id,:);
for u = 1:nunit
    y = FR(:,u);
    [b(u,:),bint,r,rint,stats(u,:)] = regress(y, [ones(size(X,1),1), X]);
end

% plot
fh = figure('Position', [150 500 300 400]);
imagesc(b(:,2:end))
xlabel('face PC'); xticklabels({'1', '2', '3', '26', '27', '28'})
ylabel('unit')
colorbar
print(fh, '-dpdf', fullfile(FigDir, sprintf('regression_%s_%s_session%d.pdf', monkey, experiment, n)));

%% manifold (face pc) 看不出规律
pcscore_human = load(fullfile(MainInterimDir, experiment, 'pcscore_human_share_rand40.mat')).pcscore_human; % (stim, face pc)
pcscore_monkey = load(fullfile(MainInterimDir, experiment, 'pcscore_monkey_share_rand40.mat')).pcscore_monkey; % (stim, face pc)

pcscore_human = pcscore_human(:,1); % pc 1
pcscore_monkey = pcscore_monkey(:,1);

for n = 1:n_files
    % load data
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(1);

    % select unit
    r = fnd.FR({1, [100 500]});
    I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % select trial
    fnd = fnd.extract_trial(fnd.getp('stim_id')<=20); % human face only
    fnd = fnd.extract_trial(fnd.getp('targ_cor')==fnd.getp('targ_cho')); % correct trials only

    % get ID
    %        stim_cat targ_cor stim_id
    % human  2        1        1-20
    % monkey 1        2        21-40
    [nunit, ntime, ntrial] = size(fnd.data{1});
    stim_id = fnd.getp('stim_id'); stim_id = stim_id(1,:);
    face_pc = nan(1, ntrial);
    for tr = 1:ntrial
        if stim_id(tr)<=20
            face_pc(tr) = pcscore_human(stim_id(tr));
        else
            face_pc(tr) = pcscore_monkey(stim_id(tr)-20);
        end
    end

    nbin = 10;
    face_pc_list = prctile(face_pc, linspace(0,100,nbin+1));
    ID = nan(1, ntrial);
    for tr = 1:ntrial
        idx = find(face_pc_list>=face_pc(tr));
        ID(tr) = idx(1);
    end
    mean_coh = nan(1, nbin);
    for b = 1:nbin
        mean_coh(b) = mean(face_pc(ID==b));
    end
    ID = repmat(ID, [nunit, 1]);

    stim_cat = fnd.getp('stim_cat');
    targ_cho = fnd.getp('targ_cho');
    targ_cor = fnd.getp('targ_cor');
    trial_classifier_result(ID, {'stim_cat', 'targ_cor', 'targ_cho'}, {stim_cat, targ_cor, targ_cho});

    % plot manifold
    data = fnd.PSTH(ID, {'boxcar', 100});

    opt.mean_coh = mean_coh;
    opt.plot = set_plot_opt('vik', max(ID(:)));
    opt.PC_range = [250 600];
    opt.Time = [100 200 300 400 500 600];
    opt.epoch = 1;

    opt.roughness = 5e-4;
    opt.dim_sign = [1 -1 1];
    opt.dim = [3 1 2];
    opt.view = [-145 34];
    fh = showPCA(fnd.tstamp, data, opt);

    % change format
    format_panel(gcf, 'fig_size', [1270 250])
    a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',6)
    h = get(gca,'xlabel'); set(h, 'FontSize', 7); h = get(gca,'ylabel'); set(h, 'FontSize', 7); h = get(gca,'zlabel'); set(h, 'FontSize', 7);
    print(fh, '-dpdf', fullfile(FigDir, sprintf('manifold_facePC_%s_%s_session%d.pdf', monkey, experiment, n)));
end


