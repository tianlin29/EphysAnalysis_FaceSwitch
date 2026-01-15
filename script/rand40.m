run('../Initialize.m');
monkey = 'Nick';
experiment = 'rand40';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

%% coordinates of rand40 faces in the shared face space
% ID of rand40 faces
fname_face40 = load('C:\Engine\rome_external\PassiveFixationMask\randomFace40\process.mat').fname;
info_human = load(fullfile(MainDir, 'data\info\rand40\Passive_fixation_faces_Human.mat')).struc; % pcscore, name
info_monkey = load(fullfile(MainDir, 'data\info\rand40\Passive_fixation_faces_Monkey.mat')).struc;

% coordinates of 189 faces in the shared face space
faces_org = load('C:\Engine\FacePCspace\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.faces_org;
scaling_factor = load('C:\Engine\FacePCspace\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.scaling_factor;
faces_org = faces_org .* scaling_factor';

% summary
data = cell(40,5); % [stim id, name in rand40, id in the shared sapce, name in the shared space of 189 faces]
for i = 1:40
    data{i,1} = i;
    data{i,2} = fname_face40{i};
    if i<=20
        data{i,3} = info_human.img_idx(strcmp(info_human.name, fname_face40{i}));
        data{i,4} = sprintf('human%03d.png', data{i,3});
    else
        data{i,3} = info_monkey.img_idx(strcmp(info_monkey.name, fname_face40{i})) + 109;
        data{i,4} = sprintf('monkey%03d.png', data{i,3}-109);
    end
    data{i,5} = faces_org(data{i,3}, :);
end
save(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat')), 'data')

%% coordinates of 4pair+rand40 faces in the shared face space
% coordinates of 189 faces in the shared face space
faces_org = load('C:\Engine\FacePCspace\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.faces_org;
scaling_factor = load('C:\Engine\FacePCspace\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.scaling_factor;
faces_org = faces_org .* scaling_factor';

% combine 4pair and rand40
idx_human = [109, 76, 89, 71];
idx_monkey = [80, 4, 38, 32];
pc_human_4pair = faces_org(idx_human,:);
pc_monkey_4pair = faces_org(idx_monkey+109,:);

data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;
pc_human_rand40 = cell2mat(data(1:20,5));
pc_monkey_rand40 = cell2mat(data(21:end,5));

pc_4pair_rand40 = [pc_human_4pair; pc_human_rand40; pc_monkey_4pair; pc_monkey_rand40];

% MDS
RSM = corrcoef(pc_4pair_rand40'); % correlation between conditions
RDM = 1-RSM; % dissimilarity matrix
MDS_2d = mdscale(RDM, 2, 'Criterion', 'metricstress'); % get the first 2 dims

ImgName_human_4pair = arrayfun(@(x) sprintf('human%03d.png', x), idx_human', 'uni', 0);
ImgName_monkey_4pair = arrayfun(@(x) sprintf('monkey%03d.png', x), idx_monkey', 'uni', 0);
ImgName_human_rand40 = data(1:20, 4);
ImgName_monkey_rand40 = data(21:end, 4);
ImgName = [ImgName_human_4pair; ImgName_human_rand40; ImgName_monkey_4pair; ImgName_monkey_rand40];

ImgPath_face189 = 'C:\Engine\FacePCspace\human_monkey_face_space\data';
color_list = [0 0 0; 44 145 224; 58 191 153; 240 169 58]/255;
fh = figure('Position', [50 100 300 300]); hold on;
colormap(gray);
plot_list = [(1:20)+4, (21:40)+4+4, 1:4, 25:28];
for i = plot_list
    % load image
    img = imread(fullfile(ImgPath_face189, ImgName{i}));

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.0006*1.6/2; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    if ismember(i, [1 2 3 4 25 26 27 28])
        rectangle('Position', [MDS_2d(i,1)-width/2, MDS_2d(i,2)-height/2, width, height], 'EdgeColor', color_list(mod(i, 24), :), 'LineWidth', 5);
    end
    image([MDS_2d(i,1)-width/2, MDS_2d(i,1)+width/2], [MDS_2d(i,2)-height/2, MDS_2d(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');

end
format_panel(fh, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', 'xlim')
axis equal

%% randomize coordinates of 4pair+rand40 faces in the shared face space
rng(4) % 4

% coordinates of 189 faces in the shared face space
faces_org = load('C:\Engine\FacePCspace\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.faces_org;
scaling_factor = load('C:\Engine\FacePCspace\human_monkey_face_space\FaceSwitch_all_faces_PCA.mat').fspace.scaling_factor;
faces_org = faces_org .* scaling_factor';

% combine 4pair and rand40
idx_human = [109, 76, 89, 71];
idx_monkey = [80, 4, 38, 32];
pc_human_4pair = faces_org(idx_human,:);
pc_monkey_4pair = faces_org(idx_monkey+109,:);

data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;
pc_human_rand40 = cell2mat(data(1:20,5));
pc_monkey_rand40 = cell2mat(data(21:end,5));

pc_4pair_rand40 = [pc_human_4pair; pc_human_rand40; pc_monkey_4pair; pc_monkey_rand40];

% plot
mn = [mean(pc_4pair_rand40(1:24,:)); mean(pc_4pair_rand40(25:end,:))]; se = [std(pc_4pair_rand40(1:24,:)); std(pc_4pair_rand40(25:end,:))];
figure('Position', [100 100 520 180]); 
subplot(1,2,1); hold on
plot(mn(1,:), '.-', 'Color', 'black')
cerrorbar(1:size(pc_4pair_rand40,2), mn(1,:), se(1,:), 'Color', 'black');
plot(mn(2,:), '.-', 'Color', 'red')
cerrorbar(1:size(pc_4pair_rand40,2), mn(2,:), se(2,:), 'Color', 'red');
plot(xlim, [0 0], ':', 'Color', 'black')
format_panel(gca, 'xlabel', '#PC', 'ylabel', 'Score', 'xtick', [1 25 50], 'axis', 'normal')
title('Unnormalized')

% randomize
pc_4pair_rand40 = pc_4pair_rand40 * normrnd(0, 2, 50, 50);

% plot
mn = [mean(pc_4pair_rand40(1:24,:)); mean(pc_4pair_rand40(25:end,:))]; se = [std(pc_4pair_rand40(1:24,:)); std(pc_4pair_rand40(25:end,:))];
subplot(1,2,2); hold on
plot(mn(1,:), '.-', 'Color', 'black')
cerrorbar(1:size(pc_4pair_rand40,2), mn(1,:), se(1,:), 'Color', 'black');
plot(mn(2,:), '.-', 'Color', 'red')
cerrorbar(1:size(pc_4pair_rand40,2), mn(2,:), se(2,:), 'Color', 'red');
plot(xlim, [0 0], ':', 'Color', 'black')
format_panel(gca, 'xlabel', '#PC', 'ylabel', 'Score', 'xtick', [1 10 25 50], 'axis', 'normal')
title('Randomize')

% save
pc = pc_4pair_rand40(:,1:10); % (stim, dimension)
save(fullfile(InterimDir, sprintf('pc_4pair_rand40.mat')), 'pc')
% save('D:\FaceSwitch\code\RNN_modeling\ContextSwitch\data\preproc\pc_4pair_rand40.mat', 'pc')

%% confirm correlation between face40 and face189
data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;

ImgPath_face40 = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim';
ImgPath_face189 = 'C:\Engine\FacePCspace\human_monkey_face_space\data';
coordinate = [[1:20,1:20]', [2*ones(20,1); 1*ones(20,1)]];
struct_data.ImgName_face40 = data(:,2);
struct_data.ImgName_face189 = data(:,4);

for face = {'face40', 'face189'}
    fh = figure('Position', [50 100 1800 200]); hold on;
    colormap(gray);
    for i = 1:40
        % load image
        img = imread(fullfile(eval(['ImgPath_', face{1}]), struct_data.(['ImgName_', face{1}]){i}));

        % determine width and height
        img_width = size(img, 2);
        img_height = size(img, 1);
        scale = 0.0006*1.7*3.5; % 0.0006 for z-scored; 0.000015 for un-z-scored
        width = img_width * scale;
        height = img_height * scale;

        % plot image at location
        image([coordinate(i,1)-width/2, coordinate(i,1)+width/2], [coordinate(i,2)-height/2, coordinate(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
    end
    axis equal
    title(face)
end

%% plot images based on 2-d face MDS
CAT = 'monkey';

data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;
switch CAT
    case 'human'
        pc = cell2mat(data(1:20, 5));
        ImgName = data(1:20, 2);
    case 'monkey'
        pc = cell2mat(data(21:end, 5));
        ImgName = data(21:end, 2);
end

RSM = corrcoef(pc'); % correlation between conditions
RDM = 1-RSM; % dissimilarity matrix
MDS_2d = mdscale(RDM, 2, 'Criterion', 'metricstress'); % get the first 2 dims

ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';
fh = figure('Position', [50 100 300 300]); hold on;
colormap(gray);
for i = 1:20
    % load image
    img = imread(fullfile(ImgPath, ImgName{i}));

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.0006*1.6/2; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    image([MDS_2d(i,1)-width/2, MDS_2d(i,1)+width/2], [MDS_2d(i,2)-height/2, MDS_2d(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
end
format_panel(fh, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', 'xlim')
axis equal
print(fh, '-dpdf', fullfile(FigDir, sprintf('MDS_image_face_%s.pdf', CAT)));

%% plot images based on 2-d face MDS (both human and monkey)
data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;
pc = cell2mat(data(:, 5));
ImgName = data(:, 2);

RSM = corrcoef(pc'); % correlation between conditions
RDM = 1-RSM; % dissimilarity matrix
MDS_2d = mdscale(RDM, 2, 'Criterion', 'metricstress'); % get the first 2 dims

ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';
fh = figure('Position', [50 100 300 300]); hold on;
colormap(gray);
for i = 1:40
    % load image
    img = imread(fullfile(ImgPath, ImgName{i}));

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.0006*1.6/2; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    image([MDS_2d(i,1)-width/2, MDS_2d(i,1)+width/2], [MDS_2d(i,2)-height/2, MDS_2d(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
end
format_panel(fh, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', 'xlim')
axis equal
print(fh, '-dpdf', fullfile(FigDir, sprintf('MDS_image_face.pdf')));


%% plot images based on 2-d neural MDS
CAT = 'human';

n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(1);
data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;
switch CAT
    case 'human'
        fnd = fnd.extract_trial(fnd.getp('stim_id')<=20);
        ID = fnd.getp('stim_id');
        ImgName = data(1:20, 2);
    case 'monkey'
        fnd = fnd.extract_trial(fnd.getp('stim_id')>20);
        ID = fnd.getp('stim_id')-20;
        ImgName = data(21:end, 2);
end

FR = fnd.FR({1, [300 500]}, ID, true); % (unit, condition)
RSM = corrcoef(FR); % correlation between conditions
RDM = 1-RSM; % dissimilarity matrix
MDS_2d = mdscale(RDM, 2, 'Criterion', 'metricstress'); % get the first 2 dims

ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';
fh = figure('Position', [50 100 300 300]); hold on;
colormap(gray);
for i = 1:20
    % load image
    img = imread(fullfile(ImgPath, ImgName{i}));

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.0006*1.6/2; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    image([MDS_2d(i,1)-width/2, MDS_2d(i,1)+width/2], [MDS_2d(i,2)-height/2, MDS_2d(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
end
format_panel(fh, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', 'xlim')
axis equal
print(fh, '-dpdf', fullfile(FigDir, sprintf('MDS_image_neural_%s.pdf', CAT)));

%% plot images based on 2-d neural MDS (both human and monkey)
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(1);
data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;

ID = fnd.getp('stim_id');
ImgName = data(:, 2);

FR = fnd.FR({1, [300 500]}, ID, true); % (unit, condition)
RSM = corrcoef(FR); % correlation between conditions
RDM = 1-RSM; % dissimilarity matrix
MDS_2d = mdscale(RDM, 2, 'Criterion', 'metricstress'); % get the first 2 dims

ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';
fh = figure('Position', [50 100 300 300]); hold on;
colormap(gray);
for i = 1:40
    % load image
    img = imread(fullfile(ImgPath, ImgName{i}));

    % determine width and height
    img_width = size(img, 2);
    img_height = size(img, 1);
    scale = 0.0006*1.6/2; % 0.0006 for z-scored; 0.000015 for un-z-scored
    width = img_width * scale;
    height = img_height * scale;

    % plot image at location
    image([MDS_2d(i,1)-width/2, MDS_2d(i,1)+width/2], [MDS_2d(i,2)-height/2, MDS_2d(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
end
format_panel(fh, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', 'xlim')
axis equal
print(fh, '-dpdf', fullfile(FigDir, sprintf('MDS_image_neural.pdf')));

%% similarity of RDM of neural data and shared face PC
CAT = 'monkey';
[RDM, MDS_2d] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(1);

    % get FR
    switch CAT
        case 'human'
            fnd = fnd.extract_trial(fnd.getp('stim_id')<=20);
            ID = fnd.getp('stim_id');
        case 'monkey'
            fnd = fnd.extract_trial(fnd.getp('stim_id')>20);
            ID = fnd.getp('stim_id')-20;
    end
    FR = fnd.FR({1, [300 500]}, ID, true); % (unit, condition)

    % get neural RDM
    RSM = corrcoef(FR); % correlation between conditions
    RDM{n} = 1-RSM; % dissimilarity matrix
    MDS_2d{n} = mdscale(RDM{n}, 2, 'Criterion', 'metricstress'); % get the first 2 dims
end

% similarity among neural RDMs
p_value_neural = nan(n_files, n_files);
for n = 1:n_files
    for m = n:n_files
        [r_obs, p_value_neural(n,m)] = mantel_test(RDM{n}, RDM{m}, 10000, 'both');
        p_value_neural(m,n) = p_value_neural(n,m);
    end
end

% get face PC
data = load(fullfile(InterimDir, sprintf('shared_face_space_rand40.mat'))).data;

switch CAT
    case 'human'
        main_pc = cell2mat(data(1:20,5));
    case 'monkey'
        main_pc = cell2mat(data(21:end,5));
end

% get MDS of face PC
face_pc = main_pc'; % (dim, condition)
RSM = corrcoef(face_pc); % correlation between conditions
RDM_face = 1-RSM; % dissimilarity matrix
MDS_2d = mdscale(RDM_face, 2, 'Criterion', 'metricstress'); % get the first 2 dims

% similarity between neural RDM and face PC RDM
p_value_face = nan(1, n_files);
for n = 1:n_files
    [r_obs, p_value_face(n)] = mantel_test(RDM_face, RDM{n}, 10000, 'both');
end

% plot RDMs
figure('Position', [100 100 970 310]);
for n = 1:n_files
    subplot(2,5,n); hold on
    imagesc(RDM{n}, [0 2]); colorbar; set(gca, 'YDir', 'reverse'); title({'Neural RDM', sprintf('Session %d', n)})
    format_panel(gca, 'xlabel', '#Face', 'ylabel', '#Face', 'xtick', [1 10 20], 'ytick', [1 10 20])
end
subplot(2,5,5); hold on
imagesc(RDM_face, [0 2]); colorbar; set(gca, 'YDir', 'reverse'); title('Face PC RDM')
format_panel(gca, 'xlabel', '#Face', 'ylabel', '#Face', 'xtick', [1 10 20], 'ytick', [1 10 20])

subplot(2,5,6)
imagesc(p_value_neural, [0 0.05]); colorbar; title({'P-value', 'between sessions'})
format_panel(gca, 'xlabel', '#Session', 'ylabel', '#Session', 'xtick', 1:4, 'ytick', 1:4)

subplot(2,5,7)
imagesc(p_value_face, [0 0.05]); colorbar; title({'P-value', 'between neuron and face PC'})
format_panel(gca, 'xlabel', '#Session', 'xtick', 1:4, 'axis', 'equal', 'ylim', [-2.5 2.5])

subplot(2,5,8)
text(0.5, 0.5, sprintf('%s %s (%s)', monkey, experiment, CAT), 'FontSize', 10); axis off
colormap jet

%%
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted_preprocessed')).fnd;

clear opt;
opt.time_win = {{1, [50 100]}, {1, [100 150]}, {1, [150 200]}, {1, [200 250]}, ...
                {2, [200 250]}, {2, [250 300]}, {2, [300 350]}, {2, [350 400]}, ...
                {3, [-250 -200]}, {3, [-200 -150]}, {3, [-150 -100]}, {3, [-100 -50]}}; % define time windows
% opt.time_win = {{1, [150 200]}, {3, [-200 -150]}}; % for report
opt.plot = set_plot_opt('vik', 40);
opt.ID_method = 'cat_1'; % one category only
run_mds = runMDS(opt);

% dissimilarity matrix
[fh_sim, fh_mds, stat] = run_mds.getRDM(fnd);
% print(fh_sim, '-dpdf', [FigDir 'fig3_similarity.pdf']);
% print(fh_mds, '-dpdf', [FigDir 'fig3_MDS.pdf']);

% face images
% category = 1; % human category
% run_mds.plotImage(stat.mds_2d, FigDir, category);

% Mantel test
[fh_sim, fh_p] = run_mds.getCorrRDM(stat.RDM); % whether dissimilarity matrix is similar across timepoints
% print(fh_sim, '-dpdf', [FigDir 'fig3_structure_is_similar.pdf']);

% FR matrix 
fh_FR = run_mds.plotFR(stat.FR); % plot FR
print(fh_FR, '-dpdf', [FigDir 'fig3_FR.pdf']);
[fh_sim, fh_p] = run_mds.getCorrFR(stat.FR); % whether FR matrix is similar across timepoints
print(fh_sim, '-dpdf', [FigDir 'fig3_FR_is_not_similar.pdf']);

%% PCA (stim id)
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(1);
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











