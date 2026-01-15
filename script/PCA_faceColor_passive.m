run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'faceColor_passiveLong'; % learnTask2, learnTask3, learnTask4, faceColor, faceColor_passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'PCA'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'PCA'); mkdir(InterimDir);

%% classifier
file_list = [3];
n_files = length(file_list);

REPEAT = 1;
for r = 1:REPEAT
    data = cell(n_files, 1); % (session, condition)
    for n = 1:n_files
        fprintf('repreat %d, session %d\n', r, n)

        % set options
        clear opt;
        opt.epoch = 1;
        opt.tstamp = {[400]}; % 对于不同位置的时刻点，使用的是不同的随机种子，因此结果会不同
        opt.t_win = 200;
        opt.glm_type = 'lasso';

        % load neural data
        fnd = load(get_file_path(monkey, experiment, file_list(n), 'FND_sorted')).fnd;
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        task_set = fnd.getp('task_set');
        targ_cor = fnd.getp('targ_cor');

        % build pair 1 classifier
        ID = targ_cor; ID(task_set~=1) = NaN;
        opt.Kfold = 1;
        data_pair1_axis = PopulationClassifier(fnd, ID, opt);

        % project pair 2 to pair 1
        ID = targ_cor; ID(task_set~=2) = NaN;
        opt.Kfold = 1;
        opt.refBeta = data_pair1_axis.data_ind; % attach pair 1 classification result
        data{n,1} = PopulationClassifier(fnd, ID, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r)), 'data', 'opt')
end

%% analyze classifier
task_set = fnd.getp('task_set'); task_set = task_set(1,:);
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:);
color_set_id = fnd.getp('color_set_id'); color_set_id = color_set_id(1,:);
color_set_full = color_set_id; color_set_full(targ_cor==2) = color_set_full(targ_cor==2) + max(color_set_full);

target = color_set_full(task_set==2); target = target';
correct = data{1}.Correct{1};

I = ~isnan(correct);
[p1, pse1] = calcGroupMean(correct(I), target(I), unique(target(I)), 'binary');

x = 0:22.5:360-22.5;
fh = figure('Position', [50 100 150 200]); hold on
plot(x, p1, '.-');
plot([45 45], ylim, ':', 'Color', 'black')
plot([45 45]+180, ylim, ':', 'Color', 'black')
format_panel(gcf, 'xlabel', 'Angle on color ring (deg)', 'ylabel', 'Category decoding accuracy')
print(fh, '-dpdf', fullfile(FigDir, sprintf('classifier_face_16color_%s_%s_session%d.pdf', monkey, experiment, n)));

%% PCA for faceColor_passiveLong (2 faces + 16 colors)
n = 3;

% load data
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);

% select unit
r = fnd.FR({1, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% get ID
task_set = fnd.getp('task_set');
targ_cor = fnd.getp('targ_cor');
color_set_id = fnd.getp('color_set_id');

ID = color_set_id; % color set 1-8
ID(task_set==1) = 9; % face image
ID(targ_cor==2) = ID(targ_cor==2)+max(ID(:));
trial_classifier_result(ID, {'targ_cor', 'color_set_id', 'task_set'}, {targ_cor, color_set_id, task_set});

% plot PCA
clear opt
opt.plot = set_plot_opt_2cond('roma', 'roma', max(ID(:))/2);
Lxy = [0.3449    0.2956;
    0.3598    0.3274;
    0.3627    0.3585;
    0.3539    0.3840;
    0.3324    0.3994;
    0.2992    0.4007;
    0.2679    0.3896;
    0.2385    0.3602;
    NaN       NaN
    0.2285    0.3112;
    0.2136    0.2794;
    0.2106    0.2482;
    0.2195    0.2227;
    0.2410    0.2073;
    0.2741    0.2061;
    0.3054    0.2172;
    0.3348    0.2465;
    NaN       NaN];
Lxy = [10*ones(size(Lxy,1),1), Lxy];
calib_file = 'C:\Users\HP\Desktop\FaceSwitch\code\FaceColor_Monkey\code_to create colored face monkey\calib_data_19-Nov-2024a.mat';
rgb = Lxy2rgb(Lxy, calib_file);
rgb([9 18],:) = [0 0 0; 0 0 0]; % change to black
rgb = rgb/255;
opt.plot.color = rgb;
opt.plot.facecolor = rgb;
opt.plot.edgecolor = rgb;

opt.epoch = 1;
opt.PC_range = [250 600];
opt.PSTH_conv = {'boxcar', 100};
opt.PC_kernel = {'bartlett', 200};

fh = showPCA_trajectory(fnd, ID, opt);
saveas(fh, fullfile(FigDir, sprintf('PCA_face_16color_%s_%s_session%d.fig', monkey, experiment, n)));
print(fh, '-dpdf', fullfile(FigDir, sprintf('PCA_face_16color_%s_%s_session%d.pdf', monkey, experiment, n)));

%% PCA for faceColor_passiveLong (16 colors)
n = 3;

% load data
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);

% select unit
r = fnd.FR({1, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd = fnd.set_unit_criteria('custom', I);

% get ID
task_set = fnd.getp('task_set');
targ_cor = fnd.getp('targ_cor');
color_set_id = fnd.getp('color_set_id');

ID = color_set_id; % color set 1-8
ID(task_set==1) = NaN; % remove face image
ID(targ_cor==2) = ID(targ_cor==2)+max(ID(:));
trial_classifier_result(ID, {'targ_cor', 'color_set_id', 'task_set'}, {targ_cor, color_set_id, task_set});

% plot PCA
clear opt
opt.plot = set_plot_opt_2cond('roma', 'roma', max(ID(:))/2);
Lxy = [0.3449    0.2956;
    0.3598    0.3274;
    0.3627    0.3585;
    0.3539    0.3840;
    0.3324    0.3994;
    0.2992    0.4007;
    0.2679    0.3896;
    0.2385    0.3602;
    0.2285    0.3112;
    0.2136    0.2794;
    0.2106    0.2482;
    0.2195    0.2227;
    0.2410    0.2073;
    0.2741    0.2061;
    0.3054    0.2172;
    0.3348    0.2465];
Lxy = [10*ones(size(Lxy,1),1), Lxy];
calib_file = 'C:\Users\HP\Desktop\FaceSwitch\code\FaceColor_Monkey\code_to create colored face monkey\calib_data_19-Nov-2024a.mat';
rgb = Lxy2rgb(Lxy, calib_file);
rgb = rgb/255;
opt.plot.color = rgb;
opt.plot.facecolor = rgb;
opt.plot.edgecolor = rgb;

opt.epoch = 1;
opt.PC_range = [250 600];
opt.PSTH_conv = {'boxcar', 100};
opt.PC_kernel = {'bartlett', 200};

fh = showPCA_trajectory(fnd, ID, opt);
