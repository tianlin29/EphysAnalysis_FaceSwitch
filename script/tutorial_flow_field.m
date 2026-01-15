run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, faceColor_passiveLong, threeExemplar
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'flowField'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'flowField'); mkdir(InterimDir);

%% ###### Load data (necessary for all steps below) #######
%% load data
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);

%% get PC space to plot
morph_level = fnd.getp('morph_level');
task_set = fnd.getp('task_set');
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');
correct = double(fnd.getp('targ_cho')==fnd.getp('targ_cor'));

ID = targ_cho;
ID(task_set == 2) = ID(task_set == 2) + max(ID(:));
ID(correct~=1) = NaN;
trial_classifier_result(ID, {'task_set', 'targ_cho', 'correct'}, {task_set, targ_cho, correct})

psth = fnd.PSTH(ID, {'gaussian', 20});

clear opt;
opt.PC_range = [250 600];
opt.epoch = 1;
[coef, detPSTH] = getPC_coef(fnd.tstamp, psth, opt);


%% ###### Flow field only (effective dynamics) #######
%% calculate flow field
clear opt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.dim = [1 2];
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;
psth = calc_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, opt);


%% show flow field
clear pltopt;
pltopt.min_sample = 1000;
pltopt.trj_color = [
    0.00  0.45  0.74;
    0.85  0.33  0.10;
    0.00  0.60  0.50;
    0.80  0.60  0.00];
pltopt.legend = {'Task 1 C1', 'Task 1 C2', 'Task 2 C1', 'Task 2 C2'};
pltopt.xlabel = 'PC1';
pltopt.ylabel = 'PC2';
fh = show_flow_field(psth, pltopt);

%% 1 vs 3 dim
opt.dim = [1 3];
psth = calc_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, opt);

pltopt.ylabel = 'PC3';
fh = show_flow_field(psth, pltopt);


%% ###### Flow field + Input (Input drive + Internal dynamics) #######
%% define external input
clear opt;
opt.latency = 100; % ms
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);

stim_dur = fnd.getp('stim_dur') * 1e3;
morph_level = fnd.getp('morph_level');
task_set = fnd.getp('task_set');
task_set(task_set(:) == 2) = -1; % 1 vs -1

input = gen_input_matrix(fnd.tstamp{1}, stim_dur(1,:), {morph_level(1,:), task_set(1,:)}, opt);

%% calculate flow field
clear opt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.dim = [1 2];
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;

psth = calc_input_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, input, opt);

%% show flow field
clear pltopt;
pltopt.min_sample = 1000;
pltopt.trj_color = [
    0.00  0.45  0.74;
    0.85  0.33  0.10;
    0.00  0.60  0.50;
    0.80  0.60  0.00];
pltopt.legend = {'Task 1 C1', 'Task 1 C2', 'Task 2 C1', 'Task 2 C2'};
pltopt.input_col = [0 0 0;.8 .8 .8];
pltopt.xlabel = 'PC1';
pltopt.ylabel = 'PC2';
fh = show_flow_field(psth, pltopt);
title('K arrow: morph, G arrow: task');

%% flow field in PC1-3 dimension
opt.dim = [1 3];

psth = calc_input_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, input, opt);

pltopt.ylabel = 'PC3';
fh = show_flow_field(psth, pltopt);
title('K arrow: morph, G arrow: task');


%% ###### Fit Linear Dynamical Systems to Flow Field #######
%% calculate flow field
clear opt pltopt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.dim = [1 2];
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;
psth = calc_flow_field(fnd.tstamp{1}, fnd.raster{1}, ID, opt);

%% fit LDS

lds_data = fit_LDS_to_flow_field(psth);

%lds_data = fit_LDS_to_flow_field_Nd(data, 2);

pltopt.min_sample = 1000;
pltopt.trj_color = [
    0.00  0.45  0.74;
    0.85  0.33  0.10;
    0.00  0.60  0.50;
    0.80  0.60  0.00];
pltopt.legend = {'Task 1 C1', 'Task 1 C2', 'Task 2 C1', 'Task 2 C2'};
pltopt.xlabel = 'PC1';
pltopt.ylabel = 'PC2';
fh = show_flow_field_w_LDS(psth, lds_data, pltopt);


%% ###### Flow field in 3D slice #######
%% calculate flow field at each 2D plane along PC 3 dimensioon
clear opt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.flow_field_dim = [1 2];
opt.slice_dim = 3;
opt.window_size = 0.2;  % ratio of trials in each window
opt.window_step = 0.02; % % ratio trials step (100 steps total)
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;

psth = calc_flow_field_slice(fnd.tstamp{1}, fnd.raster{1}, ID, opt);


%% show flow field
clear pltopt;

pltopt.min_sample = 1000;
pltopt.trj_color = [
    0.00  0.45  0.74;
    0.85  0.33  0.10;
    0.00  0.60  0.50;
    0.80  0.60  0.00];
pltopt.legend = {'Task 1 C1', 'Task 1 C2', 'Task 2 C1', 'Task 2 C2'};
pltopt.xlabel = 'PC1';
pltopt.ylabel = 'PC2';
pltopt.n_slices_to_show = 5;
fh = show_flow_field_slice(psth, pltopt);


%% fit LDS for each slice

lds_data = fit_LDS_to_flow_field(psth);

fh = show_LDS_fit_accuracy(lds_data);



