
FigDir = [MainFigDir 'SimpleFlowField/'];
mkdir(FigDir);

% create a very simple flow field plot (this is very crude toy analysis!)

%% ###### Load data (necessary for all steps below) #######

%% load data
clear opt;
S = load([PreprocDir 'example_face_switch_array.mat']);
fnd = S.fnd; % FND class

%% get PC space to plot
morph = fnd.getp('morph_level') * 1e2;

trc = trial_classifier('stim_group', {[0 Inf]}, 'plus_cat', 1);

ID = trc.stim_choice(morph, fnd.getp('targ_cho')); % correct trial only
task = fnd.getp('task_set');
tID = ID;
tID(task == 2) = tID(task == 2) + 2; % make it task dependent (task x choice PSTH)

data = fnd.PSTH(tID, {'gaussian', 20});

opt.plot = set_plot_opt('vik', 4);
opt.PC_range = [250 600];
opt.conv_kernel = fspecial('average', [1 100]);
opt.epoch = 2;

[coef, detPSTH] = getPC_coef(fnd.tstamp, data, opt);


%% ###### Flow field only (effective dynamics) #######


%% calculate flow field
clear opt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.dim = [1 2];
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;
data = calc_flow_field(fnd.tstamp{2}, fnd.raster{2}, tID, opt);


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
fh = show_flow_field(data, pltopt);

print(fh, '-dpdf', [FigDir 'flow_field_PC12.pdf']);

%% 1 vs 3 dim

opt.dim = [1 3];
data = calc_flow_field(fnd.tstamp{2}, fnd.raster{2}, tID, opt);

pltopt.ylabel = 'PC3';
fh = show_flow_field(data, pltopt);

print(fh, '-dpdf', [FigDir 'flow_field_PC13.pdf']);


%% ###### Flow field + Input (Input drive + Internal dynamics) #######

%% define external input
clear opt;

opt.latency = 100; %ms
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);

sdur = fnd.getp('stim_dur') * 1e3;
morph = fnd.getp('morph_level');
task = fnd.getp('task_set');
task(task(:) == 2) = -1; % 1 vs -1

input = gen_input_matrix(fnd.tstamp{2}, sdur(1,:), {morph(1,:), task(1,:)}, opt);

%% calculate flow field

clear opt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.dim = [1 2];
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;

data = calc_input_flow_field(fnd.tstamp{2}, fnd.raster{2}, tID, input, opt);

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
fh = show_flow_field(data, pltopt);
title('K arrow: morph, G arrow: task');

print(fh, '-dpdf', [FigDir 'input_flow_field_PC12.pdf']);

%% flow field in PC1-3 dimension
opt.dim = [1 3];

data = calc_input_flow_field(fnd.tstamp{2}, fnd.raster{2}, tID, input, opt);

pltopt.ylabel = 'PC3';
fh = show_flow_field(data, pltopt);
title('K arrow: morph, G arrow: task');

print(fh, '-dpdf', [FigDir 'input_flow_field_PC13.pdf']);


%% ###### Fit Linear Dynamical Systems to Flow Field #######

%% calculate flow field
clear opt pltopt;
opt.PC_coef = coef;
opt.detPSTH = detPSTH;
opt.dim = [1 2];
opt.conv_kernel = fspecial('gaussian', [1 20*6], 20);
opt.bin_num = 20;
data = calc_flow_field(fnd.tstamp{2}, fnd.raster{2}, tID, opt);

%% fit LDS

lds_data = fit_LDS_to_flow_field(data);

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
fh = show_flow_field_w_LDS(data, lds_data, pltopt);

print(fh, '-dpdf', [FigDir 'LDS2D_fit_to_flow_field.pdf']);

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

data = calc_flow_field_slice(fnd.tstamp{2}, fnd.raster{2}, tID, opt);


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
fh = show_flow_field_slice(data, pltopt);

print(fh, '-dpdf', [FigDir 'flow_field_slice.pdf']);

%% fit LDS for each slice

lds_data = fit_LDS_to_flow_field(data);

fh = show_LDS_fit_accuracy(lds_data);
print(fh, '-dpdf', [FigDir 'LDS_fit_accuracy_slice.pdf']);



