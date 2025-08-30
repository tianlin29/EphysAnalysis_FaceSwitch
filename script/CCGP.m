run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'ccgp_courellis2024', 'code', 'helper_functions')));
monkey = 'Nick';
experiment = 'learnTask4';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'CCGP'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'CCGP'); mkdir(InterimDir);

%% transfrom FND to suitable data format
file_path = arrayfun(@(x) get_file_path(monkey, experiment, x, 'FND_sorted'), 1:n_files, 'uni', 0);

neu_fnd.array = {};
neu_fnd.param = {};
neu_fnd.cellinfo = [];
neu_fnd.sessionID = {};
for n = 1:length(file_path)
    fprintf('session %d\n', n)

    % load data
    fnd = load(file_path{n}).fnd;
    fnd = fnd.extract_epoch(2);

    % remove units
    FR = fnd.FR({1, [50 600]});
    I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);

    % arrange neural data
    [nunit, ntime, ntrial] = size(fnd.data{1});
    raster = fnd.raster(1); raster = raster{1}; % (unit, time, trial)
    tstamp = fnd.tstamp{1};
    fr_stim = squeeze(mean(raster(:, tstamp>50 & tstamp<600, :), 2)*1000); % (unit, trial)
    for u = 1:nunit
        neu_fnd.array = cat(2, neu_fnd.array, {table(fr_stim(u,:)', 'VariableNames', {'fr_stim'})});
    end

    % arrange trial parameter
    context = fnd.getp('task_set'); context = context(1,:)';
    iscorrect = fnd.getp('targ_cho')==fnd.getp('targ_cor'); iscorrect = iscorrect(1,:)';
    response = fnd.getp('targ_cho'); response = response(1,:)';
    reward = abs(fnd.getp('morph_level'))>0.3; reward = reward(1,:)'; % actually it is easy vs. difficult
    neu_fnd.param = cat(2, neu_fnd.param, {table(iscorrect, context, reward, response)});

    % arrange unit and session info
    neu_fnd.cellinfo = cat(2, neu_fnd.cellinfo, ones(1,nunit)); % 1 ..PFC neurons
    neu_fnd.sessionID = cat(2, neu_fnd.sessionID, repmat({sprintf('session_%d', n)}, [1,nunit]));
end

save(fullfile(InterimDir, sprintf('neu_fnd_%s_%s.mat', monkey, experiment)), 'neu_fnd')

%% set options
neu_fnd = load(fullfile(InterimDir, sprintf('neu_fnd_%s_%s.mat', monkey, experiment))).neu_fnd;

% choose session
sessions = neu_fnd.sessionID; % (1, cell); session 1, session 2, ...

% choose brain area
cellinfo = neu_fnd.cellinfo; % (1, cell); 1 ..PFC neurons
idx_pfc = find(ismember(cellinfo,[1]));
cell_groups = {idx_pfc};
area_order = {'vlPFC'};

% set sample number
    % n_sample=15  ..For single neuron recording. Pick 15 trials in each
    %                condition for each brain area.
    % n_sample=Inf ..For population neuron recording. Use almost all trials
    %                in each condition for each brain area (trial number
    %                would be matched across conditions). 
n_samples = [Inf]; 

    % n_resample=1000 ..For single neuron recording. Resample for multiple
    %                   times to get mean decoding accuracy.
    % n_resample=1    ..For population neuron recording. Since all trials
    %                   would be used, no need to resample.
n_resample = 1; 

    % fix it to 1 since neuron resampling happens outside of sd, ccgp, ps
    % methods 
n_perm_inner = 1; 

%% run geometric analyses
rng(1); % random seed

n_sessions = n_files;
n_variable = 3; % 2 ..task*choice; 3 ..task*choice*stimulus difficulty
switch n_variable
    case 2
        set_function = 'define_sets_cond4';
    case 3
        set_function = 'define_sets';
end
[sd_,sd_boot,ccgp_,ccgp_boot,ps_,ps_boot] = deal(cell(length(cell_groups),n_sessions)); % {cell, session}; (dichotomy, resample)

for idx_brain = 1:length(cell_groups) % for every brain area
    for idx_session = 1:n_sessions % for every session; For Nick learnTask2, start from session 6
        fprintf('session %d\n', idx_session)
        idx_cell = intersect(cell_groups{idx_brain},find(ismember(sessions,{sprintf('session_%d', idx_session)}))); % get cells according to brain area and session

        for i_rs = 1:n_resample
            [avg_array,~] = construct_regressors(neu_fnd,n_samples(idx_brain),idx_cell,idx_session,n_variable); % {1, cell}; {1, cond}; (trial, 1)

            % decoding accuracy
            [perf,perf_boot,~] = sd(avg_array,n_perm_inner,n_samples(idx_brain),set_function);
            sd_{idx_brain,idx_session} = cat(2,sd_{idx_brain,idx_session},perf);
            sd_boot{idx_brain,idx_session}  = cat(2,sd_boot{idx_brain,idx_session},perf_boot);

            % CCGP
            [ccgp_{idx_brain,idx_session}]  = cat(2,ccgp_{idx_brain,idx_session},...
                ccgp(avg_array, n_perm_inner,false,n_samples(idx_brain),set_function));
        end
    end
end

%% run geometric analyses for null distribution
    % for each dichotomy, pick n_samples=Inf (actually about 150) trials in
    % each condition for classification. Repeat the process for
    % n_perm_inner_CCGP_null=1000 times to get the null distribution.  
n_perm_inner_CCGP_null = 1000;

n_sessions = n_files;
n_variable = 3; % 2 ..task*choice; 3 ..task*choice*stimulus difficulty
switch n_variable
    case 2
        set_function = 'define_sets_cond4';
    case 3
        set_function = 'define_sets';
end
[ccgp_boot] = deal(cell(length(cell_groups),n_sessions)); % {cell, session}; (dichotomy, resample)

for idx_brain = 1:length(cell_groups)% for every brain area
    for idx_session = 1:n_sessions % for every session
        fprintf('session %d\n', idx_session)
        idx_cell = intersect(cell_groups{idx_brain},find(ismember(sessions,{sprintf('session_%d', idx_session)}))); % get cells according to brain area and sessions

        for i_rs = 1:n_resample
            [avg_array,~] = construct_regressors(neu_fnd,n_samples(idx_brain),idx_cell,idx_session,n_variable); % {1, cell}; {1, cond}; (trial, 1)

            % CCGP
            ccgp_boot_tmp = ccgp(avg_array,n_perm_inner_CCGP_null,true,n_samples(idx_brain),set_function);
            ccgp_boot{idx_brain,idx_session} = cat(2,ccgp_boot{idx_brain,idx_session},ccgp_boot_tmp);
        end
    end
end

%% save results
save(fullfile(InterimDir, sprintf('stat_CCGP_%s_%s.mat', monkey, experiment)), 'sd_', 'ccgp_');
save(fullfile(InterimDir, sprintf('stat_CCGP_null_%s_%s.mat', monkey, experiment)), 'sd_boot', 'ccgp_boot');

%% load results
load(fullfile(InterimDir, sprintf('stat_CCGP_%s_%s.mat', monkey, experiment)));
load(fullfile(InterimDir, sprintf('stat_CCGP_null_%s_%s.mat', monkey, experiment)));

%% plot decoding accuracy
% idx_special = [1 2 3];
% lgd_special = {'pair','choice','parity'};
idx_special = [1 10 21 29];
lgd_special = {'pair','difficulty','choice','parity'};
to_plt      = vec(cellfun(@(x) mean(x,2),sd_,  'UniformOutput',false)');
null_plt    = vec(cellfun(@(x) vec(x),sd_boot, 'UniformOutput',false)');

null_plt = [];
connection_map = ones(size(to_plt,1),1);
figure('units','normalized','position',[0.1 0.2 0.6 0.5])
fn_plot_HP_boxSwarmVariable(gca,to_plt,idx_special,lgd_special,null_plt,connection_map,[],'sd')
prepaxis(gca,{'','Decoding Accuracy',''})
xticks(1:2:n_files); xlabel('#Session');
ylim([0.4 1])
set(gcf, 'Position', [0.12 0.11 0.3 0.25])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('decoding_accuracy_%s_%s.pdf', monkey, experiment)));

%% plot decoding accuracy (Nick learnTask2)
% idx_special = [1 2 3];
% lgd_special = {'pair','choice','parity'};
idx_special = [1 10 21 29];
lgd_special = {'pair','difficulty','choice','parity'};
to_plt      = vec(cellfun(@(x) mean(x,2),sd_,  'UniformOutput',false)');
null_plt    = vec(cellfun(@(x) vec(x),sd_boot, 'UniformOutput',false)');

to_plt(1:5) = [];

null_plt = [];
connection_map = ones(size(to_plt,1),1);
figure('units','normalized','position',[0.1 0.2 0.6 0.5])
fn_plot_HP_boxSwarmVariable(gca,to_plt,idx_special,lgd_special,null_plt,connection_map,[],'sd')
prepaxis(gca,{'','Decoding Accuracy',''})
xticks(1:2:n_files); xlabel('#Session'); xticklabels({'6', '8', '10', '12', '14'})
ylim([0.4 1])
set(gcf, 'Position', [0.12 0.11 0.3 0.25])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('decoding_accuracy_%s_%s.pdf', monkey, experiment)));

%% plot cross-condition generalization performance
% idx_special = [1 2 3];
% lgd_special = {'pair','choice','parity'};
idx_special = [1 10 21 29];
lgd_special = {'pair','difficulty','choice','parity'};
to_plt      = vec(cellfun(@(x) mean(x,2),ccgp_,  'UniformOutput',false)');
null_plt    = vec(cellfun(@(x) vec(x),ccgp_boot, 'UniformOutput',false)');

% null_plt = [];
connection_map = ones(size(to_plt,1),1);
figure('units','normalized','position',[0.1 0.2 0.6 0.5])
fn_plot_HP_boxSwarmVariable(gca,to_plt,idx_special,lgd_special,null_plt,connection_map,[],'ccgp')
prepaxis(gca,{'','CCGP',''})
xticks(1:2:n_files); xlabel('#Session');
ylim([0 1])
set(gcf, 'Position', [0.12 0.11 0.3 0.25])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('CCGP_%s_%s.pdf', monkey, experiment)));

%% plot cross-condition generalization performance (Nick learnTask2)
% idx_special = [1 2 3];
% lgd_special = {'pair','choice','parity'};
idx_special = [1 10 21 29];
lgd_special = {'pair','difficulty','choice','parity'};
to_plt      = vec(cellfun(@(x) mean(x,2),ccgp_,  'UniformOutput',false)');
null_plt    = vec(cellfun(@(x) vec(x),ccgp_boot, 'UniformOutput',false)');

to_plt(1:5) = [];

null_plt = [];
connection_map = ones(size(to_plt,1),1);
figure('units','normalized','position',[0.1 0.2 0.6 0.5])
fn_plot_HP_boxSwarmVariable(gca,to_plt,idx_special,lgd_special,null_plt,connection_map,[],'ccgp')
prepaxis(gca,{'','CCGP',''})
xticks(1:2:n_files); xlabel('#Session'); xticklabels({'6', '8', '10', '12', '14'})
ylim([0 1])
set(gcf, 'Position', [0.12 0.11 0.3 0.25])
print(gcf, '-dpdf', fullfile(FigDir, sprintf('CCGP_%s_%s.pdf', monkey, experiment)));

%% [stat] how does CCGP of choice/pair change with learning
ccgp_pair = cellfun(@(x) x(1), ccgp_);
ccgp_choice = cellfun(@(x) x(2), ccgp_);
fprintf('\n%s %s:\n', monkey, experiment)

% stat
fprintf('\nSpearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', ccgp_pair', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Pair CCGP did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Pair CCGP increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Pair CCGP decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

fprintf('\nSpearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:n_files)', ccgp_choice', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Choice CCGP did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Choice CCGP increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Choice CCGP decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end


