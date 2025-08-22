run('../Initialize.m');
clc;

%% get params from fnd
fnd = load(get_file_path(FormattedDir, 'Woody', 'learnTask3', 1, 'FND')).fnd;
coh = fnd.getp('morph_level')*100; % (unit, trial)
pair = fnd.getp('task_set');
targ_cor = fnd.getp('targ_cor');
targ_cho = fnd.getp('targ_cho');

%% generate example data
coh = randsample([-50 -20 -10 -5 0 5 10 20 50]', 1e4, true); % coherence
task = randsample([1 2]', length(coh), true); % task

targ_cor = double(coh > 0) + 1; % correct target: positive = 2, negative = 1
targ_cor(coh==0) = randsample([1 2]', sum(coh==0), true); % for 0% coh, random

pcho = 1./ (1 + exp(-0.15 * coh));
targ_cho = double(pcho > rand(size(pcho))) + 1;

%% create classifier
% create an instance by specifying the properties of the trial_classifier class
trc = trial_classifier(); % default properties
trc = trial_classifier('stim_group', {[0 24], [24 48], [48 Inf]}, 'plus_cat', 1);
trc_0ex = trial_classifier('stim_group', {[0 24], [24 48], [48 Inf]}, 'plus_cat', 1, 'include_0coh', false); % remove 0% coh trials
trc_0individual = trial_classifier('stim_group', {0, [0 24], [24 48], [48 Inf]}, 'plus_cat', 1); % 0% coh is an independent category

% check
properties('trial_classifier') % check all the properties of the class
methods('trial_classifier') % check all the methods of the class
help trc.choice % check information about a function

%% classify trials based on choice
ID = trc.choice(targ_cho, coh);
fprintf('##### choice #####\n');
fprintf('classify trials based on choice\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on choice (correct trials only)
ID = trc.choice_correct(targ_cho, targ_cor, coh);
fprintf('##### choice_correct #####\n');
fprintf('classify trials based on choice (correct trials only)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on stimulus category
ID = trc.category(targ_cor, coh);
fprintf('##### category #####\n');
fprintf('classify trials based on stimulus category\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence (regardless of choice)
ID = trc.stim_category(coh, targ_cor);
fprintf('##### stim_category #####\n');
fprintf('classify trials based on coherence (regardless of choice)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

ID = trc_0individual.stim_category(coh, targ_cor); % 0% coh is an independent category
fprintf('##### stim_category #####\n');
fprintf('classify trials based on coherence (regardless of choice)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence and choice (correct trial only)
% if targ_cor is not specified, 0% coh trials will be removed
ID = trc.stim_choice(coh, targ_cho);
fprintf('##### stim_choice #####\n');
fprintf('classify trials based on coherence and choice (correct trial only)\n');
fprintf('\t*correct category is undefined for 0%% coherence\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

% specify targ_cor, 0% coh trials will be remained
ID = trc.stim_choice(coh, targ_cho, targ_cor);
fprintf('##### stim_choice #####\n');
fprintf('classify trials based on coherence and choice (correct trial only)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

% 0% coh trials will be seperated by choice
ID = trc_0individual.stim_choice(coh, targ_cho, targ_cor);
fprintf('##### stim_choice #####\n');
fprintf('classify trials based on coherence and choice (correct trial only)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence and choice (both correct and error trials)
ID = trc.stim_choice_all(coh, targ_cho);
fprintf('##### stim_choice_all #####\n');
fprintf('classify trials based on coherence and choice (both correct and error trials)s\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence and choice (error trial only)
ID = trc.stim_choice_error(coh, targ_cho);
fprintf('##### stim_choice_error #####\n');
fprintf('classify trials based on coherence and choice (error trial only)\n');
fprintf('\t*correct category is undefined for 0%% coherence\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence (choice 2 only)
ID = trc.stim_plus_choice(coh, targ_cho, targ_cor);
fprintf('##### stim_plus_choice #####\n');
fprintf('classify trials based on coherence (choice 2 only)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence (choice 1 only)
ID = trc.stim_minus_choice(coh, targ_cho, targ_cor);
fprintf('##### stim_minus_choice #####\n');
fprintf('classify trials based on coherence (choice 1 only)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on coherence (concatenate choice 1 and 2)
ID = trc.stim_plus_minus_choice(coh, targ_cho, targ_cor);
fprintf('##### stim_plus_minus_choice #####\n');
fprintf('classify trials based on coherence (concatenate choice 1 and 2)\n');
trial_classifier_result(ID, {'coh', 'targ_cor', 'targ_cho'}, {coh, targ_cor, targ_cho});

%% classify trials based on custom params (e.g., coh, targ_cor, and task)
ID = trc.stim_category(coh, targ_cor);
fprintf('##### stim_category #####\n');
fprintf('classify trials based on coherence (regardless of choice)\n');
ID(task==2) = ID(task==2) + max(ID(:));
trial_classifier_result(ID, {'coh', 'targ_cor', 'task'}, {coh, targ_cor, task})

%% statistics
ntr = trc.ntrial_per_category(ID); % number of trials in each classified groups
coh_mn = trc.mean_coh(ID, coh); % average coherence level of trials in each classified groups
pref_mn = trc.mean_performance(ID, targ_cor, targ_cho); % average performance of trials in each classified groups
