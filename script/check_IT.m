run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, threeExemplar, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier', experiment); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% Population average PSTH
clear opt;
fnd = load('Y:\External_DataSet\Okazawa2025_IT_data\visual_IT_both_task_unit_mod.mat').fnd_spc;

trc = trial_classifier('stim_group', {[0 5], [5 12], [12 24], [24 48], [48 Inf]}, 'plus_cat', 2);

[ID, mean_coh] = trc.stim_choice(fnd.getp('morph'), fnd.getp('targ_cho')); % correct trial only
% ID = trc.stim_category(fnd.getp('coh'), fnd.getp('targ_cho')); % all trials

psth_data = fnd.PSTH(ID, {'boxcar', 100});

opt.cutoff = fnd.cutoff();
%  In case of RT task: opt.cutoff = fnd.cutoff(ID); % different cutoff for each condition
opt.plot = set_plot_opt('vik', 10);
opt.legend = arrayfun(@(x) num2str(x, '%.0f'), mean_coh, 'uni', 0);
opt.event_label = {fnd.alignto.event};
fh = showPopPSTH(fnd.tstamp, psth_data, opt);

%%
fnd = load('Y:\External_DataSet\Okazawa2025_IT_data\visual_IT_both_task_unit_mod.mat').fnd_ex;

% get ID
morph_level = fnd.getp('morph');
targ_cho = fnd.getp('targ_cho');
targ_cor = fnd.getp('targ_cor');

trc = trial_classifier('stim_group', {[0 Inf]}, 'plus_cat', 2, 'include_0coh', false);
[ID, mean_coh] = trc.stim_choice(morph_level, targ_cho, targ_cor); % correct trials only
trial_classifier_result(ID, {'morph_level', 'targ_cor', 'targ_cho'}, {morph_level, targ_cor, targ_cho});

% plot PSTH
psth_data = fnd.PSTH(ID, {'boxcar', 100}, [], 0);

clear opt
opt.cutoff = fnd.cutoff();
opt.plot = set_plot_opt('roma', max(ID(:)));
opt.legend = {'Choice 1', 'Choice 2'};
opt.legend_pos = [0.01 0.25 0.1 0.5];
opt.event_label = {fnd.alignto.event};
fh = showPopPSTH(fnd.tstamp, psth_data, opt);


%%
fnd = load('Y:\External_DataSet\Okazawa2025_IT_data\visual_IT_both_task_unit_mod.mat').fnd_ex;
fnd_spc = load('Y:\External_DataSet\Okazawa2025_IT_data\visual_IT_both_task_unit_mod.mat').fnd_spc;

fnd = fnd.extract_epoch(1);

% get smooth detrended raster
r = fnd.raster(1); r = r{1}; % (unit, time, trial)
r_mn = mean(r, 3);
r_detrend = r - r_mn;

kernel = fspecial('average', [1, 100]);
r_detrend = nanconv(r_detrend, kernel, 'same');

% get ID
morph_level = fnd.getp('morph'); morph_level = morph_level(1,:)';
targ_cho = fnd.getp('targ_cho'); targ_cho = targ_cho(1,:)';
targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
correct = targ_cho==targ_cor;
ID = targ_cor; ID(~correct) = NaN;

% get trial-averaged baseline (correct trials only)
[nunit, ntime, ~] = size(r);
ncond = max(ID(:));

psth = nan(nunit, ntime, ncond);
for c = 1:ncond
    psth(:,:,c) = nanmean(r_detrend(:,:,ID==c), 3);
end
