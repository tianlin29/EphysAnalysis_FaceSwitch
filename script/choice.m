run('../Initialize.m');
monkey = 'Nick';
experiment = 'learnTask2';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'choice'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'choice'); mkdir(InterimDir);

%%
D = combineTrialData(fullfile(PreprocDir, sprintf('%s_%s.mat', monkey, experiment)));

%% unsigned choice
session = cellfun(@(x) x.session_id, D);
cond = cellfun(@(x) x.cond, D);
coh = cellfun(@(x) x.coh, D);
targ_cor = cellfun(@(x) x.targ_cor, D);
resp = cellfun(@(x) x.resp, D);
cor = resp==targ_cor;

clear opt
% select data
opt.session_list = 1:n_files;
% process data
opt.log = true;
opt.constant = false;
opt.verbose = false;
% plot
switch experiment
    case 'learnTask2'
        opt.color = [0 0 0; 44 145 224; 44 145 224]/255;
    case 'learnTask3'
        opt.color = [0 0 0; 58 191 153; 58 191 153]/255;
    case 'learnTask4'
        opt.color = [0 0 0; 240 169 58; 240 169 58]/255;
end
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_unsigned_choice_2cond(cond, coh, cor, session, opt);
fh_all = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
save(fullfile(InterimDir, sprintf('unsigned_choice_%s_%s.mat', monkey, experiment)), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('unsigned_choice_summary_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('unsigned_choice_%s_%s.pdf', monkey, experiment)));

%% plot threshold of all sessions together
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

fh_summary = cell(2, 3); % {monkey, experiment}
for exp_id = 1:length(experiment_list)
    for monkey_id = 1:length(monkey_list)
        monkey = monkey_list{monkey_id};
        experiment = experiment_list{exp_id};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end
        D = combineTrialData(fullfile(PreprocDir, sprintf('%s_%s.mat', monkey, experiment)));
        [~, n_files] = get_file_path(monkey, experiment);

        session = cellfun(@(x) x.session_id, D);
        cond = cellfun(@(x) x.cond, D);
        coh = cellfun(@(x) x.coh, D);
        targ_cor = cellfun(@(x) x.targ_cor, D);
        resp = cellfun(@(x) x.resp, D);
        cor = resp==targ_cor;

        clear opt
        % select data
        opt.session_list = 1:n_files;
        % process data
        opt.log = true;
        opt.constant = false;
        opt.verbose = false;
        % plot
        switch experiment
            case 'learnTask2'
                opt.color = [0 0 0; 44 145 224; 44 145 224]/255;
            case 'learnTask3'
                opt.color = [0 0 0; 58 191 153; 58 191 153]/255;
            case 'learnTask4'
                opt.color = [0 0 0; 240 169 58; 240 169 58]/255;
        end
        opt.average = false; % do not average across learning sessions

        [~, ~, fh_summary{monkey_id, exp_id}, ~] = run_unsigned_choice_2cond(cond, coh, cor, session, opt);
    end
end

fh_all = plot_in_one_fig_threshold(fh_summary, [2 3], [300 200]*1.2);
print(fh_all, '-dpdf', fullfile(FigDir, 'threshold.pdf'));

%% [stat] accuracy of unsigned choice
stat = load(fullfile(InterimDir, sprintf('unsigned_choice_%s_%s.mat', monkey, experiment))).stat;
nses = length(stat);
acc1 = cellfun(@(x) x.acc1, stat);
acc2 = cellfun(@(x) x.acc2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', acc1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Accuracy of pair 1 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Accuracy of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
else
    error('check')
end
[rho, p_val] = corr((1:nses)', acc2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Accuracy of pair 2 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Accuracy of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
else
    error('check')
end

%% [stat] threshold of unsigned choice
stat = load(fullfile(InterimDir, sprintf('unsigned_choice_%s_%s.mat', monkey, experiment))).stat;
nses = length(stat);
thres1 = cellfun(@(x) x.thres1, stat);
thres2 = cellfun(@(x) x.thres2, stat);

fprintf('Spearman''s rank correlation coefficient:\n')
[rho, p_val] = corr((1:nses)', thres1', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Threshold of pair 1 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Threshold of pair 1 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Threshold of pair 1 decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end
[rho, p_val] = corr((1:nses)', thres2', 'Type', 'Spearman');
if p_val>=0.05
    fprintf('Threshold of pair 2 did not change with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho>0
    fprintf('Threshold of pair 2 increased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
elseif p_val<0.05 && rho<0
    fprintf('Threshold of pair 2 decreased with learning: rho = %.2f, %s\n', rho, p2str(p_val))
end

%% choice
session = cellfun(@(x) x.session_id, D);
cond = cellfun(@(x) x.cond, D);
coh = cellfun(@(x) x.coh, D);
resp = cellfun(@(x) x.resp, D);

clear opt
% select data
opt.session_list = 1:n_files;
% process data
opt.constant = false;
opt.verbose = false;
% plot
switch experiment
    case 'learnTask2'
        opt.color = [99 97 172; 254 160 64]/255;
    case 'learnTask3'
        opt.color = [99 97 172; 242 128 128]/255;
    case 'learnTask4'
        opt.color = [99 97 172; 178 34 34]/255;
end
opt.average = false; % do not average across learning sessions

[~, fh_idv, fh_summary, stat] = run_choice_2cond(cond, coh, resp, session, opt);
fh_all = plot_in_one_fig_choice_2cond(fh_idv, [3 5], [500 300]*1.5);
save(fullfile(InterimDir, sprintf('choice_%s_%s.mat', monkey, experiment)), 'stat')
print(fh_summary, '-dpdf', fullfile(FigDir, sprintf('choice_summary_%s_%s.pdf', monkey, experiment)));
print(fh_all, '-dpdf', fullfile(FigDir, sprintf('choice_%s_%s.pdf', monkey, experiment)));

%% plot choice for the first session of each pair
file_path = {'\\10.10.49.250\rawdata\FaceSwitch_Monkey\Nick\20250528\Nick20250528_01.mat', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Nick\20250704\Nick20250704_01.mat', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Nick\20250721\Nick20250721_01.mat';
    '', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Woody\20231103\Woody20231103_01.mat', '\\10.10.49.250\rawdata\FaceSwitch_Monkey\Woody\20240110\Woody20240110_01.mat'};
monkey_list = {'Nick', 'Woody'};
experiment_list = {'learnTask2', 'learnTask3', 'learnTask4'};

for monkey_idx = 1:2
    for experiment_idx = 1:3
        monkey = monkey_list{monkey_idx};
        experiment = experiment_list{experiment_idx};
        if strcmp(monkey, 'Woody') && strcmp(experiment, 'learnTask2')
            continue;
        end

        D = load(file_path{monkey_idx, experiment_idx}).trial_data;
        D = D(~isnan(cellfun(@(x) x.response, D))); % responded trial only

        session = ones(size(D));
        cond = cellfun(@(x) x.curr_task_set, D);
        coh = cellfun(@(x) x.morph_level(1), D);
        resp = cellfun(@(x) 3-x.response, D);

        clear opt
        % select data
        opt.session_list = 1;
        % process data
        opt.constant = false;
        opt.verbose = false;
        % plot
        switch experiment_idx
            case 1
                opt.color = [0 0 0; 44 145 224]/255;
            case 2
                opt.color = [0 0 0; 58 191 153]/255;
            case 3
                opt.color = [0 0 0; 240 169 58]/255;
        end
        opt.average = false; % do not average across learning sessions
        [~, fh_idv, ~, stat] = run_choice_2cond(cond, coh, resp, session, opt);
        print(fh_idv, '-dpdf', fullfile(FigDir, sprintf('choice_of_the_first_session_%s_%s.pdf', monkey, experiment)));
    end
end












