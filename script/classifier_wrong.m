run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, threeExemplar, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier', experiment); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% projection
[angle, dist] = deal(cell(n_files, 1));
for n = 1:n_files
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);

    % get smooth detrended raster
    r = fnd.raster(1); r = r{1}; % (unit, time, trial)
    r_mn = mean(r, 3);
    r_detrend = r - r_mn;

    kernel = fspecial('average', [1, 100]);
    r_detrend = nanconv(r_detrend, kernel, 'same');

    % get ID
    morph_level = fnd.getp('morph_level');
    task_set = fnd.getp('task_set');
    targ_cor = fnd.getp('targ_cor');
    targ_cho = fnd.getp('targ_cho');
    correct = targ_cor==targ_cho;

    %
    ntime = size(r,2);
    I1 = task_set(1,:)==1 & correct(1,:)==1;
    r_pair1 = r(:,:,I1); targ_cor_pair1 = targ_cor(1,I1);
    I2 = task_set(1,:)==2 & correct(1,:)==0;
    r_pair2 = r(:,:,I2); targ_cor_pair2 = targ_cor(1,I2); targ_cho_pair2 = targ_cho(1,I2);

    acc = nan(ntime, 1);
    for t = 1:ntime
        mdl = fitclinear(squeeze(r_pair1(:,t,:))', targ_cor_pair1'); % build decoder on training set
        y_hat = predict(mdl, squeeze(r_pair2(:,t,:))'); % predict using testing set
        acc(t) = mean(y_hat==targ_cor_pair2');
    end

end





%% classifier (category axis for main task; wrong trials only)
if strcmp(experiment, 'faceColor')
    file_list = [1:8, 18:23];
else
    file_list = 1:n_files;
end
n_files = length(file_list);

REPEAT = 1;
for r = 1:REPEAT
    data = cell(n_files, 2); % (session, condition)
    for n = 1:n_files
        fprintf('repreat %d, session %d\n', r, n)

        % set options
        clear opt;
        opt.epoch = 1;
        opt.tstamp = {[100,200,300,400,500,600]}; % 对于不同位置的时刻点，使用的是不同的随机种子，因此结果会不同
        opt.t_win = 200;
        opt.glm_type = 'lasso';

        % load neural data
        fnd = load(get_file_path(monkey, experiment, file_list(n), 'FND_sorted')).fnd;
        fnd = fnd.extract_epoch(2); % extract epoch 2

        % extract task variables
        task_set = fnd.getp('task_set');
        targ_cor = fnd.getp('targ_cor');
        targ_cho = fnd.getp('targ_cho');
        correct = targ_cor==targ_cho;

        % build pair 2 correct classifier
        ID = targ_cor; ID(task_set~=2 | ~correct) = NaN;
        opt.Kfold = 1;
        data_category_axis = PopulationClassifier(fnd, ID, opt);

        % project pair 2 wrong to pair 2 correct
        ID = targ_cor; ID(task_set~=2 | correct) = NaN;
        opt.Kfold = 1;
        opt.refBeta = data_category_axis.data_ind; % attach pair 1 classification result
        data{n,1} = PopulationClassifier(fnd, ID, opt);

        % build pair 1 correct classifier
        ID = targ_cor; ID(task_set~=1 | ~correct) = NaN;
        opt.Kfold = 1;
        opt.refBeta = [];
        data_category_axis = PopulationClassifier(fnd, ID, opt);

        % project pair 2 wrong to pair 1 correct
        ID = targ_cor; ID(task_set~=2 | correct) = NaN;
        opt.Kfold = 1;
        opt.refBeta = data_category_axis.data_ind; % attach pair 1 classification result
        data{n,2} = PopulationClassifier(fnd, ID, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d_wrong.mat', monkey, experiment, r)), 'data', 'opt')
end

%% plot threshold vs. session/timebin
% load data
n_repeats = 1;
n_conds = 2;
opt_classifier = load(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat1_wrong.mat', monkey, experiment))).opt;

[coh, cor, cond, session, repeat] = deal([]);
for r = 1:n_repeats
    data = load(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d_wrong.mat', monkey, experiment, r))).data; 
    for n = 1:n_files
        for c = 1:n_conds
            ntrial = length(data{n,c}.param.morph_level);
            coh = [coh; data{n,c}.param.morph_level];
            cor = [cor; data{n,c}.Correct{1}];
            cond = [cond; c*ones(ntrial,1)];
            session = [session; n*ones(ntrial,1)];
            repeat = [repeat; r*ones(ntrial,1)];
        end
    end
end

% plot
[fh_session, fh_timebin] = plot_threshold(coh, cor, cond, session, repeat, opt_classifier, experiment);
% print(fh_session, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_session.pdf', monkey, experiment)))
% print(fh_timebin, '-dpdf', fullfile(FigDir, sprintf('neural_threshold_%s_%s_timebin.pdf', monkey, experiment)))




