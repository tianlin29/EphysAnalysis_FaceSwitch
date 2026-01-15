run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, threeExemplar, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'classifier'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'classifier'); mkdir(InterimDir);

%% classifier (category axis for main task)
if strcmp(experiment, 'faceColor')
    file_list = [1:8, 18:23];
else
    file_list = 1:n_files;
end
n_files = length(file_list);

REPEAT = 10;
for r = 1:REPEAT
    data = cell(n_files, 3); % (session, condition)
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

        % project pair 1 to pair 1
        ID = targ_cor; ID(task_set~=1) = NaN;
        opt.Kfold = 5;
        data{n,1} = PopulationClassifier(fnd, ID, opt);

        % project pair 2 to pair 2
        ID = targ_cor; ID(task_set~=2) = NaN;
        opt.Kfold = 5;
        data{n,2} = PopulationClassifier(fnd, ID, opt);

        % build pair 1 classifier
        ID = targ_cor; ID(task_set~=1) = NaN;
        opt.Kfold = 1;
        data_category_axis = PopulationClassifier(fnd, ID, opt);

        % project pair 2 to pair 1
        ID = targ_cor; ID(task_set~=2) = NaN;
        opt.Kfold = 1;
        opt.refBeta = data_category_axis.data_ind; % attach pair 1 classification result
        data{n,3} = PopulationClassifier(fnd, ID, opt);
    end
    save(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r)), 'data', 'opt')
end








