run('../Initialize.m');
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor, threeExemplar, passiveLong
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'regression'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'regression'); mkdir(InterimDir);

%%
t_mid = [100,200,300,400,500,600]; t_win = 200; n_time = length(t_mid);
n_repeat = 3;
[abs_morph] = deal(cell(n_files, 3)); pred_cor = cell(n_files, 3, n_time, n_repeat);
for n = 1:n_files
    n
    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    fnd = fnd.extract_epoch(2);
    fnd = fnd.extract_trial(fnd.getp('targ_cho')==fnd.getp('targ_cor'));

    targ_cor = fnd.getp('targ_cor'); targ_cor = targ_cor(1,:)';
    task_set = fnd.getp('task_set'); task_set = task_set(1,:)';
    morph_level = fnd.getp('morph_level'); morph_level = morph_level(1,:)'; abs_morph_level = abs(morph_level);
    I1 = task_set==1; I2 = task_set==2;

    abs_morph{n,1} = abs_morph_level(I1);
    abs_morph{n,2} = abs_morph_level(I2);
    abs_morph{n,3} = abs_morph_level(I2);

    for t = 1:n_time
        FR = fnd.FR({1, [t_mid(t)-t_win/2, t_mid(t)+t_win/2]})'; % (trial, unit)

        for r = 1:n_repeat
            cvmdl = fitclinear(FR(I1,:), targ_cor(I1), 'kFold', 5);
            y_hat = cvmdl.kfoldPredict;
            pred_cor{n,1,t,r} = y_hat==targ_cor(I1); % (trial, 1)
            
            cvmdl = fitclinear(FR(I2,:), targ_cor(I2), 'kFold', 5);
            y_hat = cvmdl.kfoldPredict;
            pred_cor{n,2,t,r} = y_hat==targ_cor(I2); % (trial, 1)
            
            mdl = fitclinear(FR(I1,:), targ_cor(I1));
            y_hat = predict(mdl, FR(I2,:));
            pred_cor{n,3,t,r} = y_hat==targ_cor(I2); % (trial, 1)            
        end
    end
end

fcoh = linspace(0.03, 1, 100)'; % coh for function
[lucoh] = deal(cell(n_files, 3)); [threshold, accuracy] = deal(nan(n_files, 3, n_time, n_repeat));
[p, pse, fresp] = deal(cell(n_files, 3, n_time));
for n = 1:n_files
    for c = 1:3
        ucoh = unique(abs_morph{n,c});
        lucoh{n,c} = log(ucoh); % log coh
        mincoh = ucoh(1);
        if mincoh == 0
            mincoh = ucoh(2)/ucoh(3) * ucoh(2);
            lucoh{n,c}(1) = log(mincoh);
        end

        for t = 1:n_time
            for r = 1:n_repeat
                [p{n,c,t,r}, pse{n,c,t,r}] = calcGroupMean(pred_cor{n,c,t,r}, abs_morph{n,c}, unique(abs_morph{n,c}), 'binary');
                [alpha, ~, stat_glmfit] = glmfit([abs_morph{n,c}], pred_cor{n,c,t,r}, 'binomial', 'link', 'logit', 'Constant', 'off'); % logit(P) = a2*coh
                fresp{n,c,t,r} = glmval(alpha, fcoh, 'logit', 'Constant', 'off');
                threshold(n,c,t,r) = fminsearch(@(abs_coh) abs(glmval(alpha, abs_coh, 'logit', 'Constant', 'off') - 0.816), .2);
                accuracy(n,c,t,r) = mean(fresp{n,c,t,r});
            end
        end
    end
end

% plot psychometric function
t = 4;
clear opt
opt.color = [0 0 0; 44 145 224; 255 0 0]/255;
figure('Position', [100 100 1000 800]);
for n = 1:n_files
    subplot(5,5,n); hold on
    for c = 1:3
        plot(log(fcoh), fresp{n,c,t}, 'Color', opt.color(c,:));
        plot(lucoh{n,c}, p{n,c,t}, '.', 'markers', 7, 'Color', opt.color(c,:));
        cerrorbar(lucoh{n,c}, p{n,c,t}, pse{n,c,t}, 'Color', opt.color(c,:));
    end
    format_panel(gca, ...
        'ylim', [.5 1], 'ytick', .5:.125:1, ...
        'xlim', log([3 100]/100), 'xtick', log([3 6 12 24 48 96]'/100), 'xticklabel', {'0', '6', '', '24', '', '96'}, ...
        'xlabel', {'Stimulus strength','(% Morph)'}, 'ylabel', 'P(correct)');
end

% plot threshold
threshold_mn = mean(threshold, 4);
threshold_se = std(threshold, [], 4) / sqrt(n_repeat);
threshold_nrml = threshold_mn ./ threshold_mn(:,1,:);
figure;
for t = 1:n_time
    subplot(1,n_time,t); hold on
    for c = 1:3
        plot(1:n_files, threshold_nrml(:,c,t), '.-', 'Color', opt.color(c,:));
        cerrorbar(1:n_files, threshold_nrml(:,c,t), threshold_se(:,c,t), 'Color', opt.color(c,:));
    end
end
format_panel(gcf, 'ylim', [0 10])

% plot accuracy
accuracy_nrml = accuracy ./ accuracy(:,1,:);
figure;
for t = 1:n_time
    subplot(1,n_time,t); hold on
    for c = 1:3
        plot(1:n_files, accuracy(:,c,t), '.-', 'Color', opt.color(c,:));
    end
end
format_panel(gcf, 'ylim', [0.5 1])




