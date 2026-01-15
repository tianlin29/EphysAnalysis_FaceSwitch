run('../Initialize.m');
monkey = 'Nick';
experiment = 'learnTask2';
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, experiment, monkey); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, experiment, monkey); mkdir(InterimDir);

%% flow field
n = 1;
fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
fnd = fnd.extract_epoch(2);
fnd = fnd.extract_trial(fnd.getp('task_set')==2);

% get smooth detrended raster
r = fnd.raster(1); r = r{1}*1000; % (unit, time, trial)
[nunit, ntime, ntrial] = size(r);
r_mn = mean(r, 3);
r_detrend = r - r_mn;

kernel = fspecial('average', [1, 100]);
r_detrend = nanconv(r_detrend, kernel, 'same');

% PCA
targ_cho = fnd.getp('targ_cho');
ID = targ_cho;
ncond = max(ID(:));
psth = nan(nunit, ntime, ncond);
for c = 1:ncond
    psth(:,:,c) = mean(r_detrend(:,:,ID(1,:)==c), 3);
end

[coeff, score, latent] = pca(psth(:,:)');
npc = 2;
r_pc = reshape(coeff(:,1:npc)'*r_detrend(:,:), npc, ntime, ntrial);
psth_pc = reshape(score(:,1:npc)', npc, ntime, ncond);

% creat mesh
ngrid = 20;
x_edge = linspace(min(r_pc(:)), max(r_pc(:)), ngrid);
y_edge = linspace(min(r_pc(:)), max(r_pc(:)), ngrid);
[X, Y] = meshgrid(x_edge, y_edge);

% calculate velocity
[U_data, V_data] = deal(cell(ngrid, ngrid));
for i = 1:ngrid
    for j = 1:ngrid
        U_data{i,j} = [];
        V_data{i,j} = [];
    end
end

for t = 1:ntime-1
    t
    for tr = 1:ntrial
        x_t0 = r_pc(1,t,tr); x_t1 = r_pc(1,t+1,tr);
        y_t0 = r_pc(2,t,tr); y_t1 = r_pc(2,t+1,tr);
        col_idx = find(x_t0 >= x_edge(1:end-1) & x_t0 <= x_edge(2:end), 1);
        row_idx = find(y_t0 >= y_edge(1:end-1) & y_t0 <= y_edge(2:end), 1);

        u_t0 = x_t1 - x_t0;
        v_t0 = y_t1 - y_t0;
        U_data{row_idx, col_idx} = cat(1, U_data{row_idx, col_idx}, u_t0);
        V_data{row_idx, col_idx} = cat(1, V_data{row_idx, col_idx}, v_t0);
    end
end

[U, V] = deal(nan(ngrid, ngrid));
for i = 1:ngrid
    for j = 1:ngrid
        nsample = length(U_data{i,j});
        if nsample>=1000
            U(i,j) = mean(U_data{i,j});
            V(i,j) = mean(V_data{i,j});
        end
    end
end

% plot
opt.plot = set_plot_opt('vik', ncond);
figure; hold on
quiver(X, Y, U, V)
for c = 1:ncond
    plot(psth_pc(1,:,c), psth_pc(2,:,c), 'Color', opt.plot.color(c,:));
end












