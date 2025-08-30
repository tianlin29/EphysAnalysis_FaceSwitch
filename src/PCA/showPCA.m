function fh = showPCA(tstamp, data, opt)

def.PC_range = [250 600];
def.Time = 400;
def.roughness = 1e-3;
def.point_weighting = [1 1 1 1 1 1 1 1 1 1];
def.roughness_weighting = [1 1 1 .1 .01 .1 1 1 1];
def.dim_sign = [1 1 1];
def.dim = [1 2 3];
def.axis_lim = {};
def.view = [];
def.epoch = 1;
def.mean_coh = [];
def.plot = set_plot_opt('vik', 10);

opt = safeStructAssign(def, opt);

PSTH = data{opt.epoch};
Tstamp = tstamp{opt.epoch};

%% preprocess
% mean imputation to nan
if sum(isnan(PSTH(:)))>0
    fill_val = bsxfun(@times, isnan(PSTH), nanmean(PSTH(:,:),2));
    PSTH(isnan(PSTH(:))) = fill_val(isnan(PSTH(:)));
end

% detrend
PSTH = bsxfun(@minus, PSTH, nanmean(PSTH,3));

% PC range
ind = Tstamp >= opt.PC_range(1) & Tstamp <= opt.PC_range(2);
pcPSTH = PSTH(:, ind, :); % cell x time x cond

%% run PCA
coef = pca(pcPSTH(:,:)');
nTime = length(opt.Time);
PCscore = cell(nTime,1);
for t=1:nTime
    ind = opt.Time(t) == Tstamp;
    PCscore{t} = permute(PSTH(:,ind,:), [1 3 2])' * coef(:,1:3);
    PCscore{t} = PCscore{t} .* (ones(size(PCscore{t},1),1) * opt.dim_sign);
end

%% curve fit
cx = linspace(min(opt.mean_coh), max(opt.mean_coh), 100);
pt = cell(nTime,1);
for t=1:nTime
    pt{t} = nan(100, 3);
    for d=1:3
        [x, idx] = sort(opt.mean_coh(:));
        y = PCscore{t}(idx,d);
        Param = csaps(x, y, [opt.roughness, opt.roughness_weighting], [], opt.point_weighting);
        pt{t}(:,d) = fnval(Param, cx);
    end
end

%% plot
col = opt.plot.color;
fcol = opt.plot.facecolor;

fh = figure('color', 'w', 'position', [100 100 200 * nTime 200], 'paperpositionmode', 'auto');
for t=1:nTime
    subplot(1, nTime, t);
    hold on;
    d = opt.dim;
    plot3(pt{t}(:,d(1)), pt{t}(:,d(2)), pt{t}(:,d(3)), 'col', [.5 .5 .5], 'linew', 1);
    for n=1:size(PCscore{t},1)
        plot3(PCscore{t}(n,d(1)), PCscore{t}(n,d(2)), PCscore{t}(n,d(3)), 'o', ...
            'markeredgecol', col(n,:), 'markerfacecol', fcol(n,:), 'markers', 5);
    end
    axis square;
    xlabel(sprintf('PC %d', d(1)));
    ylabel(sprintf('PC %d', d(2)));
    zlabel(sprintf('PC %d', d(3)));
    if ~isempty(opt.axis_lim)
        xlim(opt.axis_lim{1});
        ylim(opt.axis_lim{2});
        zlim(opt.axis_lim{3});
    end
    if ~isempty(opt.view)
        view(opt.view);
    end
    title(sprintf('%d ms', opt.Time(t)));
end


