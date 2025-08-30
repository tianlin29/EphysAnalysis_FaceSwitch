function fh = comparePCA(tstamp1, tstamp2, data1, data2, opt)

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
def.plot1 = set_plot_opt('red', 10);
def.plot2 = set_plot_opt('blue', 10);
def.PC_target = 'average';
 % which data to be used for calculating PC coefficient
 % data1 : data 1 only
 % data2 : data 2 only
 % average: average PSTHs of data 1 and data 2 (assuming the cells and conditions are matched)
 % concatenate: concatenate PSTHs of data 1 and data 2

opt = safeStructAssign(def, opt);

PSTH1 = data1{opt.epoch};
PSTH2 = data2{opt.epoch};
Tstamp1 = tstamp1{opt.epoch};
Tstamp2 = tstamp2{opt.epoch};

if size(PSTH1,1) ~= size(PSTH2,1)
    error('Number of units does not match between data1 and data2');
end

%% preprocess
% mean imputation to nan
if sum(isnan(PSTH1(:)))>0
    fill_val = bsxfun(@times, isnan(PSTH1), nanmean(PSTH1(:,:),2));
    PSTH1(isnan(PSTH1(:))) = fill_val(isnan(PSTH1(:)));
end
if sum(isnan(PSTH2(:)))>0
    fill_val = bsxfun(@times, isnan(PSTH2), nanmean(PSTH2(:,:),2));
    PSTH2(isnan(PSTH2(:))) = fill_val(isnan(PSTH2(:)));
end

% detrend
PSTH1 = bsxfun(@minus, PSTH1, nanmean(PSTH1,3));
PSTH2 = bsxfun(@minus, PSTH2, nanmean(PSTH2,3));

% PC range

switch opt.PC_target
    case 'data1'
        ind = Tstamp1 >= opt.PC_range(1) & Tstamp1 <= opt.PC_range(2);
        pcPSTH = PSTH1(:, ind, :); % cell x time x cond
        pcPSTH = pcPSTH(:,:)';
    case 'data2'
        ind = Tstamp2 >= opt.PC_range(1) & Tstamp2 <= opt.PC_range(2);
        pcPSTH = PSTH2(:, ind, :); % cell x time x cond
        pcPSTH = pcPSTH(:,:)';
    case 'average'
        ind = Tstamp1 >= opt.PC_range(1) & Tstamp1 <= opt.PC_range(2);
        pcPSTH1 = PSTH1(:, ind, :); % cell x time x cond
        ind = Tstamp2 >= opt.PC_range(1) & Tstamp2 <= opt.PC_range(2);
        pcPSTH2 = PSTH2(:, ind, :); % cell x time x cond
        pcPSTH = (pcPSTH1(:,:)' + pcPSTH2(:,:)') / 2;
    case 'concatenate'
        ind = Tstamp1 >= opt.PC_range(1) & Tstamp1 <= opt.PC_range(2);
        pcPSTH1 = PSTH1(:, ind, :); % cell x time x cond
        ind = Tstamp2 >= opt.PC_range(1) & Tstamp2 <= opt.PC_range(2);
        pcPSTH2 = PSTH2(:, ind, :); % cell x time x cond
        pcPSTH = [pcPSTH1(:,:)'; pcPSTH2(:,:)'];
    otherwise
        error('No such PC target: %s', opt.PC_target);
end
%% run PCA
coef = pca(pcPSTH);
nTime = length(opt.Time);
PCscore1 = cell(nTime,1);
PCscore2 = cell(nTime,1);
for t=1:nTime
    ind = opt.Time(t) == Tstamp1;
    if ~any(ind)
        error('No time %d exists in data 1', opt.Time(t));
    end
    PCscore1{t} = permute(PSTH1(:,ind,:), [1 3 2])' * coef(:,1:3);
    PCscore1{t} = PCscore1{t} .* (ones(size(PCscore1{t},1),1) * opt.dim_sign);
    ind = opt.Time(t) == Tstamp2;
    if ~any(ind)
        error('No time %d exists in data 2', opt.Time(t));
    end
    PCscore2{t} = permute(PSTH2(:,ind,:), [1 3 2])' * coef(:,1:3);
    PCscore2{t} = PCscore2{t} .* (ones(size(PCscore2{t},1),1) * opt.dim_sign);
end

%% curve fit
cx = linspace(min(opt.mean_coh), max(opt.mean_coh), 100);
pt1 = cell(nTime,1);
pt2 = cell(nTime,1);
for t=1:nTime
    pt1{t} = nan(100, 3);
    pt2{t} = nan(100, 3);
    for d=1:3
        [x, idx] = sort(opt.mean_coh(:));
        y = PCscore1{t}(idx,d);
        Param = csaps(x, y, [opt.roughness, opt.roughness_weighting], [], opt.point_weighting);
        pt1{t}(:,d) = fnval(Param, cx);
        
        y = PCscore2{t}(idx,d);
        Param = csaps(x, y, [opt.roughness, opt.roughness_weighting], [], opt.point_weighting);
        pt2{t}(:,d) = fnval(Param, cx);
    end
end

%% plot
col1 = opt.plot1.color;
fcol1 = opt.plot1.facecolor;
col2 = opt.plot2.color;
fcol2 = opt.plot2.facecolor;

fh = figure('color', 'w', 'position', [100 100 200 * nTime 200], 'paperpositionmode', 'auto');
for t=1:nTime
    subplot(1, nTime, t);
    hold on;
    d = opt.dim;
    plot3(pt1{t}(:,d(1)), pt1{t}(:,d(2)), pt1{t}(:,d(3)), 'col', col1(end,:), 'linew', 1);
    for n=1:size(PCscore1{t},1)
        plot3(PCscore1{t}(n,d(1)), PCscore1{t}(n,d(2)), PCscore1{t}(n,d(3)), 'o', ...
            'markeredgecol', col1(n,:), 'markerfacecol', fcol1(n,:), 'markers', 5);
    end
    
    plot3(pt2{t}(:,d(1)), pt2{t}(:,d(2)), pt2{t}(:,d(3)), 'col', col2(end,:), 'linew', 1);
    for n=1:size(PCscore2{t},1)
        plot3(PCscore2{t}(n,d(1)), PCscore2{t}(n,d(2)), PCscore2{t}(n,d(3)), 'o', ...
            'markeredgecol', col2(n,:), 'markerfacecol', fcol2(n,:), 'markers', 5);
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


