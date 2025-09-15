function fh = comparePCA_trajectory(fnd1, fnd2, condID1, condID2, opt)

def.epoch = 1;
def.PC_range = [250 600];
def.t_range = [0 600];
def.anchor_t = 0:50:600;
def.interval = 10;
def.dim = [1 2 3];
def.PSTH_conv = {'boxcar', 100};
def.PC_kernel = {'bartlett', 200};
def.interp_method = {'polyfit',3};
def.PC_target = 'average';
 % which data to be used for calculating PC coefficient
 % data1 : data 1 only
 % data2 : data 2 only
 % average: average PSTHs of data 1 and data 2 (assuming the cells and conditions are matched)
 % concatenate: concatenate PSTHs of data 1 and data 2

def.xlim = [];
def.ylim = [];
def.zlim = [];
def.view_angle = [48, 28];
def.plot1 = [];
def.plot2 = [];

opt = safeStructAssign(def, opt);

ndim = 3;

Tstamp1 = fnd1.tstamp{opt.epoch};
Tstamp2 = fnd2.tstamp{opt.epoch};

PSTH1 = fnd1.PSTH(condID1, opt.PSTH_conv);
PSTH1 = PSTH1{opt.epoch};
PSTH2 = fnd2.PSTH(condID2, opt.PSTH_conv);
PSTH2 = PSTH2{opt.epoch};

PSTH1_unsmooth = fnd1.PSTH(condID1, []);
PSTH1_unsmooth = PSTH1_unsmooth{opt.epoch};
PSTH2_unsmooth = fnd2.PSTH(condID2, []);
PSTH2_unsmooth = PSTH2_unsmooth{opt.epoch};


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
PSTH1_unsmooth = bsxfun(@minus, PSTH1_unsmooth, nanmean(PSTH1_unsmooth,3));
PSTH2_unsmooth = bsxfun(@minus, PSTH2_unsmooth, nanmean(PSTH2_unsmooth,3));

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
coef = pca(pcPSTH(:,:));

T = opt.t_range(1):opt.t_range(2);

nT = length(T);

ncond1 = size(PSTH1, 3);
ncond2 = size(PSTH2, 3);

score1 = nan(nT, ncond1, ndim);
score2 = nan(nT, ncond2, ndim);
for dim=1:ndim
    ind = Tstamp1 >= opt.t_range(1) & Tstamp1 <= opt.t_range(2);
    PSTHtmp = PSTH1_unsmooth(:, ind, :);
    sc = PSTHtmp(:,:)' * coef(:, dim);
    score1(:,:,dim) = reshape(sc, [nT ncond1]);
    
    ind = Tstamp2 >= opt.t_range(1) & Tstamp2 <= opt.t_range(2);
    PSTHtmp = PSTH2_unsmooth(:, ind, :);
    sc = PSTHtmp(:,:)' * coef(:, dim);
    score2(:,:,dim) = reshape(sc, [nT ncond2]);
end

if ~isempty(opt.PC_kernel)
    for d=1:ndim
        for n=1:ncond1
            score1(:, n, d) = nanconv(score1(:, n, d), smoothing_filter(opt.PC_kernel)', 'same');
        end
        for n=1:ncond2
            score2(:, n, d) = nanconv(score2(:, n, d), smoothing_filter(opt.PC_kernel)', 'same');
        end
    end
end

%% Plot figure
plt1 = opt.plot1;
plt2 = opt.plot2;
fh = figure('Color', 'w', 'Position', [100 100 500 500], 'PaperPositionMode', 'auto', 'Units', 'normalized');
hold on;

minv = squeeze(min(min(cat(2, score1, score2),[],1),[],2));
maxv = squeeze(max(max(cat(2, score1, score2),[],1),[],2));
if isempty(opt.xlim)
    opt.xlim = [-.1 1.1] * (maxv(1) - minv(1)) + minv(1);
end
if isempty(opt.ylim)
    opt.ylim = [-.1 1.1] * (maxv(2) - minv(2)) + minv(2);
end
if isempty(opt.zlim)
    opt.zlim = [-.1 1.1] * (maxv(3) - minv(3)) + minv(3);
end
xlim(opt.xlim);
ylim(opt.ylim);
zlim(opt.zlim);
xlabel(sprintf('PC %d', opt.dim(1)));
ylabel(sprintf('PC %d', opt.dim(2)));
zlabel(sprintf('PC %d', opt.dim(3)));
view(opt.view_angle);
grid on;

trjh = [];
show_trj;
uicontrol(gcf, 'style', 'pushbutton', 'string', 'show movie', 'Units', 'normalized', ...
    'position', [.77 .85 .2 .04], 'callback', @show_movie);
uicontrol(gcf, 'style', 'text', 'string', 'Click dot to show time and manifold', 'Units', 'normalized', ...
    'position', [.02 .96 .5 .04], 'backgroundcolor', [1 1 1]);

function show_trj
    for c1=1:ncond1
        h = plot3(score1(:,c1, 1), score1(:,c1, 2), score1(:,c1, 3), ...
            'color', plt1.color(c1,:), 'linestyle', plt1.linestyle{c1}, 'linew', plt1.linewidth(c1));
        set(h, 'hittest', 'off');
        trjh = [trjh; h]; %#ok<AGROW>
        for m1=1:length(opt.anchor_t)
            ti = find(T == opt.anchor_t(m1));
            h = plot3(score1(ti,c1, 1), score1(ti,c1, 2), score1(ti,c1, 3), 'o', ...
                'color', plt1.color(c1,:), 'markerfacecolor', plt1.facecolor(c1,:), 'markersize', plt1.markersize(c1));
            trjh = [trjh; h]; %#ok<AGROW>
        end
    end
    for c1=1:ncond2
        h = plot3(score2(:,c1, 1), score2(:,c1, 2), score2(:,c1, 3), ...
            'color', plt2.color(c1,:), 'linestyle', plt2.linestyle{c1}, 'linew', plt2.linewidth(c1));
        set(h, 'hittest', 'off');
        trjh = [trjh; h]; %#ok<AGROW>
        for m1=1:length(opt.anchor_t)
            ti = find(T == opt.anchor_t(m1));
            h = plot3(score2(ti,c1, 1), score2(ti,c1, 2), score2(ti,c1, 3), 'o', ...
                'color', plt2.color(c1,:), 'markerfacecolor', plt2.facecolor(c1,:), 'markersize', plt2.markersize(c1));
            trjh = [trjh; h]; %#ok<AGROW>
        end
    end
end

function show_movie(~,~)
    delete_interp_manifold;
    for n1=1:length(trjh)
        set(trjh(n1), 'visible', 'off');
    end
    for t1=1:opt.interval:length(T)
        h1 = [];
        for c1=1:ncond1
            h2 = plot3(score1(1:t1,c1, 1), score1(1:t1,c1, 2), score1(1:t1,c1, 3), ...
                'color', plt1.color(c1,:), 'linestyle', plt1.linestyle{c1}, 'linew', plt1.linewidth(c1));
            h1 = [h1;h2]; %#ok<AGROW>
        end
        for c1=1:ncond2
            h2 = plot3(score2(1:t1,c1, 1), score2(1:t1,c1, 2), score2(1:t1,c1, 3), ...
                'color', plt2.color(c1,:), 'linestyle', plt2.linestyle{c1}, 'linew', plt2.linewidth(c1));
            h1 = [h1;h2]; %#ok<AGROW>
        end
        title(sprintf('%d ms', T(t1)), 'fontsize', 12);
        drawnow;
        pause(.1);
        for n1=1:length(h1)
            set(h1(n1), 'visible', 'off');
        end
    end
    trjh = [];
    show_trj;
end


end


