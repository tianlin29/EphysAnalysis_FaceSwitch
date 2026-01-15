function fh = showPCA_trajectory(fnd, condID, opt)

def.epoch = 1;
def.PC_range = [250 600];
def.t_range = [0 600];
def.anchor_t = 0:50:600;
def.interval = 10;
def.dim = [1 2 3];
def.PSTH_conv = {'boxcar', 100};
def.PC_kernel = {'bartlett', 200};
def.interp_method = {'polyfit',3};
def.dash_line_3d = false;

def.xlim = [];
def.ylim = [];
def.zlim = [];
def.view_angle = [48, 28];
def.plot = [];

def.PSTH = [];
def.PSTH_unsmooth = [];

opt = safeStructAssign(def, opt);

ndim = 3;

if ~isempty(opt.PSTH)
    PSTH = opt.PSTH;
    PSTH_unsmooth = opt.PSTH_unsmooth;
else
    PSTH = fnd.PSTH(condID, opt.PSTH_conv);
    PSTH = PSTH{opt.epoch};

    PSTH_unsmooth = fnd.PSTH(condID, []);
    PSTH_unsmooth = PSTH_unsmooth{opt.epoch};
end

Tstamp = fnd.tstamp{opt.epoch};


%% preprocess
% mean imputation to nan
if sum(isnan(PSTH(:)))>0
    fill_val = bsxfun(@times, isnan(PSTH), nanmean(PSTH(:,:),2));
    PSTH(isnan(PSTH(:))) = fill_val(isnan(PSTH(:)));
end

% detrend
PSTH = bsxfun(@minus, PSTH, nanmean(PSTH,3));
PSTH_unsmooth = bsxfun(@minus, PSTH_unsmooth, nanmean(PSTH_unsmooth,3));

% PC range
ind = Tstamp >= opt.PC_range(1) & Tstamp <= opt.PC_range(2);
pcPSTH = PSTH(:, ind, :); % cell x time x cond

%% run PCA
coef = pca(pcPSTH(:,:)');

T = opt.t_range(1):opt.t_range(2);
ind = Tstamp >= opt.t_range(1) & Tstamp <= opt.t_range(2);

nT = length(T);

ncond = size(PSTH, 3);

score = nan(nT, ncond, ndim);
for dim=1:ndim
    PSTH1 = PSTH_unsmooth(:, ind, :);
    sc = PSTH1(:,:)' * coef(:, dim);
    score(:,:,dim) = reshape(sc, [nT ncond]);
end

if ~isempty(opt.PC_kernel)
    for n=1:ncond
        for d=1:ndim
            score(:, n, d) = nanconv(score(:, n, d), smoothing_filter(opt.PC_kernel)', 'same');
        end
    end
end

%% Curve fitting
nAT = length(opt.anchor_t);
method = opt.interp_method;
switch method{1}
    case 'interparc'
        points = nan(nAT, ndim, 100);
        wb = waitbar_text(0);
        for m=1:nAT
            waitbar_text(m/nAT, wb);
            ti = T == opt.anchor_t(m);
            for dim=1:ndim
                points(m, :,:) = interparc(100, score(ti, :,1), score(ti, :,2), ...
                    score(ti, :,3), 'spline')';
            end
        end
        waitbar_text('close', wb);
    case 'polyfit'
        Beta = cell(nAT, ndim);
        for m=1:nAT
            ti = T == opt.anchor_t(m);
            for dim=1:ndim
                y = score(ti, :,dim);
                x = linspace(0, 1, length(y));
                Beta{m, dim} = polyfit(x(:), y(:), method{2});
            end
        end
        x = linspace(0, 1, 100);
        points = nan(nAT, ndim, length(x));
        for m=1:nAT
            for dim=1:ndim
                points(m, dim,:) = polyval(Beta{m, dim}, x);
            end
        end
    otherwise
        error('unknown n-d interpolation method');
end

%% Plot figure
plt = opt.plot;
fh = figure('Color', 'w', 'Position', [100 100 500 500], 'PaperPositionMode', 'auto', 'Units', 'normalized');
hold on;

minv = squeeze(min(min(score,[],1),[],2));
maxv = squeeze(max(max(score,[],1),[],2));
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
manh = [];
uicontrol(gcf, 'style', 'pushbutton', 'string', 'draw manifold', 'Units', 'normalized', ...
    'position', [.77 .95 .2 .04], 'callback', @draw_interp_manifold);
uicontrol(gcf, 'style', 'pushbutton', 'string', 'delete manifold', 'Units', 'normalized', ...
    'position', [.77 .9 .2 .04], 'callback', @delete_interp_manifold);
uicontrol(gcf, 'style', 'pushbutton', 'string', 'show movie', 'Units', 'normalized', ...
    'position', [.77 .85 .2 .04], 'callback', @show_movie);
uicontrol(gcf, 'style', 'text', 'string', 'Click dot to show time and manifold', 'Units', 'normalized', ...
    'position', [.02 .96 .5 .04], 'backgroundcolor', [1 1 1]);

function show_trj
    for c1=1:ncond
        if opt.dash_line_3d && strcmp(plt.linestyle{c1}, '--')
            warning('A wierd way to plot dashed line in 3d space.')
            hold on
            for i = 1:size(score,1)-1
                if mod(i, 16) == 1 % plot solid line
                    plot3(score(i:i+8,c1, 1), score(i:i+8,c1, 2), score(i:i+8,c1, 3), ...
                        'color', plt.color(c1,:), 'linestyle', '-', 'linew', plt.linewidth(c1));
                else
                    % do not plot solid line. The gap makes it a dashed line.
                end
            end
        else
            h = plot3(score(:,c1, 1), score(:,c1, 2), score(:,c1, 3), ...
                'color', plt.color(c1,:), 'linestyle', plt.linestyle{c1}, 'linew', plt.linewidth(c1));
            set(h, 'hittest', 'off');
        end
        trjh = [trjh; h]; %#ok<AGROW>
        for m1=1:length(opt.anchor_t)
            ti = find(T == opt.anchor_t(m1));
            h = plot3(score(ti,c1, 1), score(ti,c1, 2), score(ti,c1, 3), 'o', ...
                'color', plt.color(c1,:), 'markerfacecolor', plt.facecolor(c1,:), 'markersize', plt.markersize(c1), ...
                'ButtonDownFcn', {@show_interp_manifold, m1});
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
        for c1=1:ncond
            h2 = plot3(score(1:t1,c1, 1), score(1:t1,c1, 2), score(1:t1,c1, 3), ...
                'color', plt.color(c1,:), 'linestyle', plt.linestyle{c1}, 'linew', plt.linewidth(c1));
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

function draw_interp_manifold(~,~)
    for m1=1:length(opt.anchor_t)
        h1 = plot3(squeeze(points(m1, 1, :)), squeeze(points(m1, 2, :)), squeeze(points(m1, 3, :)), ...
            'linew', 2, 'color', [.5 .5 .5]);
        manh = [manh; h1]; %#ok<AGROW>
    end
end
function delete_interp_manifold(~,~)
    for n1=1:length(manh)
        set(manh(n1), 'visible', 'off');
    end
    manh = [];
end

function show_interp_manifold(~,~,m)
    h1 = plot3(squeeze(points(m, 1, :)), squeeze(points(m, 2, :)), squeeze(points(m, 3, :)), ...
        'linew', 2, 'color', [.5 .5 .5]);
    manh = [manh; h1];
    title(sprintf('%d ms', opt.anchor_t(m)));
end


end


