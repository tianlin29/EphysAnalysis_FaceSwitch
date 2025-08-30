function fh = show_single_tr_trajectory(Tstamp, raster, ID, coef, detPSTH, opt)

def.t_range = [0 600];
def.dim = [1 2 3];
def.conv_kernel = [];


def.xlim = [];
def.ylim = [];
def.zlim = [];
def.view_angle = [48, 28];
def.color = [1 0 0; 0 0 1];

opt = safeStructAssign(def, opt);

ndim = 3;


%% smooth raster
raster = raster * 1e3; % cell x time x trial

% detrend
raster = raster - detPSTH;

psth = nan(size(raster,1), size(raster, 2), nanmax(ID(:)));
for id = 1:nanmax(ID(:))
    psth(:, :, id) = nanmean(raster(:, :, ID(1,:) == id),3);
end


%% conver to PCA

T = opt.t_range(1):opt.t_range(2);
ind = Tstamp >= opt.t_range(1) & Tstamp <= opt.t_range(2);
Tstamp = Tstamp(ind);

nT = length(T);

ntri = size(raster, 3);
ncond = size(psth, 3);

score = nan(nT, ntri, ndim);
psth_score = nan(nT, ncond, ndim);
for dim=1:ndim
    PSTH1 = raster(:, ind, :);
    sc = PSTH1(:,:)' * coef(:, dim);
    score(:,:,dim) = reshape(sc, [nT ntri]);
    
    PSTH1 = psth(:, ind, :);
    sc = PSTH1(:,:)' * coef(:, dim);
    psth_score(:,:,dim) = reshape(sc, [nT ncond]);
end

if ~isempty(opt.conv_kernel)
    for d=1:ndim
        for n=1:ntri
            score(:, n, d) = nanconv(score(:, n, d), opt.conv_kernel(:), 'same');
        end
        for n=1:ncond
            psth_score(:, n, d) = nanconv(psth_score(:, n, d), opt.conv_kernel(:), 'same');
        end
    end
end

%% Plot figure
fh = figure('Color', 'w', 'Position', [100 100 500 500], 'PaperPositionMode', 'auto', 'Units', 'normalized');
hold on;

minv = squeeze(min(min(score,[],1),[],2));
maxv = squeeze(max(max(score,[],1),[],2));
if isempty(opt.ylim)
    opt.ylim = [-.1 1.1] * (maxv(1) - minv(1)) + minv(1);
end
if isempty(opt.zlim)
    opt.zlim = [-.1 1.1] * (maxv(2) - minv(2)) + minv(2);
end

xlim(Tstamp([1 end]));
ylim(opt.ylim);
zlim(opt.zlim);
xlabel('Time');
ylabel('PC 1');
zlabel('PC 2');
view(opt.view_angle);
grid on;

r = randsample(find(~isnan(ID(1,:))), 200);
for c1=r(:)'
    plot3(Tstamp, score(:,c1, 1), score(:,c1, 2), 'color', opt.color(ID(1, c1), :) * .5 + [.5 .5 .5], 'linestyle', '-', 'linew', 1);
end

for c = 1:ncond
    plot3(Tstamp, psth_score(:, c, 1), psth_score(:, c, 2), 'color', opt.color(c, :), 'linestyle', '-', 'linew', 3);
end



end


