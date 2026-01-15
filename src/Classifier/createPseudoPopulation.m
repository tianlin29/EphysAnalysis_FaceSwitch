function nfnd = createPseudoPopulation(fnd, condID, opt)

def.seed = NaN;
def.method = 'sortdur';
    % random: just randomly select trials from each unit
    % sortdur: prioritize trials with longer duration. Better if you want
    % to reduce NaN in data.
def.target_event = 1; % when sorting duration, you should specify event whose duration to be maximized
def.min_trial_per_unit = 0; % remove units with trial counts lower than this.
def.min_trial_dur = 0; % remove trials whose duration is shorter than this
def.preserve_param = {'targ_cho', 'targ_cor'};
    % these parameters will be kept in new fnd (other params won't be
    % copied).

opt = safeStructAssign(def, opt);

if isnan(opt.seed)
    opt.seed = round(rem(now,1)*1e5);
end
RandStream.setGlobalStream(RandStream('mt19937ar','Seed', opt.seed));

%% initialize

ntri = sum(~isnan(condID),2);
include_uni = ntri > opt.min_trial_per_unit;
if any(~include_uni)
    fprintf('[createPseudoPopulation] removing %d units due to low trial count.\n', sum(~include_uni));
    fnd = fnd.set_unit_criteria('append', 'on', 'custom', include_uni);
    condID = condID(include_uni,:);
end

% condID: cell x trial
fprintf('[createPseudoPopulation] extracting raster.\n');
raster = fnd.raster; % raster{1}: cell x time x trial

nwin = length(raster);
ncond = nanmax(condID(:));
ncell = size(raster{1},1);

tdur = squeeze(sum(~isnan(raster{opt.target_event}),2));

rmtr = tdur < opt.min_trial_dur;
if any(rmtr(:))
    fprintf('[createPseudoPopulation] Removing %d%% of trials due to short duration.\n', round(mean(rmtr(~isnan(condID(:))))*1e2));
    condID(rmtr) = NaN;
end

%% construct pseudo population raster

for n=1:length(opt.preserve_param)
    param.(opt.preserve_param{n}) = fnd.getp(opt.preserve_param{n}); % cell x trial
end

% create shuffle trial index
fprintf('[createPseudoPopulation] shuffling trial index.\n');
tr_ind = cell(ncond,1); % tr_ind{1} = cell x trial
mintri = nan(ncond,1);
for c=1:ncond
    mintri(c) = min(sum(condID == c,2));
    tr_ind{c} = nan(ncell, mintri(c));
    for u=1:ncell
        ind = find(condID(u,:) == c);
        switch opt.method
            case 'random'
                tr_ind{c}(u,:) = ind(randperm(mintri(c)));
            case 'sortdur'
                ind = ind(randperm(length(ind))); % shuffle once
                [~, ind2] = sort(tdur(u, ind), 'descend');
                tr_ind{c}(u,:) = ind(ind2(1:mintri(c)));
        end
    end
end
tr_ind = cell2mat(tr_ind'); % cell x trial
fprintf('[createPseudoPopulation] Total %d trials (average %d trials per condition)\n', sum(mintri), round(mean(mintri)));

% shuffle raster
fprintf('[createPseudoPopulation] creating pseudo population activity.\n');
for a=1:nwin
    for u=1:ncell
        raster{a}(u, :, 1:size(tr_ind,2)) = raster{a}(u, :, tr_ind(u,:));
    end
    raster{a}(:, :, (size(tr_ind,2)+1):end) = [];
end

% shuffle parameter
for n=1:length(opt.preserve_param)
    for u=1:ncell
        param.(opt.preserve_param{n})(u, 1:size(tr_ind,2)) = param.(opt.preserve_param{n})(u, tr_ind(u,:));
    end
    param.(opt.preserve_param{n})(:, (size(tr_ind,2)+1):end) = [];
end

%% create pseudo population fnd


fndopt.param = param;
fndopt.alignto = fnd.alignto;
fndopt.tstamp = fnd.tstamp;
fndopt.unitID = [ones(size(fnd.unitID,1), 2), (1:size(fnd.unitID,1))'];
fndopt.unit_type = fnd.unit_type;
fndopt.unit_quality = fnd.unit_quality;
fndopt.area = fnd.area;
fndopt.monkey = fnd.monkey;
fndopt.session_name = {'pseudo_session'};
fndopt.description = 'pseudo population generated using createPseudoPopulation.m';
fndopt.ntrial = ones(ncell,1) * sum(mintri);
nfnd = FND(raster, fndopt);

