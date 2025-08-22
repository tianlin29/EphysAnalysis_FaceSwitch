function fh = showPopPSTH(tstamp, psth, opt)

def.epoch = []; % e.g., [2,3] means the second and the third epoch
def.cutoff = [];
def.plot = set_plot_opt('vik', 10);
def.flip_preference = [];
def.legend = [];
def.legend_pos = [];
def.less_timepoints = false; % plot less timepoints to make the lines look smoother
def.psth_se = [];
def.event_label = {};
opt = safeStructAssign(def, opt);

% specify psth_se
psth_se = opt.psth_se;

% specify epoch
if ~isempty(opt.epoch)
    psth = psth(opt.epoch);
    if ~isempty(psth_se); psth_se = psth_se(opt.epoch); end
    tstamp = tstamp(opt.epoch);
    opt.cutoff = opt.cutoff(opt.epoch);
end

% flip preference
if ~isempty(opt.flip_preference)
    if length(opt.flip_preference)~=size(psth{1},1)
        error('The number of units does not match: flip_preference');
    end
    for n = 1:length(psth)
        for m = 1:size(psth{n},1)
            if opt.flip_preference(m)==1
                psth{n}(m,:,:) = flip(psth{n}(m,:,:), 3);
                if ~isempty(psth_se); psth_se{n}(m,:,:) = flip(psth_se{n}(m,:,:), 3); end
            end
        end
    end
end

% cutoff
for n = 1:length(psth)
    if size(psth{n},3)>1 % more than one condition
        psth{n} = squeeze(nanmean(psth{n},1))'; % (cond, time)
        if ~isempty(psth_se); psth_se{n} = squeeze(nanmean(psth_se{n},1))'; end
    else % only one condition
        psth{n} = squeeze(nanmean(psth{n},1)); % (cond, time)
        if ~isempty(psth_se); psth_se{n} = squeeze(nanmean(psth_se{n},1)); end
    end
    if ~isempty(opt.cutoff)
        if size(opt.cutoff{n},1)==1 % same cutoff across conditions
            psth{n} = psth{n}(:, opt.cutoff{n}(1):opt.cutoff{n}(2));
            if ~isempty(psth_se); psth_se{n} = psth_se{n}(:, opt.cutoff{n}(1):opt.cutoff{n}(2)); end
            tstamp{n} = tstamp{n}(opt.cutoff{n}(1):opt.cutoff{n}(2));
        else % different cutoff across conditions
            for c = 1:size(opt.cutoff{n},1)
                if ~all(isnan(opt.cutoff{n}(c,:)))
                    psth{n}(c, 1:(opt.cutoff{n}(c,1)-1)) = NaN;
                    psth{n}(c, (opt.cutoff{n}(c,2)+1):end) = NaN;
                    if ~isempty(psth_se)
                        psth_se{n}(c, 1:(opt.cutoff{n}(c,1)-1)) = NaN;
                        psth_se{n}(c, (opt.cutoff{n}(c,2)+1):end) = NaN;
                    end
                end
            end
            mn = min(opt.cutoff{n}(:,1));
            mx = max(opt.cutoff{n}(:,2));
            psth{n} = psth{n}(:, mn:mx);
            if ~isempty(psth_se); psth_se{n} = psth_se{n}(:, mn:mx); end
            tstamp{n} = tstamp{n}(mn:mx);
        end
    end
end

% smooth the data by using less timepoints
if opt.less_timepoints
    idx = 16;
    for e = 1:length(tstamp) % epoch
        n_tp = length(tstamp{e});
        psth{e} = interp1(tstamp{e}', psth{e}', linspace(tstamp{e}(1), tstamp{e}(end), round(n_tp/idx))')';
        if ~isempty(psth_se); psth_se{e} = interp1(tstamp{e}', psth_se{e}', linspace(tstamp{e}(1), tstamp{e}(end), round(n_tp/idx))')'; end
        tstamp{e} = linspace(tstamp{e}(1), tstamp{e}(end), round(n_tp/idx));
    end
end

fh = figure('color', 'w', 'pos', [100 100 200 * length(psth), 200]);
ax = plot_multiepoch_trace(tstamp, psth, psth_se, opt);
axes(ax);
xlabel('Time (ms)');
ylabel('Firing rate (sp/s)');

end
