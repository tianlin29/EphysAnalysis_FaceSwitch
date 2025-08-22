function fh = showPopPSTH_choice_dependent(tstamp, psth1, psth2, opt)

def.cutoff = [];
def.plot = set_plot_opt('vik', 10);
def.epoch = 1;
def.title = [];
opt = safeStructAssign(def, opt);

p1 = squeeze(nanmean(psth1{opt.epoch},1))';
p2 = squeeze(nanmean(psth2{opt.epoch},1))';

ncond = size(p1,1);

fh = figure('color', 'w', 'pos', [100 100 800 400]);

ax = nan(ncond,1);
for c=1:ncond
    ax(c) = subplot(2, ceil(ncond/2), c);
    hold on;
    h1 = plot(tstamp{opt.epoch}, p1(c,:), 'color', opt.plot.color(1,:));
    h2 = plot(tstamp{opt.epoch}, p2(c,:), 'color', opt.plot.color(2,:));
    if c==ncond
        legend([h1, h2], {'Choice 1', 'Choice 2'});
    end
    if ~isempty(opt.title)
        title(opt.title{c});
    end
end

format_panel(ax, 'xlabel', 'Time (ms)', 'ylabel', 'Firing rate (sp/s)');


end
