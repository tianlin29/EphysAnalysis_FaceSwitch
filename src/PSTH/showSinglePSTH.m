function fh = showSinglePSTH(tstamp, psth, opt)

def.epoch = 2;
def.unitID = [];
def.plot = set_plot_opt('vik', 5);
opt = safeStructAssign(def, opt);

ncond = size(psth{opt.epoch}, 3);
nunit = size(opt.unitID, 1);
nfig = ceil(nunit/40);
fh = cell(nfig,1);
for u = 1:nunit
    fig_idx = ceil(u/40);
    subplot_idx = mod(u, 40);
    if subplot_idx==0; subplot_idx = 40; end
    if subplot_idx==1
        fh{fig_idx} = figure('Position', [40 60 650 800]);
    end
    subplot(8,5,subplot_idx); hold on
    for c = 1:ncond
        plot(tstamp{opt.epoch}, psth{opt.epoch}(u,:,c), 'Color', opt.plot.color(c,:));
    end
    title(sprintf('ch%d u%d', opt.unitID(u,2), opt.unitID(u,3)))
end

end