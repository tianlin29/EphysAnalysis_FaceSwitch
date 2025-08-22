function fh = show_popresp_dPCA_axis(fname, opt)

def.conv = [];
def.plot = set_plot_opt('vik', 10);
def.format = [];
opt = safeStructAssign(def, opt);

S = load(fname);
data = S.data;

dpc = data.dpc;

ndim = size(dpc,1);

if ~isempty(opt.conv)
    for d=1:ndim
        for c=1:size(dpc,3)
            dpc(d,:,c) = nanconv(dpc(d,:,c), opt.conv, 'same');
        end
    end
end
if ~isempty(data.cutoff)
    dpc = dpc(:, data.cutoff(1):data.cutoff(2),:);
    tstamp = data.tstamp(data.cutoff(1):data.cutoff(2));
end

fh = figure('color', 'w', 'pos', [100 100 400, 200]);

ax = nan(ndim,1);
ymx = 0;
for d=1:ndim
    subplot(1,ndim,d);
    ax(d) = plot_multiepoch_trace({tstamp}, {squeeze(dpc(d,:,:))'}, [], opt);
    ymx = max(ymx, max(abs(ylim)));
    title(sprintf('%s: dim %d', data.target{d}{1}, data.target{d}{2}));
end
for d=1:ndim
    set(ax(d), 'ylim', [-1 1] * ymx);
end
opt.format.xlabel = 'Time (ms)';
opt.format.ylabel = 'Score';
format_panel(ax, opt.format);

end
