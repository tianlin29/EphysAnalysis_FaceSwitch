function opt = set_plot_opt_2cond(type1,type2, N)

switch type1
    case 'binary'
        col = [linspecer(2);linspecer(2)];
        N = 2;
    otherwise
        col = load_colormap(type1, N);
        col = [col;load_colormap(type2, N)];
end

opt.color = col;
opt.facecolor = col;
opt.edgecolor = col;
opt.linestyle = [repmat({'-'}, N, 1);repmat({'--'}, N, 1)];
opt.marker = [repmat({'.'}, N, 1);repmat({'.'}, N, 1)];
opt.markersize = ones(2*N,1) * 6;
opt.linewidth = [ones(N,1) * 1; ones(N,1) * 1];

end
