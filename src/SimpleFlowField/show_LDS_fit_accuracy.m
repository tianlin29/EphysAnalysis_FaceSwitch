function fh = show_LDS_fit_accuracy(lds_data)

fh = figure('pos', [100 100 300 300]);

Z = cellfun(@(x)x.Z, lds_data);
R2 = cellfun(@(x)x.R2, lds_data);

plot(Z, R2, 'ko-', 'linew', 2', 'markersize', 20);
xlabel('Z');
ylabel('R^2');
title('LDS fit accuracy');
format_panel(gca);


end