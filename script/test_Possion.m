run('../Initialize.m');
monkey = 'Woody'; % Nick, Woody
experiment = 'noCue'; % learnTask2, learnTask2_increaseSwitchFreq, learnTask3, learnTask4
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'PSTH'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'PSTH'); mkdir(InterimDir);

%%
fnd = load(get_file_path(monkey, experiment, 1, 'FND_sorted')).fnd;
r = fnd.raster(2); r = r{1}; % (unit, time, trial)
tstamp = fnd.tstamp{2};
I = tstamp>=0 & tstamp<=500; r = squeeze(r(1,I,:));

[ntime, ntrial] = size(r);
lambda_hat = sum(r(:)) / (ntime*ntrial/1000); % FR (sp/s)

ISI_i = cell(ntrial,1);
for tr = 1:ntrial
    spike_time = find(r(:,tr)==1)/1000;
    ISI_i{tr} = diff(spike_time);
end
ISI_i = cell2mat(ISI_i);
figure; hold on
histogram(ISI_i, linspace(0, max(ISI_i), 20), 'Normalization', 'count')

t = linspace(0, 0.5, 100);
y = lambda_hat * exp(-lambda_hat * t);
plot(t,y*100)
