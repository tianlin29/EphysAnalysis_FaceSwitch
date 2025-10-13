run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'track_neuron_Fraser2012')));
monkey = 'Nick'; % Nick, Woody
experiment = 'learnTask2'; % learnTask2, learnTask3, learnTask4, faceColor
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'trackNeuron'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'trackNeuron'); mkdir(InterimDir);

%% transform nev data to proper format
% loading each session takes ~1 min

% set option
n_electrode = 96;
invalid_waveform_unit = 255;
nev_opt = struct('skip_sp_waveform', true, 'skip_ns2', true, 'skip_ns3', true, 'skip_ns5', true, ...
    'strobed_events', false, 'delete_nev_data', false, 'delete_ns2_data', true, 'delete_ns3_data', true);
mat_files = arrayfun(@(x) get_file_path(monkey, experiment, x, 'sorted.mat', [], '/mnt/raid-volume/sorted_data/'), (1:n_files)', 'uni', 0);
nev_files = arrayfun(@(x) get_file_path(monkey, experiment, x, 'sorted.nev', [], '/mnt/raid-volume/sorted_data/'), (1:n_files)', 'uni', 0);

for s = 1:n_files
    fprintf('session %d\n', s)

    % load mat file and get waveform
    sort_info = load(mat_files{s}).sort_info;

    unit_ID = sort_info.unit_ID;
    mean_firing_rate = sort_info.mean_firing_rate;
    waveform = sort_info.waveform;

    mean_firing_rate = mean_firing_rate(unit_ID(:,2)~=0, :);
    waveform = waveform(unit_ID(:,2)~=0, :);
    waveform_array = num2cell(waveform, 2);

    % load nev file
    nf = nevFile(nev_files{s}, nev_opt);

    % get spike time, channel id and unit id
    spike = nf.spike(1:n_electrode, 2:end); % start from 2 to remove hash
    channel_ID = repmat((1:n_electrode)', 1, size(spike,2));
    unit_ID = repmat(1:size(spike,2), n_electrode, 1);
    spike_array = reshape(spike', [], 1);
    spike_array = cellfun(@(x) x/1000, spike_array, 'uni', 0); % ms -> s
    channel_ID_array = reshape(channel_ID', [], 1);
    unit_ID_array = reshape(unit_ID', [], 1);
    I = cellfun(@(x) isempty(x), spike_array);
    spike_array(I) = []; channel_ID_array(I) = []; unit_ID_array(I) = [];

    % summary
    valid = ~isnan(mean_firing_rate); % remove units with nan FR
    channel = channel_ID_array(valid);
    unit = unit_ID_array(valid);
    spiketimes = spike_array(valid);
    wmean = waveform_array(valid);
    save(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)), 'channel', 'unit', 'spiketimes', 'wmean')
end

%% [test] find the same unit across sessions
s = 1;

[channel, unit, spiketimes, wmean] = deal(cell(2, 1));
D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)));
channel{1} = D.channel;
unit{1} = D.unit;
spiketimes{1} = D.spiketimes;
wmean{1} = D.wmean;

D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s+1)));
channel{2} = D.channel;
unit{2} = D.unit;
spiketimes{2} = D.spiketimes;
wmean{2} = D.wmean;

% run the function of unitIdentification
auto = computeDailyAutocorrelations(spiketimes);
autoscore = computeAutoScore(auto);

base = computeDailyBaserates(spiketimes);
basescore = computeBaseScore(base);
if opt.overwrite
    correlations = computeDailyCorrelations(spiketimes);
    save(fullfile(InterimDir, sprintf('correlations_%s_%s_session%d.mat', monkey, experiment, s)), 'correlations')
else
    correlations = load(fullfile(InterimDir, sprintf('correlations_%s_%s_session%d.mat', monkey, experiment, s))).correlations;
end
wavescore = computeWaveScore(wmean);

[survival, corrscore, wavescore, autoscore, basescore, score] = iterateSurvival(channel, unit, correlations, wavescore, autoscore, basescore);

% fh = gaussianMix(channel, unit, corrscore, wavescore, autoscore, basescore, survival);
% print(fh, '-dpdf', fullfile(FigDir, sprintf('result_%s_%s_session%d.pdf', monkey, experiment, s)))

%% [test] show results of 2 sessions
fprintf('%d units in day 1 and 2 are detected as the same units.\n', sum(survival{1}(:)))

% plot which units are the same
figure;
imagesc(survival{1})
format_panel(gcf, 'ylabel', 'Unit id of day 1', 'xlabel', 'Unit id of day 2')

% plot waveform of example units in two sessions
[row, col] = find(survival{1}==1);

figure;
for n = 1:50
    subplot(5,10,n); hold on
    plot(wmean{1,1}{row(n)})
    plot(wmean{2,1}{col(n)})
    % ylim([-200 200])
    title(sprintf('unit %d', n))
end

%% find the same unit across sessions
clear opt;
opt.overwrite = 1;

survival = cell(n_files-1, 1);
for s = 1:n_files-1
    [channel, unit, spiketimes, wmean] = deal(cell(2, 1));
    D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)));
    channel{1} = D.channel;
    unit{1} = D.unit;
    spiketimes{1} = D.spiketimes;
    wmean{1} = D.wmean;

    D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s+1)));
    channel{2} = D.channel;
    unit{2} = D.unit;
    spiketimes{2} = D.spiketimes;
    wmean{2} = D.wmean;

    % run the function of unitIdentification
    auto = computeDailyAutocorrelations(spiketimes);
    autoscore = computeAutoScore(auto);

    base = computeDailyBaserates(spiketimes);
    basescore = computeBaseScore(base);
    if opt.overwrite
        correlations = computeDailyCorrelations(spiketimes);
        save(fullfile(InterimDir, sprintf('correlations_%s_%s_session%d.mat', monkey, experiment, s)), 'correlations')
    else
        correlations = load(fullfile(InterimDir, sprintf('correlations_%s_%s_session%d.mat', monkey, experiment, s))).correlations;
    end
    wavescore = computeWaveScore(wmean);

    survival(s) = iterateSurvival(channel, unit, correlations, wavescore, autoscore, basescore);
end

save(fullfile(InterimDir, sprintf('survival_%s_%s.mat', monkey, experiment)), 'survival')

%%
% D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)));

% load survivial file between all the sessions
load(fullfile(InterimDir, sprintf('survival_%s_%s.mat', monkey, experiment)))

% determine units survived from the beginning to ending
for s = 1:n_files-1
    if s==1
        sur_1_end = survival{s};
    else
        sur_1_end = sur_1_end * survival{s};
    end
end
[id_1, id_end] = find(sur_1_end==1);

[id, wave] = deal(cell(1, n_files));
D_1 = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, 1)));
id{1} = [D_1.channel(id_1), D_1.unit(id_1)];
wave{1} = D_1.wmean(id_1);
D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, n_files)));
id{end} = [D.channel(id_end), D.unit(id_end)];
wave{end} = D.wmean(id_end);

% determine units survived in the middle sessions
for s = 2:n_files-1
    for ss = 1:s-1
        if ss==1
            sur_1_end = survival{ss};
        else
            sur_1_end = sur_1_end * survival{ss};
        end
    end
    [id_1, id_end] = find(sur_1_end==1);
    
    id_1 = [D_1.channel(id_1), D_1.unit(id_1)];
    D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)));
    id_end = [D.channel(id_end), D.unit(id_end)];

    I = ismember(id_1, id{1}, 'rows');
    id{s} = id_end(I, :);
    wave{s} = D.wmean(I);
end

figure;
for n = 1:25
    subplot(5,5,n); hold on
    for s = 1:n_files
        plot(wave{s}{n})
    end
end


%%
B = [0 1;
    0 0];
A = [0 1;
    0 0;
    1 1;
    0 1];

[isPresent, loc] = ismember(A, B, 'rows');



%%
count = nan(n_files-1, 1);
id = cell(1, n_files);
for s = 1:n_files-1
    if s==1
        track = survival{s};
        [id{1}, id{2}] = find(track==1);
    else
        track = track * survival{s};
        [~, id{s+1}] = find(track==1);
    end
    count(s) = sum(track(:));
    fprintf('%d units\n', count(s))
end

track_wmean = cell(1, n_files);
for s = 1:n_files
    D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)));
    track_wmean{s} = D.wmean(id{s});
end


per = count/212;
per = [1; per];
figure; stairs(per)
xlim([1 150])
ylim([0 1])
xlabel('Survival (days)')
ylabel('Remaining')

%%

load('../interim/data.mat')
[survival, score, corrscore, wavescore, autoscore, basescore, correlations] = unitIdentification(channel, unit, spiketimes, wmean); % takes ~30 sec


%% show results of 3 sessions
fprintf('%d units in day 1 and 2 are detected as the same units.\n', sum(survival{1}(:)))
fprintf('%d units in day 2 and 3 are detected as the same units.\n', sum(survival{2}(:)))

% print unit id
intersection = survival{1} * survival{2};
[day1_id_list, day3_id_list] = find(intersection);
day2_id_list = find(any(survival{1}(day1_id_list,:), 1));
fprintf('%d units persisted in 3 days.\n', sum(intersection(:)))
fprintf('day1    → day2    → day3   \n')
for n = 1:length(day1_id_list)
    fprintf('ch%2d u%d → ch%2d u%d → ch%2d u%d\n', ...
        channel{1}(day1_id_list(n)), unit{1}(day1_id_list(n)), ...
        channel{2}(day2_id_list(n)), unit{2}(day2_id_list(n)), ...
        channel{3}(day3_id_list(n)), unit{3}(day3_id_list(n)));
end

% plot waveform in 3 days
figure;
for n = 1:25
    subplot(5,5,n); hold on
    plot(wmean{1,1}{day1_id_list(n)})
    plot(wmean{2,1}{day2_id_list(n)})
    plot(wmean{3,1}{day3_id_list(n)})
    ylim([-300 300])
    title(sprintf('unit %d', n))
end

