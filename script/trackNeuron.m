run('../Initialize.m');
addpath(genpath(fullfile(MainDir, 'external', 'track_neuron_Fraser2012')));
monkey = 'Nick'; % Nick, Woody
experiment = 'faceColor'; % learnTask2, learnTask3, learnTask4, faceColor
[~, n_files] = get_file_path(monkey, experiment);
FigDir = fullfile(MainFigDir, 'trackNeuron'); mkdir(FigDir);
InterimDir = fullfile(MainInterimDir, 'trackNeuron'); mkdir(InterimDir);


%% transform nev data to proper format
% loading each session takes ~1 min

% set option
invalid_waveform_unit = 255;
nev_opt = struct('skip_sp_waveform', true, 'skip_ns2', true, 'skip_ns3', true, 'skip_ns5', true, ...
    'strobed_events', false, 'delete_nev_data', false, 'delete_ns2_data', true, 'delete_ns3_data', true);
files = {'/mnt/raid-volume/sorted_data/FaceSwitch_Monkey/Nick/20250920/Nick20250920_01';
    '/mnt/raid-volume/sorted_data/FaceSwitch_Monkey/Nick/20250922/Nick20250922_01';
    '/mnt/raid-volume/sorted_data/FaceSwitch_Monkey/Nick/20250923/Nick20250923_01'};

nsession = length(files);
nelectrode = 96;
[channel, unit, spiketimes, wmean] = deal(cell(nsession,1));

for s = 1:nsession
    % load mat file and get waveform
    mat_file = [files{s}, '_sorted.mat'];
    sort_info = load(mat_file).sort_info;

    unit_ID = sort_info.unit_ID;
    mean_firing_rate = sort_info.mean_firing_rate;
    waveform = sort_info.waveform;

    mean_firing_rate = mean_firing_rate(unit_ID(:,2)~=0, :);
    waveform = waveform(unit_ID(:,2)~=0, :);
    waveform_array = num2cell(waveform, 2);

    % load nev file
    nev_file = [files{s}, '_sorted.nev'];
    nf = nevFile(nev_file, nev_opt);

    % get spike time, channel id and unit id
    spike = nf.spike(1:nelectrode, 2:end); % start from 2 to remove hash
    channel_ID = repmat((1:nelectrode)', 1, size(spike,2));
    unit_ID = repmat(1:size(spike,2), nelectrode, 1);
    spike_array = reshape(spike', [], 1);
    channel_ID_array = reshape(channel_ID', [], 1);
    unit_ID_array = reshape(unit_ID', [], 1);
    I = cellfun(@(x) isempty(x), spike_array);
    spike_array(I) = []; channel_ID_array(I) = []; unit_ID_array(I) = [];

    % summary
    valid = ~isnan(mean_firing_rate); % remove units with nan FR
    channel{s,1} = channel_ID_array(valid);
    unit{s,1} = unit_ID_array(valid);
    spiketimes{s,1} = spike_array(valid);
    wmean{s,1} = waveform_array(valid);
end

% save
save('../interim/data.mat', 'channel', 'unit', 'spiketimes', 'wmean')

%% find the same unit across sessions
load('../interim/data.mat')
[survival, score, corrscore, wavescore, autoscore, basescore, correlations] = unitIdentification(channel, unit, spiketimes, wmean); % takes ~30 sec

%% show results of 2 sessions
fprintf('%d units in day 1 and 2 are detected as the same units.\n', sum(survival{1}(:)))

% plot which units are the same
figure;
imagesc(survival{1})
format_panel(gcf, 'ylabel', 'Unit id of day 1', 'xlabel', 'Unit id of day 2')

% plot waveform of example units in two sessions
[row, col] = find(survival{1}==1);

figure;
for n = 1:25
    subplot(5,5,n); hold on
    plot(wmean{1,1}{row(n)})
    plot(wmean{2,1}{col(n)})
    ylim([-200 200])
    title(sprintf('unit %d', n))
end

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

