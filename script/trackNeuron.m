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

%% identify unit id
% load survivial file between all the sessions
load(fullfile(InterimDir, sprintf('survival_%s_%s.mat', monkey, experiment)))

% determine units survived from the beginning to ending
n_files = 6; % 5 sessions only
for s = 1:n_files-1
    if s==1
        sur_1_end = survival{s};
    else
        sur_1_end = sur_1_end * survival{s};
    end
end
[id_1, id_end] = find(sur_1_end==1);

id = nan(length(id_1), n_files);
id(:,1) = id_1;
id(:,end) = id_end;

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
    
    I = ismember(id_1, id(:,1), 'rows');
    id(:,s) = id_end(I);
end
save(fullfile(InterimDir, sprintf('id_%s_%s.mat', monkey, experiment)), 'id')

%% check waveform
wmean = cell(size(id));
for s = 1:n_files
    D = load(fullfile(InterimDir, sprintf('data_%s_%s_session%d.mat', monkey, experiment, s)));
    wmean(:,s) = D.wmean(id(:,s));
end

figure;
for n = 1:52
    subplot(6,10,n); hold on
    for s = 1:n_files
        plot(wmean{n,s})
    end
    title(sprintf('unit %d', n))
end
format_panel(gcf)

%% save fnd of track neurons for each session
id = load(fullfile(InterimDir, sprintf('id_%s_%s.mat', monkey, experiment))).id;

for n = 1:n_files
    fprintf('session %d\n', n)

    fnd = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    [nunit, ntime, ntrial] = size(fnd.data{1});
    fnd = fnd.set_unit_criteria('custom', ismember((1:nunit)', id(:,n)));

    new_file_name = get_file_path(monkey, experiment, n, 'FND_sorted_track');
    save(new_file_name, 'fnd')
end

%% merge data across sessions
id = load(fullfile(InterimDir, sprintf('id_%s_%s.mat', monkey, experiment))).id;

for n = 1:n_files
    fprintf('session %d\n', n)

    fnd_tmp = load(get_file_path(monkey, experiment, n, 'FND_sorted')).fnd;
    [nunit, ntime, ntrial] = size(fnd_tmp.data{1});
    fnd_tmp = fnd_tmp.set_unit_criteria('custom', ismember((1:nunit)', id(:,n)));
    fnd_tmp.setp('session', n*ones(size(fnd_tmp.getp('stimID'))))

    if n==1
        fnd = fnd_tmp;
    else
        % add spike data
        for e = 1:3
            fnd.data{e} = cat(3, fnd.data{e}, fnd_tmp.data{e});
        end

        % add param
        fields = fieldnames(fnd.param);
        for i = 1:length(fields)
            fnd.param.(fields{i}) = cat(2, fnd.param.(fields{i}), fnd_tmp.param.(fields{i}));
        end

        % change ntrial
        if n==n_files
            ntrial = size(fnd.data{1}, 3);
            fnd.ntrial = ntrial*ones(size(fnd.ntrial));
        end
    end
end

save(fullfile(MainDir, 'data', 'preproc', sprintf('%s_%s_FND.mat', monkey, experiment)), 'fnd')







