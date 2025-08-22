run('../Initialize.m');

%% load fnd
fnd = load(get_file_path(FormattedDir, 'Woody', 'learnTask3', 1, 'FND')).fnd;
fnd = load('\\10.10.49.250\formatted_data\FaceSwitch_Monkey\Woody\20230923\Woody20230923_01_FND_sorted.mat').fnd;

%% check all the basic info of a session
properties(fnd); % all the properties

% useful for analysis
unitID = fnd.unitID; % [session id, channel id, unit id]
unit_quality = fnd.unit_quality; % 1~3
tstamp = fnd.tstamp; % precision is 1 ms
param = fnd.param; % all the params
ntrial = fnd.ntrial(1);
data = fnd.data; % {1, epochs}, (unit, time, trial), dtype: uint8
[nunit, ntime, ntrial] = size(fnd.data{1});

% not useful
unit_type = fnd.unit_type; % single, multi
monkey = fnd.monkey;
area = fnd.area;
valid_unit = fnd.valid_unit; % 0 or 1
unit_criteria_description = fnd.unit_criteria_description; % no hash
session_name = fnd.session_name;
description = fnd.description; % project name
misc = fnd.misc; % invalid trials, SNR
info = fnd.info; % []
session_info = fnd.session_info;
date = fnd.date;
version = fnd.version;
compression_type = fnd.compression_type;
PSTH_trend = fnd.PSTH_trend; % []

% trial id
trial = fnd.getp('trial');

%% check methods
methods(fnd); % all the methods

%% raster
% spike
r = fnd.raster(2); r = r{1}; % (unit, time, trial)

% specify condition
condID = fnd.getp('task_set');
r = fnd.raster_cond(condID); % when 'format' is []: {1, epoch}, {1, condition}, (unit, trial, time)
r = fnd.raster_cond(condID, 2, 'cell'); % when 'format' is 'cell': {1, condition}, (unit, trial, time)
r = fnd.raster_cond(condID, 2, 'array'); % when 'format' is 'array': (unit, trial, time, condition)

%% FR (mean across time)
% mean FR within some time range for each unit
r = fnd.FR({2, [0 500]}); % (unit, trial)

% mean FR within some time range for each unit, z-scored
r = fnd.FR({2, [0 500]}, [], true); % (unit, trial)

% mean FR within some time range for each condition
condID = fnd.getp('task_set');
r = fnd.FR({2, [0 500]}, condID); % (unit, condition)

%% PSTH (mean acorss trial)
% mean FR across all trials
psth = fnd.PSTH([], {'boxcar', 100}); % (unit, time)
figure; plot(psth{2}(1,:)); % psth of unit 1

% mean FR in each condition, cut off time of <50% trials contributes, not
% that useful
condID = fnd.getp('task_set');
psth = fnd.PSTH(condID, {'boxcar', 100}, [], true); % (unit, time, condition)
figure; plot(psth{2}(5,:,1)); % psth of unit 5

% mean FR across all trials, detrend (subtract mean)
psth = fnd.PSTH(condID, {'boxcar', 100}, [], true); % (unit, time)
figure; plot(psth{2}(5,:,1)); % psth of unit 5

% get psth of specified epoch
psth = fnd.PSTH_epoch(condID, {'boxcar', 100}, [], [], 2);

%% extract trial
trialFlg = ~isnan(fnd.getp('targ_cho'));
nfnd = fnd.extract_trial(trialFlg);

%% extract unit
% set other units as invalid
fnd.set_unit_criteria('monkey', 'Woody')

% create a copy of fnd
fnd_new = fnd.set_unit_criteria('customID', [1 2 1]);

% create a copy of fnd based on mean FR
r = fnd.FR({2, [100 500]});
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
fnd_new = fnd.set_unit_criteria('custom', I);

%% 一个session中的2个实验，unit一样，merge二者的trial
session = 1;
fnd_1 = load(get_file_path(FormattedDir, 'Woody', 'passiveLong', session, 'FND', 1)).fnd; % main task
fnd_2 = load(get_file_path(FormattedDir, 'Woody', 'passiveLong', session, 'FND', 2)).fnd; % passive task
fnd = FND.mergetrial({fnd_1, fnd_2});
experiment = fnd.getp('session');

%% basic data info
fnd.get_data_info

%% 神经数据的本质
% spike
r = fnd.raster(2); r = r{1}; % (unit, time, trial)
r_tmp = r(:,:,20); % (unit, time) raster的数据为0/1，记录神经元在1 ms的时间窗口内有无发放
figure; imagesc(r_tmp)

% 发放率：在一段时间窗口内统计spike/s
% 可以对FR执行zscore，即每个unit变为mean=0, std=1
FR = fnd.FR({2, [100 300]}); % (unit, trial)
figure; imagesc(FR)

% PSTH：在一些trial内统计spike/s
% PSTH的detrend是每个unit减去每个时刻点在所有trial上的平均值
psth = fnd.PSTH_epoch([], [], [], [], 2); psth = psth{2}; % (unit, time)
figure; imagesc(psth);







