run('../Initialize.m');

%% info
% Woody:
% rand40

% Nick:
% rand40

%% load data
load_opt.remove_fixbreak = 0; % whether to remove fixbreak/nochoice/nofix trials
load_opt.remove_eyedata = 1; % whether to remove eyelink data  (set it true if you don't analyze eye data)
% set remove_invalid_tr False because calculating cond_switch need
% information about all the trials
load_opt.remove_invalid_tr = 0; % whether to remo ve invalid trials (default: true, set true unless you have a reason)

load_opt.data_path = '\\10.10.49.250\'; % just simply load data from dataserver
load_opt.xls_file = 'Psych_FaceSwitch.xlsx'; % Psych_FaceSwitch or Train_FaceSwitch
load_opt.subject_name = {'Woody'}; % Woody or Nick
load_opt.experiment_name = 'rand40'; % if empty, load all the sessions of this subject
load_opt.load_D = 1;

if length(load_opt.subject_name)~=1
    error('Cannot process multiple subjects')
end

[D, session_id, block_id] = load_behavioral_data(load_opt); % by default, I will not remove_fixbreak, remove_eyedata, remove_invalid_tr

%% preprocess
I = cellfun(@(x) ~isempty(x.result), D); D = D(I); 
session_id = session_id(I); 
block_id = block_id(I);

% task setting
subj = repmat(load_opt.subject_name, length(D), 1);
date = get_date(D);
stim_id = cellfun(@(x) x.stim_id, D);

% performance
trial_id = cellfun(@(x) x.trial_id, D);
rt = cellfun(@(x) x.rt, D);
targ_cor = flip_targ_cor(D);
resp = flip_response(D); % for positive coh, targ_cor is set to be 2
result = cellfun(@(x) x.result{1}, D, 'uni', 0);

% valid
if isfield(D{1}, 'valid')
    valid = cellfun(@(x) x.valid, D);
else
    valid = true(size(D));
end

% save data
trial_data = save_trial_data(subj, date, D, valid, ...
    trial_id, session_id, block_id, ...
    stim_id, ...
    rt, ...
    targ_cor, resp, result);
save(fullfile(PreprocDir, [load_opt.subject_name{1}, '_', load_opt.experiment_name, '.mat']), 'trial_data');

%% additional functions
function trial_data = save_trial_data(subj, date, D, valid, ...
    trial_id, session_id, block_id, ...
    stim_id, ...
    rt, ...
    targ_cor, resp, result)

check = table(trial_id, stim_id, targ_cor, resp);

trial_data = cell(size(D));
for k = 1:length(D)
    % id
    trial_data{k}.subj = subj{k};
    trial_data{k}.date = date{k};
    trial_data{k}.session_id = session_id(k);
    trial_data{k}.block_id = block_id(k);
    % stimulus
    trial_data{k}.stim_id = stim_id(k);
    % time
    trial_data{k}.rt = rt(k);
    % response
    trial_data{k}.targ_cor = targ_cor(k);
    trial_data{k}.resp = resp(k);
    trial_data{k}.result = result(k);
end

% remove useless trials
I_rt = isnan(rt);
I_resp = ~strcmp(result,'CORRECT') & ~strcmp(result,'WRONG');
I_valid = valid==0;

fprintf('When saving data, remove the following invalid trials:\n');
fprintf('%d trials rt=NaN, %d trials inresponded, %d trials denoted as invalid\n', ...
    sum(I_rt), sum(I_resp), sum(I_valid));
I = I_rt | I_resp | I_valid; % I registers which trials need to be removed
trial_data = trial_data(~I);

end


function date = get_date(D)

match = cellfun(@(x) regexp(x.ownership, '([a-zA-Z]+)(\d{8})', 'tokens'), D, 'uni', 0); % letters + 8 numbers
date = cellfun(@(x) [x{1}{1} x{1}{2}], match, 'uni', 0);

end


function targ_cor = flip_targ_cor(D)

% targ_cor is set to be 2 for positive coh, flip resp accordingly
flip_target = cellfun(@(x)strcmp(x.flip_target, 'true'), D);
targ_cor = cellfun(@(x)x.targ_cor, D);
targ_cor(flip_target) = 3 - targ_cor(flip_target);

end


function resp = flip_response(D)

% targ_cor is set to be 2 for positive coh, flip resp accordingly
flip_target = cellfun(@(x)strcmp(x.flip_target, 'true'), D);
resp = cellfun(@(x)x.response, D);
resp(flip_target) = 3 - resp(flip_target);

end

