run('../Initialize.m');

%% load data
load_opt.remove_fixbreak = 0; % whether to remove fixbreak/nochoice/nofix trials
load_opt.remove_eyedata = 1; % whether to remove eyelink data  (set it true if you don't analyze eye data)
% set remove_invalid_tr False because calculating cond_switch need
% information about all the trials
load_opt.remove_invalid_tr = 0; % whether to remo ve invalid trials (default: true, set true unless you have a reason)

load_opt.data_path = '\\10.10.49.250\'; % just simply load data from dataserver
load_opt.xls_file = 'Psych_FaceSwitch.xlsx'; % Psych_FaceSwitch or Train_FaceSwitch
load_opt.subject_name = {'Nick'}; % Woody or Nick
load_opt.experiment_name = 'faceColor'; % if empty, load all the sessions of this subject
load_opt.load_formatted_data = 1;

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
[block_type, block_type_processed] = get_block_type(D); % block type: 1 ..task 1, 2 ..task 2, 3 ..switch; block_type_processed: 1 ..single, 2 ..switch
cond = cellfun(@(x) x.curr_task_set, D);
[cond_switch, n_after_switch, n_after_nonswitch] = get_cond_switch(D); % 1 ..non-switch, 2 ..switch
morph_level = cellfun(@(x) x.morph_level(1), D);
color_level = cellfun(@(x) x.color_level(1), D);
fluc = cellfun(@(x) x.stim_morphvals'*2-1, D, 'uni', 0); % make sure column is time
stim_dur = cellfun(@(x) x.stim_dur, D);
stim_size = cellfun(@(x) x.stim_aperture(3), D);
stim_set = get_stim_set(D);

% performance
trial_id = cellfun(@(x) x.trial_id, D); % start from 0
rt = cellfun(@(x) x.rt, D);
targ_cor = flip_targ_cor(D, cond);
resp = flip_response(D, cond); % for positive coh, targ_cor is set to be 2
result = cellfun(@(x) x.result{1}, D, 'uni', 0);

% valid
if isfield(D{1}, 'valid')
    valid = cellfun(@(x)x.valid, D);
else
    valid = true(size(D));
end

if strcmp(load_opt.experiment_name, 'noCue')
    valid_trial = trial_id>400;
    valid = valid & valid_trial;
    fprintf('This is noCue experiment. The first 400 trials in each session are defined as invalid trials.\n')
end

% save data
[trial_data, check] = save_trial_data(subj, date, D, valid, ...
    trial_id, session_id, block_id, ...
    morph_level, color_level, stim_size, stim_set, ...
    rt, stim_dur, ...
    targ_cor, resp, result, ...
    block_type, block_type_processed, cond, cond_switch, n_after_switch, n_after_nonswitch, ...
    fluc);
if isempty(load_opt.experiment_name)
    save(fullfile(PreprocDir, [load_opt.subject_name{1}, '_training.mat']), 'trial_data');
else
    save(fullfile(PreprocDir, [load_opt.subject_name{1}, '_', load_opt.experiment_name, '.mat']), 'trial_data');
end

%% additional functions
function [trial_data, check] = save_trial_data(subj, date, D, valid, ...
    trial_id, session_id, block_id, ...
    morph_level, color_level, stim_size, stim_set, ...
    rt, stim_dur, ...
    targ_cor, resp, result, ...
    block_type, block_type_processed, cond, cond_switch, n_after_switch, n_after_nonswitch, ...
    fluc)

check = table(trial_id, cond, result, cond_switch, n_after_switch, targ_cor);

trial_data = cell(size(D));
for k = 1:length(D)
    % id
    trial_data{k}.subj = subj{k};
    trial_data{k}.date = date{k};
    trial_data{k}.session_id = session_id(k);
    trial_data{k}.block_id = block_id(k);
    % coherence
    trial_data{k}.morph_level = morph_level(k);
    trial_data{k}.color_level = color_level(k);
    trial_data{k}.stim_size = stim_size(k);
    trial_data{k}.stim_set = stim_set{k};
    % time
    trial_data{k}.rt = rt(k);
    trial_data{k}.stim_dur = stim_dur(k);
    % response
    trial_data{k}.targ_cor = targ_cor(k);
    trial_data{k}.resp = resp(k);
    trial_data{k}.result = result(k);
    % context
    trial_data{k}.block_type = block_type(k);
    trial_data{k}.block_type_processed = block_type_processed(k);
    trial_data{k}.cond = cond(k);
    trial_data{k}.cond_switch = cond_switch(k);
    trial_data{k}.n_after_switch = n_after_switch(k);
    trial_data{k}.n_after_nonswitch = n_after_nonswitch(k);
    % fluc
    trial_data{k}.fluc = fluc{k};
end

% remove useless trials for DDM fit
I_rt = isnan(rt);
I_cond_switch = isnan(cond_switch);
I_resp = ~strcmp(result,'CORRECT') & ~strcmp(result,'WRONG');
I_valid = valid==0;

fprintf('When saving data, remove the following invalid trials:\n');
fprintf('%d trials rt=NaN, %d trials condswitch=NaN, %d trials inresponded, %d trials denoted as invalid\n', ...
    sum(I_rt), sum(I_cond_switch), sum(I_resp), sum(I_valid));
I = I_rt | I_cond_switch | I_resp | I_valid; % I registers which trials need to be removed
trial_data = trial_data(~I);

end


function date = get_date(D)

match = cellfun(@(x) regexp(x.ownership, '([a-zA-Z]+)(\d{8})', 'tokens'), D, 'uni', 0); % letters + 8 numbers
date = cellfun(@(x) [x{1}{1} x{1}{2}], match, 'uni', 0);

end


function targ_cor = flip_targ_cor(D, cond)

% targ_cor is set to be 2 for positive coh, flip resp accordingly
flip_target1 = cellfun(@(x)strcmp(x.flip_target1, 'true'), D);
flip_target2 = cellfun(@(x)strcmp(x.flip_target2, 'true'), D);
targ_cor = cellfun(@(x)x.targ_cor, D);
targ_cor(cond==1 & flip_target1) = 3 - targ_cor(cond==1 & flip_target1);
targ_cor(cond==2 & flip_target2) = 3 - targ_cor(cond==2 & flip_target2);

end


function resp = flip_response(D, cond)

% targ_cor is set to be 2 for positive coh, flip resp accordingly
flip_target1 = cellfun(@(x) strcmp(x.flip_target1, 'true'), D);
flip_target2 = cellfun(@(x) strcmp(x.flip_target2, 'true'), D);
resp = cellfun(@(x) x.response, D);
resp(cond==1 & flip_target1) = 3 - resp(cond==1 & flip_target1);
resp(cond==2 & flip_target2) = 3 - resp(cond==2 & flip_target2);

end


function [cond_switch, n_after_switch, n_after_nonswitch] = get_cond_switch(D)

cond = cellfun(@(x) x.curr_task_set, D);
result = cellfun(@(x) x.result{1}, D, 'uni', 0);
date = cellfun(@(x) x.ownership(find(x.ownership=='\', 1, 'last')+1 : find(x.ownership=='_', 1, 'last')-1), D, 'uni', 0);

% exclude NOFIX, then calculate cond_switch
I = strcmp(result, 'NOFIX');
cond_switch_pre = [NaN; abs(diff(cond(~I)))+1];
cond_switch = nan(size(D));
cond_switch(~I) = cond_switch_pre;

% the first trial in each block is neither a switch nor a non-switch trial
[~, ~, ic] = unique(date);
I_first_trial = logical([1; diff(ic)]); % the first trial in each session
cond_switch(I_first_trial) = NaN;

% number of trials after switch
n_after_switch = nan(size(D));
count = 0;
for k = 1:length(cond_switch)
    if isnan(cond_switch(k))
        n_after_switch(k) = NaN;
    elseif cond_switch(k)==2
        n_after_switch(k) = 0; % 0 ..a switch trial
        count = 0; % reset trial count
    elseif strcmp(result{k},'CORRECT') || strcmp(result{k},'WRONG')
        count = count + 1;
        n_after_switch(k) = count;
    end
    
    if k~=length(cond_switch) && I_first_trial(k+1)==1 % the last trial in a session
        count = 0; % reset trial count
    end
end

% after how many nonswitch trials
n_after_nonswitch = nan(size(n_after_switch));
idx = find(n_after_switch==0);
n_after_nonswitch(n_after_switch==0) = n_after_switch(idx-1);

% check
check = table(cond, result, cond_switch, n_after_switch, n_after_nonswitch);
% if any(n_after_nonswitch==0)
%     error('consecutive switch trials') % if we switch between slow and fast block, this kind of trial could exist
% end

end


function [block_type, block_type_processed] = get_block_type(D)

% block type: 1 ..task 1, 2 ..task 2, 3 ..switch; block_type_processed: 1 ..single, 2 ..switch
block_type = nan(size(D));
block_type_processed = nan(size(D));

if ~isfield(D{1}, 'block_type')
    warning('No block type in the data.')
    return
end

for n = 1:length(D)
    tmp = D{n}.block_type;
    if strcmp(tmp, 'task 1')
        block_type(n) = 1;
        block_type_processed(n) = 1;
    elseif strcmp(tmp, 'task 2')
        block_type(n) = 2;
        block_type_processed(n) = 1;
    elseif strcmp(tmp, 'switch')
        block_type(n) = 3;
        block_type_processed(n) = 2;
    else
        error('check block_type')
    end
end

end


function stim_set = get_stim_set(D)

stim_set = cell(size(D));
for n = 1:length(D)
    if D{n}.curr_task_set==1
        stim_set{n} = D{n}.stim_set1;
    else
        stim_set{n} = D{n}.stim_set2;
    end
end

end
