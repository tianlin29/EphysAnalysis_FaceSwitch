function [trial_data, session_id, block_id] = load_behavioral_data(load_opt)

% function [trial_data, session_id, block_id] = load_behavioral_data(load_opt)
%
% example for using load_behavioral_data:
%     load_opt.data_path = 'D:/LocalData/'; % Please change it to your own data path
%     load_opt.xls_file = 'behavioralData_example.xlsx'; % Please check labcommon/data_analysis/behavioralData_example.xlsx
%     load_opt.subject_name = {'008'}; % For multiple subjects, use {'002','008'}
%     load_opt.experiment_name = '2DSwitch_SameTarget';
%     load_opt.remove_fixbreak = true; % whether to remove fixbreak/nochoice/nofix trials
%     load_opt.remove_eyedata = true; % whether to remove eyelink data  (set it true if you don't analyze eye data)
%     load_opt.remove_invalid_tr = true; % whether to remove invalid trials (default: true, set true unless you have a reason)
%     load_opt.load_formatted_data = true; % If false, load raw data
%     [trial_data, session_id, block_id] = load_behavioral_data(load_opt);
%     plot(block_id); % check if you successfully remove all the invalid blocks

def.data_path = '';
def.xls_file = '';
def.subject_name = '';
def.experiment_name = '';
def.remove_fixbreak = NaN; % see combineTrialData
def.remove_eyedata = NaN; % see combineTrialData
def.remove_invalid_tr = NaN; % see combineTrialData
def.load_formatted_data = true;
load_opt = safeStructAssign(def, load_opt);

if ~isa(load_opt.subject_name, 'cell')
    error('load_opt.subject_name should be CELL.')
end

if isnan(load_opt.remove_fixbreak)
    error('remove_fixbreak option is not specified');
end
if isnan(load_opt.remove_eyedata)
    error('remove_eyedata option is not specified');
end
if isnan(load_opt.remove_invalid_tr)
    error('remove_invalid_tr option is not specified');
end

%% load from excel sheet
[out, ~] = load_excel_sheet(load_opt.xls_file, 'struct', false);

formatted_path = cellfun(@(x) fullfile(load_opt.data_path, x.formatted_folder, x.project, x.subject, x.date, x.formatted_file), out, 'uni', 0);
raw_path = cellfun(@(x) fullfile(load_opt.data_path, x.raw_folder, x.project, x.subject, x.date, x.raw_file), out, 'uni', 0);

subject_list = cellfun(@(x) x.subject(1:end-1), out, 'uni', 0);
experiment_list = cellfun(@(x) x.experiment, out, 'uni', 0);
valid = cellfun(@(x) str2double(x.valid), out);
invalid_block = cellfun(@(x) x.invalid_block, out, 'uni', 0);

%% pick out certain subjects and experiments
if isempty(load_opt.subject_name)
    load_opt.subject_name = subject_list;
end
if isempty(load_opt.experiment_name)
    load_opt.experiment_name = experiment_list;
end

% load multiple subjects
for k = 1:length(load_opt.subject_name)
    if k==1
        I = strcmp(subject_list, load_opt.subject_name{k}) & strcmp(experiment_list, load_opt.experiment_name) & valid;
    else
        I = I | (strcmp(subject_list, load_opt.subject_name{k}) & strcmp(experiment_list, load_opt.experiment_name) & valid);
    end
end
if all(~I)
    warning('No data satisfied the specified criteria');
    trial_data = [];
    session_id = [];
    block_id = [];
    return;
end

% print
if load_opt.load_formatted_data
    path_chosen = formatted_path(I);
else
    path_chosen = raw_path(I);
end

fprintf('Loading data from:\n')
for k = 1:length(path_chosen)
    fprintf('%s\n', path_chosen{k});
end

[trial_data, session_id, block_id] = combineTrialData(path_chosen, ...
    load_opt.remove_fixbreak, load_opt.remove_eyedata, load_opt.remove_invalid_tr);

%% remove invalid blocks
invalid_block = invalid_block(I);

rmv = zeros(size(trial_data));

for n_ses = 1:length(invalid_block)
    if ~isempty(invalid_block{n_ses})
        for k = 1:2:length(invalid_block{n_ses})
            I = session_id==n_ses & block_id==str2double(invalid_block{n_ses}(k));
            rmv(I) = 1;
            if sum(I)==0
                error('Certain block was not found.')
            end
            fprintf('Remove block %d from session %d (%d trials).\n', str2double(invalid_block{n_ses}(k)), n_ses, sum(I));
        end
    end
end

trial_data = trial_data(~rmv);
session_id = session_id(~rmv);
block_id = block_id(~rmv);

end
