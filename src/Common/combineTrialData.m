
function [trial_data, session_id, block_id] = combineTrialData(file_names, remove_fixbreak, remove_eyedata, remove_invalid_tr)
% 
% function [trial_data, session_id, block_id] = combineTrialData(file_names, remove_fixbreak, remove_eyedata, remove_invalid_tr)
% 
% Combines trial_data structures (ROME file format) across multiple data files. The trial_data structures in file_names do not 
% need to be identical.
% 
% remove_fixbreak       1                   remove 'FIXBREAK', 'NOCHOICE', and 'NOFIX'
%                       [],0,or undefined   do not remove the above trials
% 
% remove_eyedata        1                   remove edf data
%                       2                   remove all eye data
%                       [],0,or undefined   do not remove the data
% 
% remove_invalid_tr     1                   remove trials that are considered 'invalid' by preformat program, 
%                                           including trials with saccade timing problem, stimulus timing 
%                                           problem, and/or trials under training mode
%                       [],0,or undefined   do not remove the above trials 

% RK, 06/2012
% RK & CH, updated to allow removal of non-edf eye position data

    if nargin<2 || isempty(remove_fixbreak)
        remove_fixbreak = false;
    end
    if nargin<3 || isempty(remove_eyedata)
        remove_eyedata = false;
    end
    if nargin<4 || isempty(remove_invalid_tr)
        remove_invalid_tr = true;
    end
    
    if ~iscell(file_names)
        if contains(file_names, '*')
            f = dir(file_names);
            file_names = arrayfun(@(x)[x.folder '/' x.name], f, 'uni', 0);
        else
            file_names = {file_names};
        end
    end
    
    hw = waitbar_text(0);
    for i = 1 : length(file_names)
        D = load(file_names{i});
        D = addSource(D.trial_data, file_names{i});
        if isrow(D)
            D = D';
        end
        if remove_eyedata>0
            if iscell(D)
                if isfield(D{1},'eyelink')
                    for j = 1 : length(D)
                        if isfield(D{j}.eyelink, 'eye_data')
                            D{j}.eyelink = rmfield(D{j}.eyelink, 'eye_data');
                        end
                    end
                    %D = cellfun(@(x)rmfield(x,'eyelink'), D , 'UniformOutput', false);
                end
                if remove_eyedata==2 && isfield(D{1},'eye_data')
                    D = cellfun(@(x)rmfield(x,'eye_data'), D , 'UniformOutput', false);
                end
            elseif ~iscell(D)
                if isfield(D(1),'eyelink')
                    for j = 1 : length(D)
                        if isfield(D{j}.eyelink, 'eye_data')
                            D(j).eyelink = rmfield(D(j).eyelink, 'eye_data');
                        end
                    end
                    %D = arrayfun(@(x)rmfield(x,'eyelink'), D , 'UniformOutput', false);
                end
                if remove_eyedata==2 && isfield(D{1},'eye_data')
                    D = arrayfun(@(x)rmfield(x,'eye_data'), D , 'UniformOutput', false);
                end
                D = cat(1,D{:});
            end
        end
        if remove_fixbreak
            D(cellfun(@(x)isempty(x.result) || any(strcmp(x.result{1}, {'FIXBREAK', 'NOCHOICE', 'NOFIX'})), D)) = [];
        end
        if remove_invalid_tr && isfield(D{1}, 'valid')
            D(cellfun(@(x)~x.valid, D)) = [];
        end
        if i==1
            trial_data = D;
            session_id = ones(length(D),1);
            if isfield(D{1},'block_id')
                block_id = cellfun(@(x)x.block_id, D);
            elseif isfield(D{1},'trial_type')
                block_id = cal_block_id(D);
            else % neural data has neither block_id nor trial_type
                block_id = ones(size(D));
            end
        else
            trial_data = [trial_data; D]; %#ok<AGROW>
            session_id = [session_id; ones(length(D),1) * i];  %#ok<AGROW>
            if isfield(D{1},'block_id')
                block_id = [block_id; cellfun(@(x)x.block_id, D)];  %#ok<AGROW>
            elseif isfield(D{1},'trial_type')
                block_id = [block_id; cal_block_id(D)];  %#ok<AGROW>
            else % neural data has neither block_id nor trial_type
                block_id = [block_id; ones(size(D))];  %#ok<AGROW>
            end
        end
        waitbar_text(i/length(file_names), hw);
    end
    waitbar_text('close', hw);
end


function block_id = cal_block_id(D)
    trial_type_logical = cellfun(@(x)isfield(x, 'trial_type') && isequal(x.trial_type, 'last'), D);
    block_id = [0; cumsum(trial_type_logical(1:end-1))] + 1;
end

function trial_data = addSource(trial_data, source)
    for i = 1 : length(trial_data)
        if iscell(trial_data)
            if ~isempty(trial_data{i})
                trial_data{i}.ownership = source;
            end
        else
            if ~isempty(trial_data(i))
                trial_data(i).ownership = source;
            end
        end
    end
end

