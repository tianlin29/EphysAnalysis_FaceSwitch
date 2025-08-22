run('../Initialize_face_switch.m');
addpath C:\Engine\EphysPreprocess\FIRA2FND

monkey = 'Nick';
experiment = 'learnTask4';

%% copy file
[~, n_files] = get_file_path(monkey, experiment);
for n = 1:n_files
    % FND
    file_path_local = get_file_path(monkey, experiment, n, 'FND_sorted');
    file_path_remote = strrep(file_path_local, 'D:\LocalData\formatted_data', '\\10.10.49.250\formatted_data');
    mkdir(fileparts(file_path_local)) % create folder in local dir
    copyfile(file_path_remote, fileparts(file_path_local))

    % eye
    file_path_local = get_file_path(monkey, experiment, n, 'eye');
    file_path_remote = strrep(file_path_local, 'D:\LocalData\formatted_data', '\\10.10.49.250\formatted_data');
    copyfile(file_path_remote, fileparts(file_path_local))
end

%% alignment
% event name: fp_on, cue_on (actually not existed), targ_on, stim_on, stim_off, fp_off, response, response_edf, feedback
% go delay: from stim off to fp off is 500~1000 ms
[~, n_files] = get_file_path(monkey, experiment);
for n = 1:n_files
    fprintf('session %d\n', n);

    fira_fname = get_file_path(monkey, experiment, n, 'FIRA_sorted', [], '\\10.10.49.250\formatted_data');
    new_fnd_fname = get_file_path(monkey, experiment, n, 'FND_sorted_preprocessed');

    % align to stim off
    clear opt
    opt.include_hash = false;
    opt.align_str = {{'event', 'stim_on';      'start_offset', -100; 'end_offset', 100; 'limits', {{'targ_on',0}, {'stim_on',100}}}
                     {'event', 'stim_on';      'start_offset', -100; 'end_offset', 700; 'limits', {{'targ_on',0}, {'stim_off',0}}}
                     {'event', 'stim_off';     'start_offset', -100; 'end_offset', 700; 'limits', {{'stim_on',0}, {'fp_off',0}}}
                     {'event', 'fp_off';       'start_offset', -100; 'end_offset', 200; 'limits', {{'stim_off',0}, {'response_edf',0}}}
                     {'event', 'response_edf'; 'start_offset', -600; 'end_offset', 200; 'limits', {{'stim_on',0}, {'response_edf',200}}}};
    fnd = FIRA2FND_FaceSwitchFix(fira_fname, opt);
    % fnd = remove_low_FR_unit(fnd);
    save(new_fnd_fname, 'fnd')
end

% check FR to stim off
unit = 40; trial = 3;
r = fnd.raster(2); r = r{1}; % (unit, time, trial)
[nunit, ntime, ntrial] = size(r);
tmp = squeeze(r(1,:,:)); I = isnan(tmp);
stim_dur = fnd.getp('stim_dur'); stim_dur = stim_dur(1,:);
tstamp = fnd.tstamp{2};

[t_raster, t_stim] = deal(nan(ntrial,1));
for tr = 1:ntrial
    r1 = squeeze(r(1,:,tr));
    idx = sum(~isnan(r1));
    if idx<ntime
        t_raster(tr) = tstamp(idx);
    end

    t_stim(tr) = round(stim_dur(tr)*1000);
end

figure; hold on
plot(t_raster, 'Color', 'black')
plot(t_stim, 'Color', 'red')
dif = t_raster - t_stim;

%% functions
function new_fnd = remove_low_FR_unit(fnd)

r = fnd.FR({2, [100 500]}, [], false, false, true); % (unit, trial)
I = nanmean(r, 2)>=1; % mean FR should >= 1 Hz
new_fnd = fnd.set_unit_criteria('custom', I);
fprintf('%d units have FR < 1 sp/s. Removed.\n', sum(~I))

end
