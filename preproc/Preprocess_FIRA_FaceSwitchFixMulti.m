run('../Initialize.m');
addpath C:\Engine\EphysPreprocess\FIRA2FND

%% preprocess (rand40)
monkey = 'Nick';
experiment = 'rand40';
[~, n_files] = get_file_path(monkey, experiment);

for n = 1:n_files
    fprintf('%d\n', n);
    fira_fname = get_file_path(monkey, experiment, n, 'FIRA_sorted', [], '\\10.10.49.250\formatted_data');
    nfnd_fname = get_file_path(monkey, experiment, n, 'FND_sorted_preprocessed', [], 'D:\LocalData\formatted_data');

    % new alignment
    clear opt
    opt.include_hash = false;
    opt.align_str = {{'event', 'stim_on';      'start_offset', -100; 'end_offset', 700; 'limits', {{'targ_on',0}, {'response_edf',0}}}
                     {'event', 'stim_off';     'start_offset', -200; 'end_offset', 500; 'limits', {{'stim_on',0}, {'response_edf',0}}}
                     {'event', 'fp_off';       'start_offset', -200; 'end_offset', 500; 'limits', {{'stim_on',0}, {'response_edf',200}}}
                     {'event', 'response_edf'; 'start_offset', -1000;'end_offset', 200; 'limits', {{'stim_on',0}, {'response_edf',200}}}};
    fnd = FIRA2FND_FaceSwitchFixMulti(fira_fname, opt);

    % remove low FR units
    FR = fnd.FR({1, [-100 700]}); % (unit, trial)
    I = nanmean(FR, 2)>=1; % mean FR should >= 1 Hz
    fnd = fnd.set_unit_criteria('custom', I);
    fprintf('%d low FR units are removed.\n', sum(~I))

    % remove units changed over time
    I = ~isnan(fnd.misc.SNR);
    fnd = fnd.set_unit_criteria('custom', I);
    fprintf('%d units changed over time.\n', sum(~I))

    % save
    save(nfnd_fname, 'fnd')
end