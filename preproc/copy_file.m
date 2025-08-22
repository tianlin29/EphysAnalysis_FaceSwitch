run('../Initialize.m');
addpath C:\Engine\EphysPreprocess\FIRA2FND

monkey = 'Woody';
experiment = 'learnTask3';

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
