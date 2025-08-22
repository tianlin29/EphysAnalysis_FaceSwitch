clear;

%% define common path
MainDir = pwd;
PreprocDir = fullfile(MainDir, 'data', 'preproc'); % preprocessed rome data
MainFigDir = fullfile(MainDir, 'figures');
MainInterimDir = fullfile(MainDir, 'data', 'interim');

% choose formatted fnd data folder
FormattedDir_List = {
    'D:\LocalData\formatted_data\', ...
    '\\10.10.49.250\formatted_data\', ...
    '/mnt/raid-volume/formatted_data/', ... % when run code on dataserver
    };
if exist(FormattedDir_List{1}, 'dir')
    FormattedDir = FormattedDir_List{1};
elseif exist(FormattedDir_List{2}, 'dir')
    warning('loading data from dataserver')
    FormattedDir = FormattedDir_List{2};
elseif exist(FormattedDir_List{3}, 'dir')
    FormattedDir = FormattedDir_List{3};
else
    error('Check FormattedDir.')
end

%% add function path
path(pathdef); % reset to default path
addpath(genpath(fullfile(MainDir, 'src'))) % genpath adds all the subfolders

%% set default figure style
set(groot, 'defaultFigureColor', 'w')
set(groot, 'defaultFigurePaperPositionMode', 'auto')
set(groot, 'defaultAxesFontSize', 7)
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesTickDirMode', 'manual');
set(groot, 'defaultAxesTickLen', [.03 .03]);
set(groot, 'DefaultAxesXTickLabelRotation', 0);

%% suppress warning
warning('OFF', 'MATLAB:MKDIR:DirectoryExists');
warning('OFF', 'MATLAB:dispatcher:nameConflict');

%% plot
format_panel('mode','manuscript');
% opengl software

