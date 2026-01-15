function [file_path, n_files] = get_file_path(monkey, experiment, session, suffix, subsession, FormattedDir)

project = 'FaceSwitch_Monkey';

if nargin<6
    FormattedDir = 'D:\LocalData\formatted_data';
end
if nargin<5 || isempty(subsession)
    subsession = 1;
end
if nargin<4
    suffix = 'FND_sorted_preprocessed'; % could be eye, FIRA, FIRA_sorted, FND, FND_sorted, FND_sorted_preprocessed
end

switch monkey
    case 'Woody'
        session_list.memoryGuided = {'20230905'};
        session_list.noCue = {'20230923', '20230925', '20230926', '20231012', '20231013'}; % delete
        session_list.targ2stim_corrected = {'20240418'};
        session_list.rotation = {'20231024', '20231025'}; % delete
        session_list.learnTask3 = {'20231103', '20231104', '20231107', '20231108', '20231109', ...
            '20231110', '20231113', '20231114', '20231115', '20231116'}; % 10 sessions
        session_list.switchFreq = {'20231129', '20231130', '20231201', '20231202', '20231204'};
        session_list.passiveLong = {'20231212', '20231213', '20231214'}; % main face task is subsession 1, passive task is subsession 2
        session_list.intermediate = {'20231225', '20231226', '20231227', '20231228', '20240105'};
        session_list.learnTask4 = {'20240110', '20240111', '20240112', '20240113', '20240115'};
        session_list.threeExemplar = {'20240304', '20240305'}; % delete
        session_list.rand40 = {'20240329', '20240330', '20240401', '20240402'};
        session_list.RSVP_task1_3 = {'20240607', '20240608', '20240613'}; % main face task is subsession 1, passive task is subsession 2

    case 'Nick'
        session_list.memoryGuided = {'20250517'};
        session_list.learnTask2 = {'20250528', '20250529', '20250530', '20250531', '20250601', ...
            '20250604', '20250605', '20250606', '20250607', '20250608', ...
            '20250609', '20250610', '20250611', '20250612', '20250613'}; % 15 sessions
        session_list.learnTask2_increaseSwitchFreq = {'20250703'};
        session_list.learnTask3 = {'20250704', '20250705', '20250706', '20250707', '20250708', ...
            '20250709', '20250710', '20250711', '20250712', '20250713', ...
            '20250714', '20250715', '20250716'}; % these 3 sessions are actually increase switch frequency
        session_list.learnTask3_increaseSwitchFreq = {};
        session_list.learnTask4 = {'20250721', '20250722', '20250723', '20250724', '20250725', ...
            '20250726', '20250727', '20250728'}; % 8 sessions
        session_list.learnTask5 = {};
        session_list.passiveLong = {'20250807', '20250808', '20250809'}; % main face task is subsession 1, passive task is subsession 2
        session_list.RSVP_task1_2 = {'20250812', '20250813'}; % session 1: pure mask from cycle 2; main face task is subsession 1, passive task is subsession 2
        session_list.rand40 = {'20250815', '20250816', '20250817', '20250818'};
        session_list.faceColor = {'20250819', '20250820', '20250821', '20250822', '20250823', ...
            '20250825', '20250826', '20250827', '20250829', '20250830', ...
            '20250901', '20250903', '20250904', '20250905', '20250906', ...
            '20250908', '20250909', '20250910', '20250911', '20250912', ...
            '20250913', '20250914', '20250915'};
        session_list.threeExemplar = {'20250922', '20250923'};
        session_list.faceColor_passiveLong = {'20251019', '20251020', '20251021'};
end

n_files = length(session_list.(experiment));
if nargin<3 % no need to load file
    file_path = [];
else % load file
    if contains(suffix, '.')
        file_name = sprintf('%s%s_%02d_%s', monkey, session_list.(experiment){session}, subsession, suffix);
    else
        file_name = sprintf('%s%s_%02d_%s.mat', monkey, session_list.(experiment){session}, subsession, suffix);
    end
    file_path = fullfile(FormattedDir, project, monkey, session_list.(experiment){session}, file_name);
end

end