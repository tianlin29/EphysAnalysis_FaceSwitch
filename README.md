Initialize.m
- 三个储存数据的地方：local、formatted data、dataserver本地（不用精确到project folder）
- format_panel参数预设

src/Common
- 重要：
- get_file_path
- set_plot_opt_2cond

- 不太重要：
- combineTrialData
- load_behavioral_data

src/PSTH
- showPopPSTH 存在更改
- showPopPSTH_choice_dependent
- plotPSTH single ??
- distinctPalette ??


*****************
TODO:
所有fnd的预处理：删除低FR的unit，删除SNR为NaN的unit，确保没有targ cho为NaN的trial


*****************
preproc/ 预处理
- Preprocess_*.m Gouki从FIRA转为FND

preproc_face_switch/
- Preprocess_ROME_FaceSwitch.m 处理main task; 处理training data
- Preprocess_ROME_FaceSwitchFixMulti.m 处理rand40

src/FiringRate/ (无用)
- showPopFR.m

src_face_switch/Common/
- get_file_path: 给定猴子和任务名称，加载各种数据

script_face_switch/learnTask2/
- unsigned choice






