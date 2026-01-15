classdef runMDS
    properties
        def = struct(...
            'ID_method', 'stim_id', ...
            'time_win', {{{1, [50 100]}, {1, [100 150]}, {1, [150 200]}, {1, [200 250]}, ...
                         {2, [0 50]}, {2, [50 100]}, {2, [100 150]}, {2, [150 200]}, ...
                         {3, [-50 0]}, {3, [0 50]}, {3, [50 100]}, {3, [100 150]}}}, ...
            'plot', set_plot_opt('vik', 40), ...
            'pcscore_human', [], ...
            'pcscore_monkey', []);
    end

    methods
        %% define options
        function run_mds = runMDS(opt) 
            if nargin > 0
                run_mds.def = safeStructAssign(run_mds.def, opt);
            end
        end

        %% getRDM
        function [fh_RDM, fh_MDS, stat] = getRDM(run_mds, fnd)
            opt = run_mds.def;

            % define ID
            switch opt.ID_method
                case 'stim_id'
                    ID = fnd.getp('stim_id');

                case 'cat_1'
                    ID = fnd.getp('stim_id');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=1) = NaN;

                case 'cat_2'
                    ID = fnd.getp('stim_id');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=2) = NaN;
                    ID = ID-20;

                case 'coh'
                    ID = get_ID(fnd, [], 'coh');

                case 'coh_cat_1'
                    ID = get_ID(fnd, [], 'coh');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=1) = NaN;
                    ID = ID - (min(ID(:))-1);

                case 'coh_cat_2'
                    ID = get_ID(fnd, [], 'coh');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=2) = NaN;

                case 'coh_task'
                    ID = get_ID(fnd, [], 'coh');
                    task_set = fnd.getp('task_set');
                    ID(task_set==2) = ID(task_set==2) + max(ID(:));

                case 'coh_task_correct'
                    ID = get_ID(fnd, [], 'coh');
                    corr = fnd.getp('targ_cor')==fnd.getp('targ_cho');
                    targ_cor = fnd.getp('targ_cor');
                    ID(~corr | targ_cor~=1) = NaN;
                    ID = ID - (min(ID(:))-1);
                    task_set = fnd.getp('task_set');
                    ID(task_set==2) = ID(task_set==2) + max(ID(:));
                    
            end

            if min(ID(:))~=1
                error('check ID');
            end

            % get dissimilarity between stimuli, then perform MDS on dissimilarity matrix
            ntime = length(opt.time_win);
            [FR, RSM, RDM, MDS_2d] = deal(cell(1, ntime));
            for t = 1:ntime
                FR{t} = fnd.FR(opt.time_win{t}, ID, true); % (unit, condition)
                RSM{t} = corrcoef(FR{t}); % correlation between conditions
                RDM{t} = 1-RSM{t}; % dissimilarity matrix
                MDS_2d{t} = mdscale(RDM{t}, 2, 'Criterion', 'metricstress'); % get the first 2 dims
            end

            stat.FR = FR; % (unit, condition)
            stat.RSM = RSM; % (condition, condition)
            stat.RDM = RDM; % (condition, condition)
            stat.MDS_2d = MDS_2d; % (condition, 2)

            % plot dissimilarity matrix (actually similarity matrix)
            fh_RDM = run_mds.plotRSM(stat);

            % plot MDS dots
            fh_MDS = run_mds.plotMDS(stat);
        end

        %% getRDM_PSTH
        function [fh_RDM, fh_MDS, stat] = getRDM_PSTH(run_mds, fnd)
            opt = run_mds.def;

            % define ID
            switch opt.ID_method
                case 'stim_id'
                    ID = fnd.getp('stim_id');

                case 'cat_1'
                    ID = fnd.getp('stim_id');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=1) = NaN;

                case 'cat_2'
                    ID = fnd.getp('stim_id');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=2) = NaN;
                    ID = ID-20;

                case 'coh'
                    ID = get_ID(fnd, [], 'coh');

                case 'coh_cat_1'
                    ID = get_ID(fnd, [], 'coh');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=1) = NaN;
                    ID = ID - (min(ID(:))-1);

                case 'coh_cat_2'
                    ID = get_ID(fnd, [], 'coh');
                    targ_cor = fnd.getp('targ_cor');
                    ID(targ_cor~=2) = NaN;

                case 'coh_task'
                    ID = get_ID(fnd, [], 'coh');
                    task_set = fnd.getp('task_set');
                    ID(task_set==2) = ID(task_set==2) + max(ID(:));

                case 'coh_task_correct'
                    ID = get_ID(fnd, [], 'coh');
                    corr = fnd.getp('targ_cor')==fnd.getp('targ_cho');
                    targ_cor = fnd.getp('targ_cor');
                    ID(~corr | targ_cor~=1) = NaN;
                    ID = ID - (min(ID(:))-1);
                    task_set = fnd.getp('task_set');
                    ID(task_set==2) = ID(task_set==2) + max(ID(:));
                    
            end

            if min(ID(:))~=1
                error('check ID');
            end

            % prepare PSTH
            psth = fnd.PSTH(ID, {'boxcar', 100}, [], true); % detrended
            tstamp = fnd.tstamp;

            all_psth = cat(2, psth{1}(:,tstamp{1}>=100 & tstamp{1}<400,:), ...
                psth{2}(:,tstamp{2}>=0 & tstamp{2}<300,:), ...
                psth{3}(:,tstamp{3}>=-300 & tstamp{3}<0,:));

            % get dissimilarity between stimuli, then perform MDS on dissimilarity matrix
            ntime = size(all_psth, 2);
            [FR, RSM, RDM, MDS_2d] = deal(cell(1, ntime));
            wb = waitbar_text(0);
            for t = 1:ntime
                waitbar_text(t/ntime, wb);
                FR{t} = squeeze(all_psth(:,t,:)); % (unit, condition)
                RSM{t} = corrcoef(FR{t}); % correlation between conditions
                RDM{t} = 1-RSM{t}; % dissimilarity matrix
                MDS_2d{t} = mdscale(RDM{t}, 2, 'Criterion', 'metricstress'); % get the first 2 dims
            end
            waitbar_text('close', wb);

            stat.FR = FR; % (unit, condition)
            stat.RSM = RSM; % (condition, condition)
            stat.RDM = RDM; % (condition, condition)
            stat.MDS_2d = MDS_2d; % (condition, 2)

            % plot dissimilarity matrix (actually similarity matrix)
            fh_RDM = run_mds.plotRSM(stat);

            % plot MDS dots
            fh_MDS = run_mds.plotMDS(stat);
        end
  
        %% plotFR
        function fh_FR = plotFR(run_mds, FR)
            % plot FR
            opt = run_mds.def;
            epoch_name_list = {'Stim on', 'Stim off', 'FP off', 'Resp'};
            [redblue, clims] = defineColor();
            ntime = length(FR);

            fh_FR = figure('Position', [50 300 100*ntime 100]);
            for t = 1:ntime
                subplot(1,ntime,t); hold on; colormap(redblue);
                imagesc(FR{t}, clims); axis ij; axis equal; xlim([1 size(FR{1},2)]); ylim([1 size(FR{1},1)]);
                title({epoch_name_list{opt.time_win{t}{1}}, sprintf('%d~%d ms', opt.time_win{t}{2}(1), opt.time_win{t}{2}(2))});
            end
            % colorbar;
            format_panel(fh_FR, 'xlabel', 'stimulus ID', 'ylabel', 'unit', ...
                'label_clean', {0,0,0}, 'ticklabel_clean', {0,0,0})
        end

        %% plotFRDist
        function fh = plotFRDist(run_mds, FR_dist)
            % plot FR
            opt = run_mds.def;
            [redblue, ~] = defineColor();
            ntime = length(FR_dist);

            fh = figure('Position', [50 300 130*ntime 200]);
            for t = 1:ntime
                subplot(1,ntime,t); hold on; colormap(redblue);
                imagesc(FR_dist{t}); axis ij; axis equal; xlim([0 size(FR_dist{1},2)]); ylim([0 size(FR_dist{1},1)]);
                title(sprintf('epoch %d: %d~%d ms', opt.time_win{t}{1}, opt.time_win{t}{2}(1), opt.time_win{t}{2}(2)));
            end
            colorbar;
            format_panel(fh, 'xlabel', 'stim ID', 'ylabel', 'stim ID', ...
                'label_clean', {0,0,0}, 'ticklabel_clean', {0,0,0})
        end

        %% plotRSM
        function fh_RSM = plotRSM(run_mds, stat)
            % plot dissimilarity matrix (actually similarity matrix)
            RSM = stat.RSM;
            opt = run_mds.def;
            epoch_name_list = {'Stim on', 'Stim off', 'FP off', 'Resp'};
            [redblue, clims] = defineColor();
            ntime = length(RSM); ncond = size(RSM{1},1);

            fh_RSM = figure('Position', [50 300 100*ntime 100]);
            for t = 1:ntime
                subplot(1,ntime,t); hold on; colormap(redblue);
                imagesc(RSM{t}, clims); axis ij; axis equal; xlim([1 ncond]); ylim([1 ncond]);
                title({epoch_name_list{opt.time_win{t}{1}}, sprintf('%d~%d ms', opt.time_win{t}{2}(1), opt.time_win{t}{2}(2))});
            end
            % colorbar;
            format_panel(fh_RSM, 'xlabel', 'stimulus ID', 'ylabel', 'stimulus ID', ...
                'label_clean', {0,0,0}, 'ticklabel_clean', {0,0,0})
        end

        %% plotMDS
        function fh_MDS = plotMDS(run_mds, stat)
            MDS_2d = stat.MDS_2d;
            opt = run_mds.def;
            epoch_name_list = {'Stim on', 'Stim off', 'FP off', 'Resp'};
            [redblue, ~] = defineColor();
            ntime = length(MDS_2d); ncond = size(MDS_2d{1},1);

            fh_MDS = figure('Position', [50 100 100*ntime 100]);
            for t = 1:ntime
                subplot(1,length(opt.time_win),t); hold on; colormap(redblue);
                for i = 1:ncond
                    scatter(MDS_2d{t}(i,1), MDS_2d{t}(i,2), 8, opt.plot.color(i,:), 'filled');
                end
                title({epoch_name_list{opt.time_win{t}{1}}, sprintf('%d~%d ms', opt.time_win{t}{2}(1), opt.time_win{t}{2}(2))});
            end
            format_panel(fh_MDS, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', ...
                'label_clean', {0,0,0}, 'ticklabel_clean', {0,0,0})
        end

        %% getCorrFacePC
        function [fh_RDM_facePC, fh_corr, fh_pval] = getCorrFacePC(run_mds, RDM)
            [redblue, clims] = defineColor();
            opt = run_mds.def;

            % get RDM of face PC
            if contains(opt.ID_method, '1') % human category
                RDM_facePC = corrcoef(opt.pcscore_human'); % correlation between conditions
            elseif contains(opt.ID_method, '2') % monkey category
                RDM_facePC = corrcoef(opt.pcscore_monkey');
            end

            % get correlation between timepoints
            [R, P] = deal(1, length(RDM));
            for i = 1:length(RDM)
                [r_, p_] = mantel_test(RDM{i}, RDM_facePC, 1000, 'both');
                R(i) = r_;
                P(i) = p_;
            end

            % plot
            fh_RDM_facePC = figure('Position', [50 100 150 150]); colormap(redblue);
            imagesc(RDM_facePC, clims); colorbar; title('Face RDM')
            format_panel(fh_RDM_facePC, 'xlabel', 'stimulus ID', 'ylabel', 'stimulus ID', ...
                'label_clean', {0,0,0}, 'ticklabel_clean', {0,0,0})
            xtickangle(0)

            fh_corr = figure('Position', [50 100 150 150]); colormap(redblue);
            imagesc(R, clims); colorbar; title({'Correlation between', 'neural RDM and face RDM'})
            format_panel(fh_corr, 'xlabel', 'Time window', 'ylim', xlim)
            xtickangle(0)

            fh_pval = figure('Position', [300 100 150 150]); colormap(redblue);
            imagesc(P<0.05); colorbar; title('P<0.05 (Mantel test)')
            format_panel(fh_pval, 'xlabel', 'Time window', 'ylim', xlim)
            xtickangle(0)
        end



        %% getCorrRDM
        function [fh_corr, fh_pval] = getCorrRDM(run_mds, RDM)
            [redblue, clims] = defineColor();

            % get correlation between timepoints
            [R, P] = deal(nan(length(RDM), length(RDM)));
            for i = 1:length(RDM)
                for j = i:length(RDM)
                    [r_, p_] = mantel_test(RDM{i}, RDM{j}, 1000, 'both');
                    R(i,j) = r_; R(j,i) = r_;
                    P(i,j) = p_; P(j,i) = p_;
                end
            end

            % plot
            fh_corr = figure('Position', [50 100 150 150]); colormap(redblue);
            imagesc(R, clims); colorbar; title('Correlation of RDM')
            format_panel(fh_corr, 'xlabel', 'Time window', 'ylabel', 'Time window')
            xtickangle(0)

            fh_pval = figure('Position', [300 100 150 150]); colormap(redblue);
            imagesc(P<0.05); colorbar; title('P<0.05 (Mantel test)')
            format_panel(fh_pval, 'xlabel', 'Time window', 'ylabel', 'Time window')
            xtickangle(0)
        end

        %% getCorrFR
        function [fh_corr, fh_pval] = getCorrFR(run_mds, FR)
            [redblue, clims] = defineColor();

            % get correlation between timepoints
            [R, P] = deal(nan(length(FR), length(FR)));
            for i = 1:length(FR)
                for j = i:length(FR)
                    [r_, p_] = mantel_test(FR{i}, FR{j}, 1000, 'both');
                    R(i,j) = r_; R(j,i) = r_;
                    P(i,j) = p_; P(j,i) = p_;
                end
            end

            % plot
            fh_corr = figure('Position', [50 100 150 150]); colormap(redblue);
            imagesc(R, clims); colorbar; title('Correlation of FR')
            format_panel(fh_corr, 'xlabel', 'Time window', 'ylabel', 'Time window')
            xtickangle(0)
            
            fh_pval = figure('Position', [300 100 150 150]); colormap(redblue);
            imagesc(P<0.05); colorbar; title('P<0.05 (Mantel test)')
            format_panel(fh_pval, 'xlabel', 'Time window', 'ylabel', 'Time window')
            xtickangle(0)
        end

        %% plotImage
        function plotImage(run_mds, mds_2d, FigDir, category)
            fname = load('C:\Engine\rome_external\PassiveFixationMask\randomFace40\process.mat').fname;
            ImgPath = 'C:\Engine\rome_external\PassiveFixationMask\randomFace40\stim\';
            opt = run_mds.def;

            if category==2
                tmp = 20;
            else
                tmp = 0;
            end
            ncond = size(mds_2d{1},1);
            for t = 1:length(opt.time_win)
                fh = figure('Position', [50 100 300 300]); hold on;
                colormap(gray);
                for i = 1:ncond
                    % load image
                    img = imread([ImgPath, fname{i+tmp}]);

                    % determine width and height
                    img_width = size(img, 2);
                    img_height = size(img, 1);
                    scale = 0.0006; % 0.0006 for z-scored; 0.000015 for un-z-scored
                    width = img_width * scale;
                    height = img_height * scale;

                    % plot image at location
                    image([mds_2d{t}(i,1)-width/2, mds_2d{t}(i,1)+width/2], [mds_2d{t}(i,2)-height/2, mds_2d{t}(i,2)+height/2], flipud(img), 'CDataMapping', 'scaled');
                end
                title(sprintf('epoch %d: %d~%d ms', opt.time_win{t}{1}, opt.time_win{t}{2}(1), opt.time_win{t}{2}(2)));
                format_panel(fh, 'xlabel', 'Dim 1', 'ylabel', 'Dim 2', 'xlim')
                axis equal
                print(fh, '-dpdf', [FigDir, sprintf('zscored_images_epoch%d_%d_%d_ms.pdf', opt.time_win{t}{1}, opt.time_win{t}{2}(1), opt.time_win{t}{2}(2))]);
            end
        end
    end
end

%% define imagesc color
function [redblue, clims] = defineColor()
m = 256;
redblue = [linspace(0,1,m/2)', linspace(0,1,m/2)', ones(m/2,1); ...
    ones(m/2,1), linspace(1,0,m/2)', linspace(1,0,m/2)'];
clims = [-1 1];
end

%% Mantel test between two dissimilarity matrices
function [r_obs, p_value] = mantel_test(D1, D2, n_perm, tail)
    % Input: 
    %   D1, D2: Two n×n symmetric distance matrices
    %   n_perm: Number of permutations (default 1000)
    %   tail: Test type 'both' (two-tailed), 'right' (positive correlation), 'left' (negative correlation)
    % Output:
    %   r_obs: Observed correlation coefficient
    %   p_value: Significance p-value

    if nargin < 3, n_perm = 1000; end
    if nargin < 4, tail = 'both'; end

    if size(D1,1)~=size(D1,2)
        % in case of FR (unit, condition), not a square matrix
        vec_D1 = D1(:);
        vec_D2 = D2(:);

        % calculate original correlation coefficient (Pearson)
        r_obs = corr(vec_D1, vec_D2, 'Type', 'Pearson');

        % permutation test
        [m, n] = size(D1);
        r_perm = zeros(n_perm, 1);
        for i = 1:n_perm
            m_perm_idx = randperm(m);
            n_perm_idx = randperm(n);
            D2_perm = D2(m_perm_idx, n_perm_idx); % permute both rows and columns
            vec_D2_perm = D2_perm(:);
            r_perm(i) = corr(vec_D1, vec_D2_perm, 'Type', 'Pearson');
        end
    else
        % in case of RDM (condition, condition), extract lower triangular parts (excluding diagonal)
        vec_D1 = D1(tril(true(size(D1)), -1));
        vec_D2 = D2(tril(true(size(D2)), -1));

        % calculate original correlation coefficient (Pearson)
        r_obs = corr(vec_D1, vec_D2, 'Type', 'Pearson');

        % permutation test
        n = size(D1, 1);
        r_perm = zeros(n_perm, 1);
        for i = 1:n_perm
            perm_idx = randperm(n);
            D2_perm = D2(perm_idx, perm_idx); % permute both rows and columns at the same time
            vec_D2_perm = D2_perm(tril(true(size(D2_perm)), -1));
            r_perm(i) = corr(vec_D1, vec_D2_perm, 'Type', 'Pearson');
        end
    end

    % calculate p-value
    switch tail
        case 'right'
            p_value = sum(r_perm >= r_obs) / n_perm;
        case 'left'
            p_value = sum(r_perm <= r_obs) / n_perm;
        case 'both'
            p_value = sum(abs(r_perm) >= abs(r_obs)) / n_perm;
    end
end