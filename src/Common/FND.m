classdef FND < matlab.mixin.Copyable
% flexible neural data format (class)
%   Copyable class alows copy function

    %% public properties
    properties (Access = public)
        % meta data
        alignto % 1 x nepoch structu
        unitID % unit x 3 [session, channel, unit]
        unit_type % hash, multi, single
        unit_quality % 1.. high quality, 3.. low quality, NaN .. hash or multi
        monkey % unit x 1 cell
        area % unit x 1 cell
        valid_unit % unit x 1 logical
        unit_criteria_description
        session_name % nsession x 1 cell
        tstamp % 1 x nepoch cell
        description
        misc % misc.(param_name): unit x 1 cell (other neuron info)
        info % variable reserved to put any info
        session_info % info about the recorded session (if the data is only one session)
        
        % parameter
        param % param.(param_name): neuron x trial (cell or double)
        ntrial % number of trial: unit x 1
        
    end
    
    %% protected properties
    properties (SetAccess = protected)
        % meta data
        date
        version
        compression_type % raw, zip
        
        % spike data
        data % data{n} = neuron x time x trial
             % this is uint8 format to reduce memory load
             % 1 .. spike, 0 .. no spike, 255 .. no data (NaN)
        PSTH_trend % PSTH_trend{n} = neuron x time
             % used to detrend PSTH if required
             % set_detrend_PSTH set this value
             % if not set, PSTH function will just run set_detrend_PSTH
    end
    
    %% @@@@ main methods @@@@
    methods
        %% constructor
        function fnd = FND(data, opt)
            % constructer function.
            % data should follow FND data format (see above) but it could
            % be double (1.. spike, 0.. no spike, NaN.. no data)
            % opt can include all other parameters in FND (e.g., opt.param, opt. unitID) 
            
            % copy data in opt to fnd
            fi = fieldnames(opt);
            for n=1:length(fi)
                fnd.(fi{n}) = opt.(fi{n});
            end

            % compress spike data
            if isempty(fnd.compression_type)
                fnd.compression_type = 'raw';
            end
            fnd.compress_data(data);
            
            % set unit info
            nunit = size(data{1},1);
            ntri = size(data{1},3);
            fnd.valid_unit = true(nunit,1);
            fnd.unit_criteria_description = 'all';
            
            if isempty(fnd.description)
                fnd.description = 'no info';
            end
            
            % validating information (if the number of info matches data)
            fnd.unitID = fnd.check_unit_info(fnd.unitID, nunit, 3, 'unitID');
            fnd.unit_type = fnd.check_unit_info(fnd.unit_type, nunit, 1, 'unit_type');
            fnd.unit_quality = fnd.check_unit_info(fnd.unit_quality, nunit, 1, 'unit_quality');
            fnd.monkey = fnd.check_unit_info(fnd.monkey, nunit, 1, 'monkey');
            fnd.area = fnd.check_unit_info(fnd.area, nunit, 1, 'area');
            fnd.ntrial = fnd.check_unit_info(fnd.ntrial, nunit, 1, 'ntrial');
            
            if isstruct(fnd.param)
                fi = fieldnames(fnd.param);
                for n=1:length(fi)
                    if size(fnd.param.(fi{n}),1) ~= nunit || size(fnd.param.(fi{n}),2) ~= ntri
                        error('param %s do not match spike data.', fi{n});
                    end
                end
            end
            if isstruct(fnd.misc)
                fi = fieldnames(fnd.misc);
                for n=1:length(fi)
                    if length(fnd.misc.(fi{n})) ~= nunit
                        error('misc.%s do not match spike data.', fi{n});
                    end
                end
            end
            if ~isempty(fnd.alignto)
                if length(fnd.alignto) ~= length(data)
                    error('alignto and spike data do not have the same epochs.');
                end
            end
            if ~isempty(fnd.tstamp)
                if length(fnd.tstamp) ~= length(data)
                    error('tstamp and spike data do not have the same epochs.');
                end
            end
            
            % set trial num
            if any(isnan(fnd.ntrial))
                warning('ntrial is deficient. filling automatically.');
                for n=1:length(fnd.ntrial)
                    if isnan(fnd.ntrial(n))
                        mxtri = 0;
                        for a=1:length(fnd.data)
                            d = permute(fnd.data{a}(n,:,:), [2 3 1]);
                            t = find(any(d~=255,1), 1, 'last');
                            if ~isempty(t)
                                mxtri = max(mxtri, t);
                            end
                        end
                        fnd.ntrial(n) = mxtri;
                    end
                end
            end
            
            % FND format info
            fnd.date = datestr(now);
            fnd.version = '1.1';
        end
        
        %% set_unit_criteria
        function nfnd = set_unit_criteria(fnd, varargin)
            % set unit selection criteria
            % fnd.set_unit_criteria('quality', 'no hash', 'monkey', 'IS', ..)
            %   -> valid_unit option will be modified by this function.
            %      once criteria are set, other FND functions (e.g. PSTH)
            %      return only neural data of valid units
            %      everytime you set criteria, only the criteria you put
            %      will be at work (previous criteria will be cleared)
            % fnd_new = fnd.set_unit_criteria('quality', ...
            %   -> if output is specified, it will return a new FND format
            %      (copy), which only contains valid unit (invalidated unit
            %      is unrecoverable from the returned fnd_new).
            %
            %   'append': if 'off'[default] previous criteria are cleared. 'on' or 'off' 
            %   'quality': no hash, single, hash only
            %   'monkey': monkey name
            %   'area': area name
            %   'percent missing trial': percentage of missing trials
            %               (some units are recorded only in part of trials)
            %   'custom': custom criterion (put N x 1 logical)
            %   'customID': choose custom unit by unitID
            %               (N x 3 matrix; [sessionID, ch, unitID])
            prev_u = fnd.valid_unit;
            fnd.valid_unit(:) = true;
            fnd.unit_criteria_description = '';
            for n1=1:nargin/2
                switch varargin{n1*2-1}
                    case 'append'
                        switch varargin{n1*2}
                            case 'on'
                                fnd.valid_unit = fnd.valid_unit & prev_u;
                            case 'off'
                                % do no thing
                            otherwise
                                error('no such append option type: %s', varargin{n1*2});
                        end
                    case 'quality'
                        switch varargin{n1*2}
                            case 'no hash'
                                fnd.valid_unit = fnd.valid_unit & ~strcmp(fnd.unit_type, 'hash');
                                fnd.unit_criteria_description = [fnd.unit_criteria_description, 'no hash, '];
                            case 'single'
                                fnd.valid_unit = fnd.valid_unit & ismember(fnd.unit_quality, [1,2,3]);
                                fnd.unit_criteria_description = [fnd.unit_criteria_description, 'single, '];
                            case 'hash only'
                                fnd.valid_unit = fnd.valid_unit & strcmp(fnd.unit_type, 'hash');
                                fnd.unit_criteria_description = [fnd.unit_criteria_description, 'hash only, '];
                            otherwise
                                error('no such quality type: %s', varargin{n1*2});
                        end
                    case 'monkey'
                        fnd.valid_unit = fnd.valid_unit & strcmp(fnd.monkey, varargin{n1*2});
                        fnd.unit_criteria_description = [fnd.unit_criteria_description, varargin{n1*2} ', '];
                    case 'area'
                        fnd.valid_unit = fnd.valid_unit & strcmp(fnd.area, varargin{n1*2});
                        fnd.unit_criteria_description = [fnd.unit_criteria_description, varargin{n1*2} ', '];
                    case 'percent missing trial'
                        flg = false(size(fnd.data{1},1), size(fnd.data{1},3));
                        for a=1:length(fnd.data)
                            flg = flg | permute(any(fnd.data{a}~=255,2), [1 3 2]);
                        end
                        rec_ntri = sum(flg,2);
                        per_missing = (fnd.ntrial - rec_ntri)./fnd.ntrial * 1e2;
                        fnd.valid_unit = fnd.valid_unit & per_missing < varargin{n1*2};
                        fnd.unit_criteria_description = [fnd.unit_criteria_description, 'small missing trials, '];
                    case 'custom'
                        fnd.valid_unit = fnd.valid_unit & varargin{n1*2};
                        fnd.unit_criteria_description = [fnd.unit_criteria_description, 'custom, '];
                    case 'session'
                        ID = varargin{n1*2};
                        fnd.valid_unit = fnd.valid_unit & ismember(fnd.unitID(:,1), ID);
                        fnd.unit_criteria_description = [fnd.unit_criteria_description, 'session, '];
                    case 'customID'
                        ID = varargin{n1*2};
                        fnd.valid_unit = fnd.valid_unit & ismember(fnd.unitID, ID, 'rows');
                        fnd.unit_criteria_description = [fnd.unit_criteria_description, 'customID, '];
                    otherwise
                        error('no such criterion: %s', varargin{n1*2-1});
                end
            end
            if nargout > 0
                % if output is specified. A new fnd with the specified
                % selection criteria will be returned. The old fnd does not
                % apply the criteria (remain unchanged).
                nfnd = fnd.extract_valid_unit;
                fnd.set_unit_criteria('custom', prev_u);
            end
        end
        
        %% extract_trial
        function nfnd = extract_trial(fnd, trialFlg, rmv_empty_unit)
            % extract specified trials from fnd and return a new fnd
            %   trialFlg = unit x trial flag (1 or 0)
            fnd = fnd.extract_valid_unit;
            if size(fnd.data{1},1) ~= size(trialFlg,1) || size(fnd.data{1},3) ~= size(trialFlg,2)
                error('trialFlg size does not match fnd data size');
            end
            common_trial = all(all(trialFlg == trialFlg(1,:),1), 2); % if all trials are common across unit
            trialFlg = logical(trialFlg);
            maxtr = max(sum(trialFlg,2));
            nfnd = copy(fnd);
            p = fieldnames(nfnd.param);
            for n=1:length(p)
                if common_trial
                    nfnd.param.(p{n})(:, ~trialFlg(1,:)) = [];
                else
                    pm = nfnd.param.(p{n});
                    if iscell(pm)
                        npm = cell(size(pm,1), maxtr);
                    else
                        npm = nan(size(pm,1), maxtr);
                    end
                    for u=1:size(trialFlg,1)
                        npm(u,1:sum(trialFlg(u,:))) = pm(u, trialFlg(u,:));
                    end
                    nfnd.param.(p{n}) = npm;
                end
            end
            
            for a=1:length(nfnd.data)
                if common_trial
                    nfnd.data{a}(:,:, ~trialFlg(1,:)) = [];
                else
                    d1 = nfnd.data{a};
                    nd1 = ones([size(d1,1), size(d1,2), maxtr], 'uint8') * 255;
                    for u=1:size(trialFlg,1)
                        nd1(u,:,1:sum(trialFlg(u,:))) = d1(u,:, trialFlg(u,:));
                    end
                    nfnd.data{a} = nd1;
                end
            end
            nfnd.ntrial = sum(trialFlg,2);
            
            if ~exist('rmv_empty_unit', 'var') || rmv_empty_unit
                nfnd = set_unit_criteria(nfnd, 'custom', nfnd.ntrial~=0);
            end
        end
        
        %% extract_epoch
        function nfnd = extract_epoch(fnd, epoch)
            % extract specified epoch from fnd and return a new fnd
            %   epoch = vector
            if length(fnd.data) < max(epoch)
                error('There is no epoch %d', max(epoch));
            end
            nfnd = copy(fnd);
            nfnd.alignto = nfnd.alignto(epoch);
            nfnd.tstamp= nfnd.tstamp(epoch);
            nfnd.data = nfnd.data(epoch);
            if ~isempty(nfnd.PSTH_trend)
                nfnd.PSTH_trend = nfnd.PSTH_trend(epoch);
            end
        end
        
        %% getp (get parameter)
        function out = getp(fnd, param_name, rmv_redundant_trial)
            % out = getp(fnd, param_name)
            %   -> get parameter value
            %      example: out = fnd.getp('targ_cho');
            %      it will return only the parameters of valid units
            %      if you instead use fnd.param.targ_cho, you will get
            %      info from all the units, so should be careful.
            if ~exist('param_name', 'var')
                if isstruct(fnd.param)
                    fprintf('Parameter:\n');
                    fi = fieldnames(fnd.param);
                    for n=1:length(fi)
                        fprintf('\t%s\n', fi{n});
                    end
                end
                return;
            end
            
            fnd = fnd.extract_valid_unit;
            if contains(param_name, 'misc.')
                param_name = strrep(param_name, 'misc.', '');
                out = fnd.misc.(param_name);
                
            else
                out = fnd.param.(param_name);
            end
            
            if ~exist('rmv_redundant_trial', 'var') || rmv_redundant_trial
                if size(out,2) > max(fnd.ntrial)
                    out = out(:, 1:max(fnd.ntrial));
                end
            end
        end
        
        function param = getp_all(fnd, rmv_redundant_trial)
            % out = fnd.getp_all()
            %   get all parameters
            fnd = fnd.extract_valid_unit;
            param = fnd.param;
            if ~exist('rmv_redundant_trial', 'var') || rmv_redundant_trial
                fi = fieldnames(param);
                for n=1:length(fi)
                    if size(param.(fi{n}),2) > max(fnd.ntrial)
                        param.(fi{n}) = param.(fi{n})(:, 1:max(fnd.ntrial));
                    end
                end
            end            
        end
        
        %% setp (set parameter)
        function setp(fnd, param_name, param1)
            % setp(fnd, param_name, param)
            %   -> set parameter value
            %      example: out = fnd.setp('targ_cho', targ_cho);
            %      it will set parameters for all units
            if isscalar(param1)
                fnd.param.(param_name) = param1 * ones(size(fnd.data{1},1), size(fnd.data{1},3));
                return;
            end
            if size(param1,1) ~= size(fnd.data{1},1)
                error('Input parameters has invalid unit count');
            end
            if size(param1,2) ~= size(fnd.data{1},3)
                error('Input parameters has invalid trial count');
            end
            fnd.param.(param_name) = param1;
        end
        
        %% raster
        function r = raster(fnd, epochs, rmv_redundant_trial)
            % r = raster(fnd, epochs, rmv_redundant_trial)
            %   -> return raster (1.. spk, 0.. no spk, NaN.. no data)
            % r = 1 x nepoch cell
            % r{1} = unit x time x trial double
            fnd = fnd.extract_valid_unit;
            
            if ~exist('epochs', 'var') || isempty(epochs)
                epochs = 1:length(fnd.data);
            end
            nepoch = length(epochs);
            r = cell(1, nepoch);
            for n=1:nepoch
                r{n} = double(fnd.data{epochs(n)});
                r{n}(r{n}==255) = NaN;
            end
            
            if ~exist('rmv_redundant_trial', 'var') || rmv_redundant_trial
                if size(r{1},3) > max(fnd.ntrial)
                    for n=1:nepoch
                        r{n} = r{n}(:,:,1:max(fnd.ntrial));
                    end
                end
            end
            
        end
        
        %% raster_cond
        function r = raster_cond(fnd, condID, epochs, format)
            % r = raster_cond(fnd, condID)
            %   -> return raster for each condition
            % r = 1 x nepoch cell
            % format: 'cell' (default)
            %   r{1} = 1 x ncond cell
            %   r{1}{1} = unit x trial x time double
            % format: 'array'
            %   r{1} = unit x trial x time x condition double
            %
            % condID can be generated using trial_classifier.m
            fnd = fnd.extract_valid_unit;
            ncell = size(fnd.data{1},1);
            ncond = nanmax(condID(:));
            
            if ~exist('epochs', 'var') || isempty(epochs)
                epochs = 1:length(fnd.data);
            end
            nepoch = length(epochs);
            
            r = cell(1, nepoch);
            wb = waitbar_text(0);
            for n=1:nepoch
                if exist('format', 'var') && strcmp(format, 'array')
                    maxtr = 0;
                    for c=1:ncond
                        maxtr = max(maxtr, max(sum(condID==c,2)));
                    end
                    r{n} = nan(ncell, maxtr, size(fnd.data{epochs(n)},2), ncond);
                else
                    r{n} = cell(1, ncond);
                end
                for c = 1:ncond
                    waitbar_text((c + (n-1)*ncond)/(nepoch * ncond), wb);
                    if exist('format', 'var') && strcmp(format, 'array')
                        for u = 1:ncell
                            raster = permute(double(fnd.data{epochs(n)}(u,:,condID(u,:)==c)), [2 3 1]);
                            raster(raster==255) = NaN;
                            r{n}(u, 1:size(raster,2), :, c) = raster';
                        end
                    else
                        r{n}{c} = nan(ncell, max(sum(condID==c,2)), size(fnd.data{epochs(n)},2));
                        for u = 1:ncell
                            raster = permute(double(fnd.data{epochs(n)}(u,:,condID(u,:)==c)), [2 3 1]);
                            raster(raster==255) = NaN;
                            r{n}{c}(u, 1:size(raster,2), :) = raster';
                        end
                    end
                end
            end
            waitbar_text('close', wb);
        end
        
        %% PSTH
        function r = PSTH(fnd, condID, conv, cutoff, detrend, epoch)
            % r = PSTH(fnd, condID, conv, cutoff, detrend)
            %   -> return PSTH
            % r = 1 x nepoch cell
            % r{1} = unit x time x condition double (spk/s)
            %
            % condID can be generated using trial_classifier.m
            %    if empty ([]), average all trials
            % conv: convolusion kernel (default: [])
            %       {'boxcar', 100}: boxcar function with 100ms width
            %       {'gaussian', 20}: gaussian kernel with 20ms SD
            % cutoff: if you want to cut off time if less than 50% of
            %         trials contribute (default: false)
            % detrend: if you want to detrend PSTH (default: false)
            % epoch: if you want to analyze only part of the epochs (default: [])
            %        [1, 2]: caculate PSTH of epoch 1 and 2
            fnd = fnd.extract_valid_unit;
            if ~exist('conv', 'var')
                conv = [];
            end
            if ~exist('cutoff', 'var')
                cutoff = false; % no cutoff
            end
            if ~exist('detrend', 'var')
                detrend = false;
            end
            if ~exist('epoch', 'var')
                epoch = [];
            end
            if detrend
                fnd.set_detrend_PSTH();
            end
            
            ncell = size(fnd.data{1},1);
            trial_num = size(fnd.data{1},3);
            if ~exist('condID', 'var') || isempty(condID)
                condID = ones(ncell, trial_num);
            end
            
            ncond = nanmax(condID(:));
            nepoch = length(fnd.data);
            
            r = cell(1, nepoch);
            wb = waitbar_text(0);
            for n=1:nepoch
                if ~isempty(epoch) && ~ismember(n, epoch)
                    continue
                end
                
                r{n} = nan(ncell, size(fnd.data{n},2), ncond);
                for u = 1:ncell
                    waitbar_text((u + (n-1)*ncell)/(nepoch * ncell), wb);
                    raster = squeeze(double(fnd.data{n}(u,:,:)));
                    raster(raster==255) = NaN;
                    for c = 1:ncond
                            % convert spike flag to firing rate
                        r{n}(u,:,c) = nanmean(raster(:,condID(u,:)==c),2) * 1e3;
                        if detrend
                            r{n}(u,:,c) = r{n}(u,:,c) - fnd.PSTH_trend{n}(u,:);
                        end
                        if cutoff
                            % cut off if less than half of trials contribute
                            ind = mean(isnan(raster(:,condID(u,:)==c)),2) > .5;
                            r{n}(u,ind,c) = NaN;
                        end
                        if ~isempty(conv)
                            r{n}(u,:,c) = fnd.nanconv_data(r{n}(u,:,c), conv);
                        end
                    end
                end
            end
            waitbar_text('close', wb);
        end
        
        %% cutoff
        function co = cutoff(fnd, condID, colimit)
            % function co = cutoff(fnd, condID, NaNrate)
            %    -> return cutoff range (index for tstamp)
            %   co = 1 x nepoch cell
            %    if condID is specified, co{1} = ncond x 2 (start, end)
            %    if not specified, co{1} = 1 x 2 (start, end)
            %
            %  NaNrate: cut off criterion. If less than colimit*100 %
            %  trials exist, that timing will be cut off (default: 0.5)
            if ~exist('NaNrate', 'var')
                colimit = .5;
            end
            fnd = fnd.extract_valid_unit;
            
            nepoch = length(fnd.data);
            co = cell(1, nepoch);
            
            if ~exist('condID', 'var') || isempty(condID) % condition independent
                for n=1:nepoch
                    tot_tri = squeeze(sum(any(fnd.data{n}~=255,2),3)); % cell x 1
                    ind = squeeze(nanmean(bsxfun(@rdivide, sum(fnd.data{n}~=255, 3), tot_tri),1));
                    if all(ind < colimit)
                        error('No trial left with this cut off criterion');
                    end
                    co{n} = [find(ind >= colimit, 1), find(ind >= colimit, 1, 'last')];
                end
            else
                ncond = nanmax(condID(:));
                for n=1:nepoch
                    co{n} = nan(ncond, 2);
                    for c = 1:ncond
                        ind = nan(size(fnd.data{n},1), size(fnd.data{n},2));
                        for u=1:size(fnd.data{n},1)
                            d = fnd.data{n}(u,:, condID(u,:)==c);
                            tot_tri = squeeze(sum(any(d ~= 255, 2),3));
                            ind(u,:) = squeeze(sum(d ~= 255, 3)/tot_tri);
                        end
                        ind_mn = squeeze(nanmean(ind,1));
                        if all(ind_mn < colimit)
                            error('No trial left with this cut off criterion: cond %d', c);
                        end
                        if all(isnan(ind_mn))
                            continue;
                        end
                        co{n}(c,:) = [find(ind_mn >= colimit, 1), find(ind_mn >= colimit, 1, 'last')];
                    end
                end
            end
        end
        
        %% set_detrend_PSTH
        function out = set_detrend_PSTH(fnd, OVERWRITE)
            % set detrending PSTH to PSTH_trend field.
            % typically, PSTH function will just call it if detrend option
            % is on. But you want to call this function before that, if you
            % extract part of trials during your analysis (e.g. for cross
            % validation) but you want to detrend PSTH using all trials.
            % Then you call it before extracting trials (e.g. extract_trial
            % function).
            ncell = size(fnd.data{1},1);
            nepoch = length(fnd.data);
            if ~isempty(fnd.PSTH_trend) && (~exist('OVERWRITE', 'var') || ~OVERWRITE)
                fprintf('detrend PSTH is already set.\n');
                if nargout > 0
                    out = fnd.PSTH_trend;
                end
                return;
            end
            fnd.PSTH_trend = cell(1, nepoch);
            
            wb = waitbar_text(0);
            for n=1:nepoch
                fnd.PSTH_trend{n} = nan(ncell, size(fnd.data{n},2));
                for u = 1:ncell
                    waitbar_text((u + (n-1)*ncell)/(nepoch * ncell), wb);
                    raster = squeeze(double(fnd.data{n}(u,:,:)));
                    raster(raster==255) = NaN;
                    fnd.PSTH_trend{n}(u,:) = nanmean(raster,2)' * 1e3;
                end
            end
            waitbar_text('close', wb);
            if nargout > 0
                out = fnd.PSTH_trend;
            end
        end
        
        %% FR
        function r = FR(fnd, t_range, condID, zscoring, rmv_redundant_trial, ignore_nan)
            % r = FR(fnd, t_range, condID, zscoring)
            %   -> return mean firing rate
            % r = cell x trial double (spk/s)
            %
            % t_range: {epochID, [startT, endT]}
            % condID can be generated using trial_classifier.m
            % zscoring: if zscoring firing rate or not (default: false)
            fnd = fnd.extract_valid_unit;
            
            ncell = size(fnd.data{1},1);
            if ~iscell(t_range)
                if length(fnd.data) == 1
                    t_range = {1, t_range};
                else
                    error('t_range should specify event epoch');
                end
            end
            ev = t_range{1};
            t = t_range{2};
            
            tind = fnd.tstamp{ev} >= t(1) & fnd.tstamp{ev} < t(2);
            
                % convert to firing rate
            raster = double(fnd.data{ev}(:, tind, :));
            raster(raster==255) = NaN;
            
            if exist('ignore_nan', 'var') && ignore_nan
                r = squeeze(nanmean(raster,2)) * 1e3; % cell x trial
            else
                r = squeeze(mean(raster,2)) * 1e3; % cell x trial
                    % if there is any missing value (NaN) in this time window,
                    % FR for that trial will be NaN
            end            
            if exist('zscoring', 'var') && zscoring
                for n = 1:ncell
                    r(n, ~isnan(r(n,:))) = zscore(r(n, ~isnan(r(n,:))));
                end
            end
            
            if exist('condID', 'var') && ~isempty(condID)
                ncond = nanmax(condID(:));
                rn = nan(ncell, ncond);
                for n = 1:ncell
                    for c = 1:ncond
                        rn(n,c) = nanmean(r(n, condID(n,:) == c),2);
                    end
                end
                r = rn;
            elseif ~exist('rmv_redundant_trial', 'var') || rmv_redundant_trial
                if size(r,2) > max(fnd.ntrial)
                    r = r(:, 1:max(fnd.ntrial));
                end
            end
        end
        
        %% set_data
        function set_data(fnd, data)
            fnd.compress_data(data);
        end
                
        %% get_data_info
        function i = get_data_info(fnd, command, raw)
            % function i = info(fnd, command)
            %   -> return fnd data info (valid_unit only)
            %      if command is not specified, return a struct including
            %      all info
            %  command: nepoch, nunit, unitID, unit_quality, unit_type
            i.nallunit = size(fnd.unitID,1);
            fi = fieldnames(fnd.param);
            i.nalltrial = size(fnd.getp(fi{1}),2);
            
            if ~exist('raw', 'var') || ~raw
                % if raw, include invalid units
                fnd = fnd.extract_valid_unit;
            end
            i.nepoch = length(fnd.data);
            i.nunit = size(fnd.unitID,1);
            i.ntrial = size(fnd.getp(fi{1}),2);
            i.unitID = fnd.unitID;
            i.unit_quality = fnd.unit_quality;
            i.unit_type = fnd.unit_type;
            i.monkey = fnd.monkey;
            i.area = fnd.area;
            if exist('command', 'var') && ~isempty(command)
                i = i.(command);
            end
        end
        
        
        %% disp
        function disp(fnd)
            % default display option of FND function
            fprintf('FND data class\n');
            fprintf('----BASIC INFO-------\n');
            fprintf([fnd.description '\n']);
            
            % info
            fprintf('\tMonkey:');
            M = fnd.monkey(~cellfun(@isempty, fnd.monkey));
            for s = unique(M)'
                fprintf('%s, ', s{1});
            end
            fprintf('\n');
            fprintf('\tArea:');
            A = fnd.area(~cellfun(@isempty, fnd.area));
            if ~isempty(A)
                for s = unique(A)'
                    fprintf('%s, ', s{1});
                end
            end
            fprintf('\n');
            
            % unit
            fprintf(['\tUnit criteria: ' fnd.unit_criteria_description, '\n']);
            if all(fnd.valid_unit)
                fprintf('\t%d sessions\n', length(unique(fnd.unitID(:,1))));
                fprintf('\t%d units\n', length(fnd.valid_unit));
            else
                fprintf('\t%d/%d valid sessions\n', length(unique(fnd.unitID(fnd.valid_unit,1))), ...
                    length(unique(fnd.unitID(:,1))));
                fprintf('\t%d/%d valid units\n', sum(fnd.valid_unit), length(fnd.valid_unit));
            end
            
            % alignment
            for a = 1:length(fnd.alignto)
                fprintf('\tAlign %d: %s, from %d to %d', a, fnd.alignto(a).event, fnd.alignto(a).start_offset, fnd.alignto(a).end_offset);
                fprintf(' (raster: %d units x %d ms x %d trials)\n', size(fnd.data{a}));
            end
            
            % parameter
            if isstruct(fnd.param)
                fprintf('\tParameter:');
                fi = fieldnames(fnd.param);
                for n=1:length(fi)
                    fprintf('%s, ', fi{n});
                end
                fprintf('\n');
            end
            
            % other info
            fprintf('\tSize:');
            fnd.getsize();
            
            fprintf('\tCreated: %s\n', fnd.date);
            fprintf('\tClass version: %s\n', fnd.version);
            
            fprintf('\n-----USAGE--------\n');
            fnd.show_help();
            
            fprintf('\n\n');
        end
        
    end
    %% @@@@ public static methods @@@@
    methods (Static)
        %% mergedata
        function fnd = mergedata(fnds)
            fnd = copy(fnds{1});
            mergefield = {'unit_type', 'session_name', 'unit_quality', 'monkey', 'area', 'valid_unit', 'ntrial'};
            pfield = fieldnames(fnd.param);
            mfield = fieldnames(fnd.misc);
            for n=2:length(fnds)
                if ~isequal(fnd.alignto, fnds{n}.alignto)
                    error('alignto does not match between FNDs');
                end
                uID = fnds{n}.unitID;
                uID(:,1) = uID(:,1) + length(fnd.session_name);
                fnd.unitID = [fnd.unitID; uID];
                for d=1:length(fnd.data)
                    fnd.data{d} = catdata(fnd.data{d}, fnds{n}.data{d});
                end
                for f=1:length(mergefield)
                    fnd.(mergefield{f}) = [fnd.(mergefield{f}); fnds{n}.(mergefield{f})];
                end
                for p=1:length(pfield)
                    fnd.param.(pfield{p}) = catdata(fnd.param.(pfield{p}), fnds{n}.param.(pfield{p}));
                end
                for p=1:length(mfield)
                    fnd.misc.(mfield{p}) = [fnd.misc.(mfield{p}); fnds{n}.misc.(mfield{p})];
                end
                fnd.PSTH_trend = [fnd.PSTH_trend; fnds{n}.PSTH_trend];
            end
            
            function dn = catdata(d1,d2)
                s1 = size(d1);
                s2 = size(d2);
                m = length(s1);
                if s1(m) > s2(m)
                    nd = s2;
                    nd(m) = s1(m) - s2(m);
                    if iscell(d2)
                        d2 = cat(m, d2, cell(nd));
                    elseif isinteger(d2)
                        d2 = cat(m, d2, uint8(ones(nd)*255));
                    else
                        d2 = cat(m, d2, nan(nd));
                    end
                elseif s1(m) < s2(m)
                    nd = s1;
                    nd(m) = s2(m) - s1(m);
                    if iscell(d1)
                        d1 = cat(m, d1, cell(nd));
                    elseif isinteger(d1)
                        d1 = cat(m, d1, uint8(ones(nd)*255));
                    else
                        d1 = cat(m, d1, nan(nd));
                    end
                end
                dn = cat(1, d1, d2);
            end
        end
    end
    
    %% @@@@ internal methods @@@@
    methods (Access = protected)
        %% compress_data
        function compress_data(fnd, data)
            compress_code = uint8(127);
            if nargin < 2
                data = fnd.data;
            end
            if isa(data{1}, 'double')
                for a=1:length(data)
                    data{a}(isnan(data{a})) = 255;
                    data{a} = uint8(data{a});
                end
            end
            if data{1}(1) ~= compress_code && isequal(fnd.compression_type, 'zip')
                error('not implemented well yet');
                try
                    for a=1:length(data)
                        s = size(data{a});
                        M=[compress_code; uint8(length(s));typecast(s(:),'uint8');data{a}(:)];
                        f=java.io.ByteArrayOutputStream();
                        g=java.util.zip.DeflaterOutputStream(f);
                        g.write(M);
                        g.close;
                        data{a} = typecast(f.toByteArray,'uint8');
                        f.close;
                    end
                catch ME
                    if ~isempty(strfind(ME.message,'java.lang.OutOfMemoryError'))
                        URL = 'https://www.mathworks.com/help/matlab/matlab_external/java-heap-memory-preferences.html';
                        error(['Java heap space out of memory\n' ...
                            'see %s'], URL);
                    end
                    rethrow(ME);
                end
            end
            fnd.data = data;
        end
        
        %% extract_valid_unit
        function nfnd = extract_valid_unit(fnd)
            nfnd = copy(fnd);
            v = ~nfnd.valid_unit;
            nfnd.unitID(v,:) = [];
            % [sesID, ~, nfnd.unitID(:,1)] = unique(nfnd.unitID(:,1));
            % do not change session ID because it could be confusing if we
            % want to compare multiple datasets
            % nfnd.session_name = nfnd.session_name(sesID);
            nfnd.unit_type(v) = [];
            nfnd.unit_quality(v) = [];
            nfnd.valid_unit(v) = [];
            nfnd.monkey(v) = [];
            nfnd.area(v) = [];
            nfnd.ntrial(v) = [];
            fi = fieldnames(nfnd.param);
            for n=1:length(fi)
                nfnd.param.(fi{n})(v,:) = [];
            end
            if ~isempty(nfnd.misc)
                fi = fieldnames(nfnd.misc);
                for n=1:length(fi)
                    if size(nfnd.misc.(fi{n}),1) == length(v)
                        nfnd.misc.(fi{n})(v,:) = [];
                    end
                end
            end
            for n=1:length(nfnd.data)
                if ~isempty(nfnd.PSTH_trend)
                    nfnd.PSTH_trend{n}(v,:) = [];
                end
                nfnd.data{n}(v,:,:) = [];
            end
            nfnd.date = datestr(now);
        end
        
        
        %% getsize
        function getsize(fnd)
           props = properties(fnd); 
           totSize = 0;

           for ii=1:length(props)
              currentProperty = fnd.(props{ii}); %#ok<NASGU>
              s = whos('currentProperty'); 
              totSize = totSize + s.bytes; 
           end
           fprintf('%1.2f MB\n', totSize/1e3/1e3);
        end
    end
    

    %% @@@@ internal methods (static) @@@@
    methods (Static, Access = protected)
        %% nanconv_data
        function r = nanconv_data(r, conv)
            % conv ... {type, value}
            % type ... boxcar, gaussian
            % value ... boxcar width, gaussian SD
            if iscell(conv)
                switch conv{1}
                    case 'boxcar'
                        kernel = fspecial('average', [1, conv{2}]);
                    case 'gaussian'
                        kernel = fspecial('gaussian', [1 conv{2}*6], conv{2});
                    case 'bartlett'
                        kernel = 1 - abs(linspace(-1, 1, conv{2})); % Symmetric triangular window
                        kernel = kernel / sum(kernel); % Normalize for unity gain
                    otherwise
                        error('No such filter: %s', conv{1});
                end
            else
                kernel = conv;
            end
            r = nanconv(r, kernel, 'same');
        end
        
        %% show_help
        function show_help()
            fprintf(' r = fnd.raster():  return unit x time x trial raster\n');
            fprintf(' r = fnd.raster_cond(condID): return 1 x cond cell (each contains raster)\n');
            fprintf(' r = fnd.PSTH(condID): return unit x time x cond PSTH\n');
            fprintf(' co = fnd.cutoff(condID, colimit): return cut off range of PSTH\n');
            fprintf(' r = fnd.FR(time): return unit x trial mean firing rate\n');
            fprintf(' fnd = fnd.set_unit_criteria(varargin)\n');
            fprintf(' out = fnd.getp(param_name): return task parameters\n');
            fprintf(' fnd.setp(param_name, param1): set task parameters\n');
            fprintf(' fnd = fnd.extract_trial(trial_flg): extract a part of trials\n');
            fprintf(' fnd.get_data_info(param): return data info\n');
            fprintf(' fnd.set_detrend_PSTH: set detrending PSTH\n');
        end
        
        %% check_unit_info
        function info = check_unit_info(info, nunit, len, name)
            if isempty(info)
                info = nan(nunit, len);
            elseif size(info,1) ~= nunit
                error('%s and spike data do not match', name);
            end
        end
    end
end

