
function [hf, ax] = plotPSTHsingle(fig_title, raster, unit_id, alignto, options)
%
% function plotPSTHsingle(fig_title, raster, unit_id, alignto, options)
%
% shows PSTHs of many simultaneously recorded units in a compact form for quick navigation
%
% Inputs:
%   fig_title   title of the figure, usually the name of the experiment or FIRA file
%   raster      a cell array or rasters aligned to different events. raster{a} must be the
%               output of getRasterArray and be of the following form raster{a}=[units x trials x time]
%   unit_id     unit IDs returned by getRasterArray, it must be a matrix with rows
%               corresponding to units. each row must be [channel_num unit_num].
%   alignto     a struct array. these are the alignto structures passed to getRatserArray
%   options
%       color_lib   a cell array that defines the color of different trial groups
%       style_lib   a cell array that defines the line style of different trial groups
%       legend      a cell array that contains labels of different trial groups
%       show_se     show the standard error of the PSTHs. default is false because calculation of
%                   standard errors is time consuming. explore your data with show_se=false
%                   but for final production set it to true.
%       zoom_mode   defines what will be shown when a panel is clicked. if 'all_alignments', all panels
%                   of a unit are shown. if 'single_alignment', only the clicked panel will be zoomed in.
%
% Outputs:
%   hf          figure handle
%   ax          axis handles
%
% Example:
%   Draw array PSTHs aligned to target onset, dots onset, saccade, and feedback for the
%   right and left choices in the direction discrimination task. The code assumes that
%   FIRA is already loaded and basic variables like valid_trials, cor, trg_cho, and
%   trg_right are already defined.
%
%             %trial groups
%         trials_cor{1} = find(valid_trials & (cor==1|coh==0) & trg_cho==trg_right);   %correct right choices
%         trials_cor{2} = find(valid_trials & (cor==1|coh==0) & trg_cho~=trg_right);   %correct left choices
%
%             %alignto structure
%         align_struct = {    %align to target onset
%                         {'event', 'targ_on'; 'start_offset', -300; 'end_offset', 500; 'limits', {{'in_fp_wd',0}, {'dots_on_screen',0}}}
%                             %align to dots onset
%                         {'event', 'dots_on_screen'; 'start_offset', -300; 'end_offset', 1000; 'limits', {{'targ_on',0}, {'out_fp_wd',-50}}}
%                             %align to saccade
%                         {'event', 'out_fp_wd'; 'start_offset', -1000; 'end_offset', 300; 'limits', {{'dots_on_screen',200}, {'feedback',0}}}
%                             %align to feedback
%                         {'event', 'feedback'; 'start_offset', -250; 'end_offset', 500; 'limits', {{'out_fp_wd',100}, {'feedback',500}}}
%         };
%             %recast align_struct into a structure
%         for i = 1 : length(align_struct)
%             alignto(i) = cell2struct(align_struct{i}(:,2), align_struct{i}(:,1), 1);
%         end
%
%             %get the raster of individual trials
%         raster = cell(1, length(alignto));
%         hw = waitbar_text(0, [], 'getting rasters for all alignments');
%         for a = 1 : length(alignto)
%             [raster{a}, unit_id] = getRasterArray(FIRA, alignto(a), trials_cor, alignto(a).limits{1}, alignto(a).limits{2});
%             waitbar_text(a/length(alignto), hw);
%         end
%         waitbar_text('close', hw);
%
%             %show the PSTHs
%         plotPSTHArray(file_name, raster, unit_id, alignto, {'r','b'},{'-','--'}, false);
%


%
% RK, 1/2/2017
%


cur_page = 1;
opt = struct('max_unit_per_page', 14, ...
    'max_align_per_unit', 5, ...
    'conv_kernel', fspecial('average', [1 50]), ...
    'color_lib', {repmat(distinctPalette',[1 5])}, ...
    'style_lib', {repmat({'-'},[1 80])}, ...
    'legend', [], ...
    'show_se', false, ...
    'zoom_mode', 'all_alignments');

if nargin>=5 && ~isempty(options)
    opt = safeStructAssign(opt, options);
    if mod(opt.max_unit_per_page,2)~=0
        opt.max_unit_per_page = opt.max_unit_per_page + 1;
    end
end

%unit number and unit names
unit_num = size(unit_id, 1);
unit_name = cell(unit_num, 1);
for u = 1 : unit_num
    unit_name{u} = sprintf('u%02d-%d', unit_id(u,1), unit_id(u,2));
end

%smooth the rasters
if opt.show_se
    hw = waitbar_text(0, [], 'smoothing rasters');
    for a = 1 : length(alignto)
        raster{a} = smoothRasterArray(raster{a}, opt.conv_kernel);
        waitbar_text(a/length(alignto), hw);
    end
    waitbar_text('close', hw);
end

%calculate the PSTHs
psth = cell(size(alignto));     %psth{a} = [units x time x condition]
psth_se = cell(size(psth));
curve_cutoffs = cell(size(alignto));    %curve_cutoffs{a} = [low high; low high; ...]
axis_cutoffs = nan(length(alignto),2);
t = cell(size(alignto));
condition_num = nan(size(alignto));
hw = waitbar_text(0, [], 'calculating PSTHs');
for a = 1 : length(alignto)
    if a>opt.max_align_per_unit
        break;
    end
    %number of conditions and time line
    condition_num(a) = length(raster{a});
    t{a} = alignto(a).start_offset:alignto(a).end_offset;

    psth{a} = nan(size(raster{a}{1},1), size(raster{a}{1},3), condition_num(a));   %[units x time x condition]
    psth_se{a} = nan(size(psth{a}));
    curve_cutoffs{a} = nan(condition_num(a), 2);
    for i = 1 : condition_num(a)
        %mean response of each unit, averaged across trials
        psth{a}(:,:,i) = squeeze(nanmean(raster{a}{i},2)) * 1e3;    %convert to sp/s
        n = squeeze(sum(~isnan(raster{a}{i}),2));
        if opt.show_se==true
            %standard error of each unit across the trials
            sd = squeeze(nanstd(raster{a}{i},[],2));
            psth_se{a}(:,:,i) = sd./sqrt(n) * 1e3;      %convert to sp/s
        else
            for u = 1 : unit_num
                psth{a}(u,:,i) = nanconv(squeeze(psth{a}(u,:,i)), opt.conv_kernel, 'same');
            end
        end
        %for plotting the PSTHs we limit the curves to periods when at least 50% of
        %trials contributed to the average
        trial_num = size(raster{a}{i},2);
        nav = mean(n,1);
        curve_cutoffs{a}(i,:) = [find(nav>=max(nav)/2,1,'first'), find(nav>=max(nav)/2,1,'last')];
        % curve_cutoffs{a}(i,:) = [find(n(1,:)>=trial_num/2,1,'first'), find(n(1,:)>=trial_num/2,1,'last')];
    end

    axis_cutoffs(a,:) = [floor(min(t{a}(curve_cutoffs{a}(:,1)))/50)*50  %earliest cutoff rounded to 50ms
        ceil(max(t{a}(curve_cutoffs{a}(:,2)))/50)*50]; %latest cutoff rounded to 50ms
    waitbar_text(a/length(alignto), hw);
end
waitbar_text('close', hw);


% Create the GUIs but hide them until all PSTHs are drawn
pos = [0 50 1280 750];
screen_rect = get(0, 'ScreenSize');
pos(3) = pos(3)+0.5*max(0,screen_rect(3)-pos(3));
pos(4) = pos(4)+0.5*max(0,screen_rect(3)-50-pos(3));
hf = figure('Visible', 'off', 'Color', 'w', ...
    'Position', pos, 'PaperPositionMode', 'auto', ...
    'PaperUnits', 'points', 'PaperSize', pos(3:4), ...
    'CloseRequestFcn', @closeFigCallback, ...
    'MenuBar', 'none', ...
    'ToolBar', 'figure');
ax = nan(opt.max_unit_per_page, opt.max_align_per_unit);
hnext = uicontrol(hf, 'Style', 'pushbutton', 'String', '>', 'Position', [pos(3)/2+5 10 30 15], 'Callback', @nextCallback);
hprev = uicontrol(hf, 'Style', 'pushbutton', 'String', '<', 'Position', [pos(3)/2-35 10 30 15], 'CallBack', @prevCallback);
hnext_mult = uicontrol(hf, 'Style', 'pushbutton', 'String', '>>', 'Position', [pos(3)/2+40 10 30 15], 'Callback', @nextMultCallback);
hprev_mult = uicontrol(hf, 'Style', 'pushbutton', 'String', '<<', 'Position', [pos(3)/2-70 10 30 15], 'CallBack', @prevMultCallback);
hmenu(1) = uimenu(hf, 'Label', 'Save');
uimenu(hmenu(1), 'Label', 'Single units PSTHs', 'CallBack', @printAllUnitFigs);
hmenu(2) = uimenu(hf, 'Label', 'GoTo');
uimenu(hmenu(2), 'Label', 'Page', 'CallBack', @gotoPageCallback);
drawPage(cur_page);
% change units to normalized so components resize automatically.
set([hf,hnext,hprev,hnext_mult,hprev_mult], 'Units', 'Normalized');
% now make the GUI visible.
set(hf, 'Visible', 'on');



%% drawPage
    function drawPage(page)
        ax_h = 0.92/(opt.max_unit_per_page/2);
        ax_w = diff(axis_cutoffs,1,2);
        ax_w = 0.35*ax_w/nansum(ax_w);
        if length(alignto)==1
            ax_x_gap = 0.05;
        else
            ax_x_gap = 0.05/(length(alignto)-1);
        end

        figure(hf);
        set(hf, 'Name', sprintf('%s, Page %d, units %d-%d',fig_title,page,(page-1)*opt.max_unit_per_page+1,min(page*opt.max_unit_per_page,unit_num)));

        %if the figure is opened from a .fig file, ax need to be refreshed
        if any(~ishandle(ax(:)))
            hf = gcf;
            getAxesFromFigure;
        end

        %clear the figure
        for jj = 1 : numel(ax)
            if ishandle(ax(jj)) && ax(jj)~=0
                cla(ax(jj),'reset');
                set(ax(jj), 'Visible', 'off');
                % delete(ax(jj));
            end
        end

        %make the requested page
        hleg = nan(condition_num(1), 1);
        for uu = (page-1)*opt.max_unit_per_page+1 : page*opt.max_unit_per_page
            if uu>unit_num
                break;
            end
            %where should the PSTHs of this unit be shown on the figure?
            row = ceil((uu-(page-1)*opt.max_unit_per_page)/2);
            col = 1-mod((uu-(page-1)*opt.max_unit_per_page),2);
            ax_x = col*0.5 + 0.05;
            ax_y = 0.98-ax_h*row;

            %make the PSTHs aligned to different events for this unit. make sure the
            %time axis is not unequally stretched for different alignments (equalize the spacing of x-axis ticks)
            ylim = nan(length(alignto), 2);
            for aa = 1 : length(alignto)
                if aa>opt.max_align_per_unit
                    break;
                end

                ax(uu,aa) = subplot('position', [ax_x+ax_x_gap*(aa-1)+sum(ax_w(1:aa-1)) ax_y ax_w(aa) ax_h*0.9]);
                set(ax(uu,aa), 'ButtonDownFcn', @axisCallback, 'UserData', [uu,aa]);

                set(gca, 'Visible', 'on');
                hold on;
                for ii = 1 : condition_num(aa)
                    h = plot(t{aa}(curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2)), ...
                        psth{aa}(uu,curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2),ii), ...
                        'Color', opt.color_lib{ii}, 'LineStyle', opt.style_lib{ii});
                    set(h, 'ButtonDownFcn', {@plotRasters,uu,aa,ii});
                    if uu==(page-1)*opt.max_unit_per_page+1 && aa==1
                        hleg(ii) = h;
                    end
                    if opt.show_se==true
                        h = patch([t{aa}(curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2)), fliplr(t{aa}(curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2)))], ...
                            [psth{aa}(uu,curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2),ii)+psth_se{aa}(uu,curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2),ii), ...
                            fliplr(psth{aa}(uu,curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2),ii)-psth_se{aa}(uu,curve_cutoffs{aa}(ii,1):curve_cutoffs{aa}(ii,2),ii))], ...
                            opt.color_lib{ii}, 'EdgeColor', opt.color_lib{ii}, 'EdgeAlpha', 0.5, 'FaceColor', opt.color_lib{ii}, 'FaceAlpha', 0.5);
                        set(h, 'ButtonDownFcn', {@plotRasters,uu,aa,ii});
                    end
                end
                xtick = axis_cutoffs(aa,1):100:axis_cutoffs(aa,2);
                set(gca, 'XLim', axis_cutoffs(aa,:), 'XTick', xtick, 'TickDir', 'out');
                ylim(aa,:) = get(gca, 'YLim');
                if row==opt.max_unit_per_page/2 || uu==unit_num || uu==unit_num-1
                    set(gca, 'XTickLabel', makeTickLabel(xtick,200));
                    if aa==1
                        xlabel('Time (ms)');
                    end
                else
                    set(gca, 'XTickLabel', []);
                end
                if (col==0 && aa==length(alignto)) || (col==1 && aa==1)
                    set(gca, 'YTickLabel', []);
                end
                if (row==1 || row==opt.max_unit_per_page/2) && ...
                        ((col==0 && aa==1) || (col==1 && aa==length(alignto)))
                    ylabel('Firing rate (sp/s)');
                end
            end

            %equalize the ylim of the PSTHs aligned to different events
            ylim = [min(ylim(:,1)) max(ylim(:,2))];
            ylim = ylim + diff(ylim)*0.03*[-1 +1];
            for aa = 1 : length(alignto)
                plot(ax(uu,aa), [0 0], ylim, '--k');
                set(ax(uu,aa), 'YLim', ylim);
                if aa==1
                    text(axis_cutoffs(aa,1)+0.1*diff(axis_cutoffs(aa,:)), ylim(1)+0.8*diff(ylim), unit_name{uu}, 'Parent', ax(uu,aa), 'FontSize', 9, 'Color', 'k', 'FontWeight', 'Bold');
                elseif aa>1 && aa<length(alignto)
                    set(ax(uu,aa), 'YColor', 'w', 'YTickLabel', []);
                elseif aa==length(alignto)
                    set(ax(uu,aa), 'YAxisLocation', 'right');
                end
                if row==1
                    text(0, ylim(2)*1.02, strrep(alignto(aa).event,'_',' '), 'Parent', ax(uu,aa), 'FontSize', 9, 'HorizontalAlignment', 'Center');
                end
            end
        end
        if ~isempty(opt.legend)
            if ishandle(ax(1,1))
                if length(opt.legend)==length(hleg)
                    [h, hobj] = legend(ax(1,1), hleg, opt.legend);
                    set(h, 'FontSize', 8);
                    legendShrink(h, hobj, 0.6);
                    pos = get(h, 'Position');
                    set(h, 'Position', [0.47 pos(2:4)]);
                    legend(ax(1,1), 'boxoff');
                    % uistack(ax(1,1), 'Layer', 'Top');
                else
                    warning('plotPSTHArray:LegendMismatch', 'legend ignored because text does not match the number of conditions');
                end
            end
        end
        drawnow;

        %update page counter
        cur_page = page;
    end



%% getAxesFromFigure
% returns the list of axes in the figure currently open. handy when the figures is opened
% by loading a .fig file.
    function getAxesFromFigure
        hchild = get(hf,'Children');
        haxes = hchild(strcmp(get(hchild,'Type'),'axes'));
        haxes = reshape(haxes, [min(length(alignto),opt.max_align_per_unit) length(haxes)/min(length(alignto),opt.max_align_per_unit)])';
        haxes = rot90(haxes,2);
        ax(1:size(haxes,1),1:size(haxes,2)) = haxes;
    end



%% nextCallback
%callback function for next-page button
    function nextCallback(~, ~)
        if cur_page<ceil(unit_num/opt.max_unit_per_page)
            drawPage(cur_page+1);
        end
    end



%% prevCallback
%callback function for previous-page button
    function prevCallback(~, ~)
        if cur_page>1
            drawPage(cur_page-1);
        end
    end



%% nextCallback
%callback function for next-page button
    function nextMultCallback(~, ~)
        if cur_page<ceil(unit_num/opt.max_unit_per_page)
            drawPage(min(cur_page+5,ceil(unit_num/opt.max_unit_per_page)));
        end
    end



%% prevCallback
%callback function for previous-page button
    function prevMultCallback(~, ~)
        if cur_page>1
            drawPage(max(cur_page-5,1));
        end
    end



%% axisCallback
%draw the PSTHs in a bigger window
    function axisCallback(source, ~)
        axis_data = get(source,'UserData');
        uu = axis_data(1);  %which unit
        aa = axis_data(2);  %which alignment
        %if the figure is opened from a .fig file, ax need to be refreshed
        if any(~ishandle(ax(:)))
            getAxesFromFigure;
        end
        %plot larger PSTHs
        if ~isfield(opt,'zoom_mode') || strcmpi(opt.zoom_mode,'all_alignments')
            %if all alignments are copied (default function) use the cloneAllAlignments function
            hf_unit = figure('Name', unit_name{uu}, 'Color', 'w', 'Position', [100 100 640 280], 'PaperUnits', 'points', 'PaperSize', [640 280], 'PaperPositionMode', 'auto');
            cloneAllAlignments(uu, hf_unit);
        elseif strcmpi(opt.zoom_mode,'single_alignment')
            %if a single alignments is zoomed into, draw only the clicked alignment
            figure('Name', sprintf('%s, %s',unit_name{uu},alignto(aa).event), ...
                'Color', 'w', 'Position', [100 100 500 375], 'PaperUnits', 'points', ...
                'PaperSize', [500 375], 'PaperPositionMode', 'auto');
            ax2 = subplot(1,1,1);
            copyobj(allchild(source), ax2);
            %adjust axes
            set(ax2, 'XLim', get(ax(uu,aa),'XLim'), 'XTick', get(ax(uu,aa),'XTick'), ...
                'YLim', get(ax(uu,aa),'YLim'), 'YTick', get(ax(uu,aa),'YTick'), 'TickDir', 'out');
            ylabel('Firing rate (sp/s)');
            xlabel('Time (ms)');
            %restore callback functions
            h1 = findall(ax(uu,aa), 'Type', 'Line');
            h2 = findall(ax2, 'Type', 'Line');
            for ii = 1 : length(h1)
                set(h2(ii), 'ButtonDownFcn', get(h1(ii),'ButtonDownFcn'));
            end
            if opt.show_se==true
                h1 = findall(ax(uu,aa), 'Type', 'Patch');
                h2 = findall(ax2, 'Type', 'Patch');
                for ii = 1 : length(h1)
                    set(h2(ii), 'ButtonDownFcn', get(h1(ii),'ButtonDownFcn'));
                end
            end
        end
    end



%% plotRasters
%callback function that creates the corresponding raster plots when a psth is double-clicked
    function plotRasters(~, ~, uu, aa, ii)
        %if the figure is opened from a .fig file, ax need to be refreshed
        if any(~ishandle(ax(:)))
            getAxesFromFigure;
        end
        %now draw rasters
        figure('Name', sprintf('%s, %s, PSTH %d',unit_name{uu},alignto(aa).event,ii), ...
            'Color', 'w', 'Position' , [750 100 400 300], 'PaperUnits', 'points', 'PaperSize', [400 300], 'PaperPositionMode', 'auto');
        hold on;
        r = squeeze(raster{aa}{ii}(uu,:,:));
        if opt.show_se==true
            r = r * 1e3;
        else
            r(r==0) = 0.9;
            r(r==1) = 0;
            r(isnan(r)) = 1;
        end
        trial_num = size(raster{aa}{ii}(uu,:,:),2);
        plot([0 0], [0 trial_num+1], 'r-', 'LineWidth', 2);
        h = imagesc(t{aa}, 1:trial_num, r);
        alpha(h, 0.5);
        text(0, trial_num*1.07, strrep(alignto(aa).event,'_',' '), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
        max_trial_per_cond = max(cellfun(@(x)size(x,2),raster{aa}));
        set(gca, 'XLim', get(ax(uu,aa),'XLim'), 'YLim', [0 max_trial_per_cond+1], 'TickDir', 'out');
        % axis image
        xlabel('Time (ms)');
        ylabel('Trial #' );
        if opt.show_se==true
            colormap(inferno);   %plasma, magma, inferno, viridis
            h = colorbar;
            hpos = get(h, 'Position');
            apos = get(gca, 'Position');
            set(h, 'Position', hpos.*[1 1 0.5 1], 'TickDir', 'out');
            set(gca, 'Position', apos);
            ylabel(h, 'Firing rate (sp/s)');
        else
            %             colormap(hot);
            colormap(gray);
        end
    end



%% closeFigCallback
    function closeFigCallback(~, ~)
        selection = questdlg('Close Array PSTH Viewer?', 'Close', 'Yes', 'No', 'No');
        if strcmp(selection,'Yes')
            %close the main dialog box
            delete(hf);
        end
    end



%% cloneAllAlignments
% plots the PSTHs of unit uu into figure handle hf_unit
    function cloneAllAlignments(uu, hf_unit)
        figure(hf_unit);
        ax2 = nan(1, length(alignto));
        ax_h = 0.85;
        ax_w = diff(axis_cutoffs,1,2);
        ax_w = 0.78*ax_w/nansum(ax_w);
        ax_x_gap = 0.1/(length(alignto)-1);
        row = ceil((uu-(cur_page-1)*opt.max_unit_per_page)/2);
        for aa = 1 : length(alignto)
            if aa>opt.max_align_per_unit
                break;
            end
            ax2(aa) = subplot('position', [0.07+ax_x_gap*(aa-1)+sum(ax_w(1:aa-1)) 0.1 ax_w(aa) ax_h]);
            copyobj(allchild(ax(uu,aa)), ax2(aa));
            %adjust axes
            xtick = get(ax(uu,aa),'XTick');
            set(ax2(aa), 'XLim', get(ax(uu,aa),'XLim'), 'XTick', xtick, 'XTickLabel', makeTickLabel(xtick,200), ...
                'YLim', get(ax(uu,aa),'YLim'), ... %'YTick', get(ax(uu,aa),'YTick'), ...
                'YColor', get(ax(uu,aa),'YColor'), 'YAxisLocation', get(ax(uu,aa),'YAxisLocation'), 'TickDir', 'out');
            if aa==1
                ylabel('Firing rate (sp/s)');
                xlabel('Time (ms)');
            end
            %restore callback functions
            h1 = findall(ax(uu,aa), 'Type', 'Line');
            h2 = findall(ax2(aa), 'Type', 'Line');
            for ii = 1 : length(h1)
                set(h2(ii), 'ButtonDownFcn', get(h1(ii),'ButtonDownFcn'));
            end
            if opt.show_se==true
                h1 = findall(ax(uu,aa), 'Type', 'Patch');
                h2 = findall(ax2(aa), 'Type', 'Patch');
                for ii = 1 : length(h1)
                    set(h2(ii), 'ButtonDownFcn', get(h1(ii),'ButtonDownFcn'));
                end
            end
            if row~=1
                ylim = get(ax(uu,aa),'YLim');
                text(0, ylim(2)*1.02, strrep(alignto(aa).event,'_',' '), 'Parent', ax2(aa), 'FontSize', 9, 'HorizontalAlignment', 'Center');
            end
        end
    end



%% printAllUnitFigs
%save the PSTHs of all units, making one png file per unit
    function printAllUnitFigs(~, ~)
        %ask where to save the png files
        [base_name,save_path,filter_ind] = uiputfile({'*.png';'*.pdf';'*.eps';'*.tif'}, 'Save PSTHs of all units', fig_title);
        if ~ischar(base_name)
            return;
        end
        switch filter_ind
            case 1, driver='-dpng'; extension='.png';
            case 2, driver='-dpdf'; extension='.pdf';
            case 3, driver='-epsc'; extension='.eps';
            case 4, driver='-dtiff';    extension='.tif';
            otherwise,  error('Unrecognized image format');
        end
        [~,base_name] = fileparts(base_name);
        %save the current page in the population PSTH
        cur_page_bkp = cur_page;
        %define the order of pages, put the current page at the beginning to avoid
        %replotting the PSTHs that are already on the screen (this saves time)
        page_order = [cur_page setdiff(1:ceil(unit_num/opt.max_unit_per_page),cur_page)];
        %plot the units one by one, paging through all units
        hf_unit = figure('Color', 'w', 'Position', [100 100 640 280], 'PaperUnits', 'points', 'PaperSize', [640 280], 'PaperPositionMode', 'auto');
        for page = page_order
            if page~=cur_page
                drawPage(page);
            end
            for uu = (page-1)*opt.max_unit_per_page+1 : page*opt.max_unit_per_page
                if uu>unit_num
                    break;
                end
                clf(hf_unit);
                set(hf_unit, 'Name', unit_name{uu});
                cloneAllAlignments(uu, hf_unit);
                print(hf_unit, driver, '-r300', fullfile(save_path,[base_name,'_',unit_name{uu},extension]));
            end
        end
        close(hf_unit);
        %go back to the page that was viewed before the units were printed
        if cur_page_bkp~=cur_page
            drawPage(cur_page_bkp);
        end
    end



%% gotoPageCallback
    function gotoPageCallback(~, ~)
        answer = inputdlg('Go to page: ', '', 1, {'1'});
        if isempty(answer)
            return;
        end
        page = str2double(answer{1});
        if page<ceil(unit_num/opt.max_unit_per_page)
            drawPage(page);
        end
    end
end
