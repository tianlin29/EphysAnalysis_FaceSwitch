function [fh_summary, fh_all_psych] = plot_psychometric_function(coh, cor, cond, session, repeat, opt_classifier, experiment)

n_timebins = length(opt_classifier.tstamp{1});
n_repeats = max(repeat); % repetition is not considered here
n_files = max(session);
n_conds = max(cond);

[fh_summary, fh_all_psych] = deal(cell(n_timebins,1));
for t = 1:n_timebins
    clear opt
    % select data
    opt.session_list = 1:n_files;
    % process data
    opt.log = true;
    opt.constant = false;
    opt.verbose = false;
    % plot
    switch experiment
        case {'learnTask2', 'threeExemplar'}
            opt.color = [0 0 0; 44 145 224; 255 0 0]/255;
        case 'learnTask3'
            opt.color = [0 0 0; 58 191 153; 255 0 0]/255;
        case 'learnTask4'
            opt.color = [0 0 0; 240 169 58; 255 0 0]/255;
        case 'faceColor'
            opt.color = [0 0 0; 1 0 0; 1 0 0];
    end
    opt.linewidth = [0.5, 0.5, 0.5];
    opt.average = false; % do not average across learning sessions
    opt.normalize_threshold = 0;

    [~, fh_idv, fh_summary{t}, stat] = run_unsigned_choice_3cond(cond, coh, cor(:,t), session, opt);
    fh_all_psych{t} = plot_in_one_fig_unsigned_choice_2cond(fh_idv, [5 5], [500 500]*1.5);
end

end