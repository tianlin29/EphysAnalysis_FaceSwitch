function [coh, cor, cond, session, repeat, opt_classifier] = merge_across_repetition(InterimDir, monkey, experiment, n_repeats, n_files, n_conds)

opt_classifier = load(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat1.mat', monkey, experiment))).opt;

[coh, cor, cond, session, repeat] = deal([]);
for r = 1:n_repeats
    data = load(fullfile(InterimDir, sprintf('classifier_%s_%s_repeat%d.mat', monkey, experiment, r))).data; 
    for n = 1:n_files
        for c = 1:n_conds
            ntrial = length(data{n,c}.param.morph_level);
            coh = [coh; data{n,c}.param.morph_level];
            cor = [cor; data{n,c}.Correct{1}];
            cond = [cond; c*ones(ntrial,1)];
            session = [session; n*ones(ntrial,1)];
            repeat = [repeat; r*ones(ntrial,1)];
        end
    end
end

end