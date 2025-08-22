function input = gen_input_matrix(tstamp, stim_dur, param, opt)
% function input = gen_input_matrix(tstamp, stim_dur, param, opt)

%
% input = time x trial x param
% 

def.latency = 100; %ms
def.conv_kernel = [];

opt = safeStructAssign(def, opt);

ntri = length(stim_dur);
nparam = length(param);

input = zeros(length(tstamp), ntri, nparam);

for n=1:ntri
    idx = tstamp >= opt.latency & tstamp <= (stim_dur(n) + opt.latency);
    for p = 1:nparam
        input(idx, n, p) = param{p}(n);
        if ~isempty(opt.conv_kernel)
            input(:, n, p) = nanconv(input(:, n, p), opt.conv_kernel(:), 'same');
        end
    end
end



