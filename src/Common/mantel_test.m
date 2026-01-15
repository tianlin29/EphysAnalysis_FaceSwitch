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
        % get distance matrix
        D1 = pdist2(D1, D1);
        D2 = pdist2(D2, D2);
    end

    % extract lower triangular parts (excluding diagonal)
    vec_D1 = D1(tril(true(size(D1)), -1));
    vec_D2 = D2(tril(true(size(D2)), -1));

    % calculate original correlation coefficient (Pearson)
    r_obs = corr(vec_D1, vec_D2, 'Type', 'Pearson');

    % permutation test
    n = size(D1, 1);
    r_perm = zeros(n_perm, 1);
    for i = 1:n_perm
        perm_idx = randperm(n);
        D2_perm = D2(perm_idx, perm_idx); % permute both rows and columns
        vec_D2_perm = D2_perm(tril(true(size(D2_perm)), -1));
        r_perm(i) = corr(vec_D1, vec_D2_perm, 'Type', 'Pearson');
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