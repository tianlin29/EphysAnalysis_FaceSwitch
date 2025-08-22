function out = fit_LDS_to_flow_field(data)
% Fits an affine 2D linear dynamical system (LDS) to a 2D flow field
%
% Model: dx/dt = A * x + b
%
% Input:
%   data.Xgrid: N x M matrix of x-coordinates
%   data.Ygrid: N x M matrix of y-coordinates
%   data.U:     N x M x nslice matrix of flow in x-direction
%   data.V:     N x M x nslice matrix of flow in y-direction
%   data.Cnt:   N x M matrix of trial counts (weights)
%
% Output:
%   out{i}: i = nslice, struct with fields:
%   A:   2x2 recurrent matrix
%   b:   2x1 offset vector
%   pU:  N x M matrix of predicted flow in x-direction
%   pV:  N x M matrix of predicted flow in y-direction
%   R2:  1x2 vector of R² for U and V fits

nslice = size(data.U, 3);

for i = 1:nslice
    % Flatten grid and flow fields for slice i
    X = [data.Xgrid(:), data.Ygrid(:)];  % Position: (N*M) x 2
    U_slice = data.U(:,:,i);             % N x M
    V_slice = data.V(:,:,i);             % N x M
    V = [U_slice(:), V_slice(:)];        % Velocity: (N*M) x 2
    W = data.Cnt(:,:,i);                     % Trial counts (weights)
    W = W(:);

    % Remove invalid entries
    valid = isfinite(W) & (W > 0) & all(isfinite(X),2) & all(isfinite(V),2);
    X_valid = X(valid, :);
    V_valid = V(valid, :);
    W_valid = W(valid);

    % Augment X with a constant term for bias (affine term b)
    X_aug = [X_valid, ones(size(X_valid,1),1)];  % (N*M) x 3

    % Weighted least squares: solve V = X_aug * [A'; b']
    WX = X_aug .* W_valid;   % element-wise multiply weights
    WV = V_valid .* W_valid;

    % Solve for [A'; b']
    params = WX \ WV;   % 3 x 2 matrix: first two rows = A', last row = b'
    A = params(1:2, :)';  % 2x2
    b = params(3, :)';    % 2x1

    % Predict flow at all original grid locations
    X_all = [data.Xgrid(:), data.Ygrid(:)];
    X_all_aug = [X_all, ones(size(X_all,1),1)];
    V_pred = (X_all_aug * params);  % (N*M) x 2

    % Reshape to original grid
    pU = reshape(V_pred(:,1), size(data.Xgrid));
    pV = reshape(V_pred(:,2), size(data.Ygrid));

    % Goodness of fit (Weighted R² for each component)
    valid_all = all(isfinite(V),2) & isfinite(W) & (W > 0);
    V_emp_valid = V(valid_all,:);
    V_pred_valid = V_pred(valid_all,:);
    W_valid_all = W(valid_all);

    mean_V = sum(W_valid_all .* V_emp_valid) ./ sum(W_valid_all);

    % Squared residuals and total variance (vector norms)
    residual_norm_sq = sum(W_valid_all .* sum((V_emp_valid - V_pred_valid).^2, 2));
    total_norm_sq    = sum(W_valid_all .* sum((V_emp_valid - mean_V).^2, 2));

    R2 = 1 - (residual_norm_sq / total_norm_sq);

    out{i}.A = A;
    out{i}.b = b;
    out{i}.pU = pU;
    out{i}.pV = pV;
    out{i}.R2 = R2;
    if isfield(data, 'Z')
        out{i}.Z = data.Z(i);
    end
end

if nslice == 1
    out = out{1};
end



