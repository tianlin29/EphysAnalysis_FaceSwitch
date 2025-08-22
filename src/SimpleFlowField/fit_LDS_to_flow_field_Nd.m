function out = fit_LDS_to_flow_field_Nd(data, N)
% Fit a latent linear dynamical system (LDS) with linear basis of dimension N
%
% Model:
%   phi(x) ∈ R^N: zero-padded linear basis
%   dx/dt ≈ C (A phi(x) + b) = M * [phi(x); 1]
%
% Inputs:
%   data: struct with fields Xgrid, Ygrid, U, V, Cnt (all N x M)
%   N: latent dimensionality (>=2)
%
% Outputs:
%   M: 2 x (N+1) matrix = C [A b]
%   pU, pV: predicted flow fields (N x M)
%   R2: combined weighted R² of the fit

% Flatten inputs
X = [data.Xgrid(:), data.Ygrid(:)];
V = [data.U(:), data.V(:)];
W = data.Cnt(:);

valid = isfinite(W) & (W > 0) & all(isfinite(X),2) & all(isfinite(V),2);
X = X(valid,:);
V = V(valid,:);
W = W(valid);

% Linear basis phi(x): zero-padded x
X_aug = zeros(length(W), N);
X_aug(:,1:min(2,N)) = X(:,1:min(2,N));  % Use x, y in first two dims
V_aug = zeros(length(W), N);
V_aug(:,1:min(2,N)) = V(:,1:min(2,N));  % Use x, y in first two dims

% Augment with bias
% X_aug = [X_aug, ones(size(X_aug,1),1)];  % n x (N+1)

% Weighted least squares: solve V = Phi_aug * M'
WPhi = X_aug .* W;   % n x (N+1)
WV   = V_aug .* W;         % n x 2

params = WPhi \ WV;      % (N+1) x 2 matrix = C [A b]
A = params(1:2, :)';
b = params(3, :)';

% Predict full flow field
X_all = [data.Xgrid(:), data.Ygrid(:)];
Phi_all = zeros(size(X_all,1), N);
Phi_all(:,1:min(2,N)) = X_all(:,1:min(2,N));
V_pred = [Phi_all, ones(size(Phi_all,1),1)] * params;  % n x 2

% Reshape to original grid
pU = reshape(V_pred(:,1), size(data.Xgrid));
pV = reshape(V_pred(:,2), size(data.Ygrid));

% Goodness of fit (Weighted R² for each component)
V_emp = [data.U(:), data.V(:)];
W_all = data.Cnt(:);
valid_all = all(isfinite(V_emp),2) & isfinite(W_all) & (W_all > 0);

V_emp = V_emp(valid_all,:);
V_pred_valid = V_pred(valid_all,:);
W_valid = W_all(valid_all);

mean_V = sum(W_valid .* V_emp) ./ sum(W_valid);

% Squared residuals and total variance (vector norms)
residual_norm_sq = sum(W_valid .* sum((V_emp - V_pred_valid).^2, 2));
total_norm_sq    = sum(W_valid .* sum((V_emp - mean_V).^2, 2));

R2 = 1 - (residual_norm_sq / total_norm_sq);

out.A = A;
out.b = b;
out.pU = pU;
out.pV = pV;
out.R2 = R2;



