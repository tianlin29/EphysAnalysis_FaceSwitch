function roughness_score = get_roughness(u,v,x,y)

n_files = length(u);
roughness_score = nan(n_files, 1);
for n = 1:n_files
    [du_dx, du_dy] = gradient(u{n}, x{n}(1,:), y{n}(:,1));
    [dv_dx, dv_dy] = gradient(v{n}, x{n}(1,:), y{n}(:,1));
    grad_mag = sqrt(du_dx.^2 + du_dy.^2 + dv_dx.^2 + dv_dy.^2);
    
    % 2. 计算曲率
    [d2u_dx2, d2u_dy2] = gradient(du_dx, x{n}(1,:), y{n}(:,1));
    [d2v_dx2, d2v_dy2] = gradient(dv_dx, x{n}(1,:), y{n}(:,1));
    
    curvature = sqrt(d2u_dx2.^2 + d2u_dy2.^2 + d2v_dx2.^2 + d2v_dy2.^2);
    
    % 3. 归一化处理
    grad_norm = (grad_mag - min(grad_mag(:))) / (max(grad_mag(:)) - min(grad_mag(:)));
    curv_norm = (curvature - min(curvature(:))) / (max(curvature(:)) - min(curvature(:)));
    
    % 4. 组合指标（可根据需要调整权重）
    alpha = 0.6;  % 梯度权重
    beta = 0.4;   % 曲率权重
    
    roughness = alpha * grad_norm + beta * curv_norm;
    
    % 5. 转换为顺滑度评分
    roughness_score(n) = mean(roughness(:), 'omitnan');
end

end