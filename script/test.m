a = [-0.1 0.5 0.2];
b = [0.2 -0.4 -0.1];
a * b'

%%
% 示例1：创建2个10维的正交向量
dim = 10;
num_vectors = 2;
orthogonal_inputs = create_orthogonal_vectors(dim, num_vectors);

input1 = orthogonal_inputs(:, 1);
input2 = orthogonal_inputs(:, 2);

fprintf('输入向量1:\n');
disp(input1');
fprintf('输入向量2:\n');
disp(input2');

% 验证正交性
dot_product = dot(input1, input2);
fprintf('内积: %.10f\n', dot_product);

% 验证范数（应该是单位向量）
fprintf('向量1的范数: %.6f\n', norm(input1));
fprintf('向量2的范数: %.6f\n', norm(input2));

%%
function orthogonal_vectors = create_orthogonal_vectors(dim, num_vectors)
%CREATE_ORTHOGONAL_VECTORS 创建一组正交向量
%   输入参数:
%       dim - 向量维度
%       num_vectors - 需要创建的向量数量
%   输出参数:
%       orthogonal_vectors - dim × num_vectors 的正交矩阵

    % 参数验证
    if num_vectors > dim
        error('向量数量不能超过维度数');
    end
    
    % 生成随机矩阵
    random_matrix = randn(dim, num_vectors);
    
    % 使用QR分解获得正交向量
    [Q, R] = qr(random_matrix);
    
    % 提取前num_vectors列作为正交向量
    orthogonal_vectors = Q(:, 1:num_vectors);
    
    % 可选：确保向量是单位向量（归一化）
    % for i = 1:num_vectors
    %     orthogonal_vectors(:, i) = orthogonal_vectors(:, i) / norm(orthogonal_vectors(:, i));
    % end
end