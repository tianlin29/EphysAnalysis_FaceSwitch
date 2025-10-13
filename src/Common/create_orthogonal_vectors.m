function orthogonal_vectors = create_orthogonal_vectors(dim, num_vectors)

if num_vectors > dim
    error('Number of vectors cannot be larget than number of dimensions.');
end

random_matrix = randn(dim, num_vectors);

% get orthogonal vectors by QR decomposition
[Q, R] = qr(random_matrix);
orthogonal_vectors = Q(:, 1:num_vectors);

% make sure norm=1 
% for i = 1:num_vectors
%     orthogonal_vectors(:, i) = orthogonal_vectors(:, i) / norm(orthogonal_vectors(:, i));
% end

end