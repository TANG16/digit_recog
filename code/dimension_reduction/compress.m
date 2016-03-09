function Y = compress(X, V, mu, k)

%% X_c: N-by-D
X_c = bsxfun(@minus, X, mu);

%% V: D-by-D -> D-by-d
V = V(:, 1:k);

%% Principle components
Y = (V' * X_c')';

end