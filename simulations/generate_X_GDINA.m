function [X, X_A] = generate_X_GDINA(N, nu_true, Q, delta_true)
% This function generate data under the GDINA model given Q-matrix and model
% parameters nu_true and delta_true for sample size N


[M, K] = size(Q);
% generate multinomial counts
counts = mnrnd(N, nu_true);
X_A = zeros(N, 1);
n = 1;
% A stores empirical CDF for categories: 1,2,3,...,2^K
for a = 1:length(nu_true)
    X_A(n:(n+counts(a)-1)) = a-1;
    n = n+counts(a);
end
% permute A, size N by 1
X_A = X_A(randperm(N));
% convert each number in A (ranging from 0 to 2^K-1) to binary form, N * K
% each row of A stores the attribute profile for each individual
X_A = binary(X_A, K); % N * K

attr_combo = get_I(K, K); % (2^K-1) * K
expand_XA = prod(bsxfun(@power, reshape(X_A, [N 1 K]), reshape(attr_combo, [1 2^K-1 K])), 3); % N * (2^K-1)
expand_XA = [ones(N, 1), expand_XA]; % M matrix in the GDINA paper


%%%%%%%%%% generate GDINA parameters %%%%%%%%%%%
p_correct = expand_XA * delta_true';

X = double(rand(size(p_correct)) < p_correct);

end