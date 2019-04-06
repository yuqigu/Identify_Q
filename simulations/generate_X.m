function [X, A] = generate_X(N, p, Q, c, g)
% This function generate data under the DINA model given Q-matrix and model
% parameters p, c(=1-s), g for sample size N

% @param N               : sample size
% @param p               : 2^K x 1 vector of profile probabilities
% @param Q               : M x K item-attribute relationship
% @param c               : M x 1 vector of 1-slipping param
% @param g               : M x 1 guessing param
%
% @return X              : N x M matrix of responses
% @return A              : N x K matrix of attributes for individuals

[M, K] = size(Q);
counts = mnrnd(N, p);
A = zeros(N, 1);
n = 1;
% A stores empirical CDF for categories: 1,2,3,...,2^K
for a = 1:length(p)
    A(n:(n+counts(a)-1)) = a-1;
    n = n+counts(a);
end
% permute A, size N by 1
A = A(randperm(N));
% convert each number in A (ranging from 0 to 2^K-1) to binary form, size N by K
% each row of A stores the attribute profile for each individual
A = binary(A, K);

% ideal response matrix, size N by K
xi = prod(bsxfun(@power, reshape(A, [N 1 K]), reshape(Q, [1 M K])), 3);
% Conditional probability for each individual correctly answering each
% question, size N by K
p_correct = bsxfun(@plus, g', bsxfun(@times, xi, (c-g)'));

X = double(rand(size(p_correct)) < p_correct);

end