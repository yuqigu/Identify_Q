function [T] = get_T_GDINA(Q, I, prob_posi)
% Given a Q-matrix and positive response probabilities, output the correct
% T-matrix

% @param Q               : m x k Q-matrix
% @param A               : ? x k attribute matrix (profiles by row)
% @param I               : ? x k cellarray of item combos by row
% @param c               : m x 1 vector c = 1-s, s = slipping parameters
% @param g               : m x 1 vector of guessing parameters
%
% @return T              : ? x k T-matrix

[~, K] = size(Q);
% A = binary(0:(2^K-1), K); % 2^K * K
% attr_combo = get_I(K, K); % (2^K-1) * K
% expand_A = prod(bsxfun(@power, reshape(A, [size(A,1) 1 K]), reshape(attr_combo, [1 2^K-1 K])), 3); 
% expand_A = [ones(size(A,1), 1), expand_A]; % M matrix in the GDINA paper
% ideal response for all single items and all attribute profiles, M * 2^K
% prob_posi = delta * expand_A';

T = zeros(size(I, 1), 2^K);
for i = 1:size(I, 1)
    T(i, :) = prod(prob_posi(I(i,:)==1,:), 1);
end


end

 