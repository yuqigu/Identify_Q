function [Rp] = get_rp(A, Q, c, g, I_full)
% produce R-matrix for all rows denoted by I, storing proportions of 
% positive response for each response pattern in I

[M, K] = size(Q);

ideal_resp = prod( bsxfun(@power, reshape(A, [1 size(A,1) K]), reshape(Q, [M 1 K])), 3);

prob_posi = bsxfun(@plus, g, bsxfun(@times, ideal_resp, c - g)); % M * 2^K prob of positive response for the capable
prob_zero = 1 - prob_posi; % M * 2^K prob of negative response for the capable

  
Rp = zeros(2^M, size(A, 1));
for i = 1:size(I_full, 1)
    Rp(i, :) = prod(prob_zero(I_full(i,:)==0,:),1) .* prod(prob_posi(I_full(i,:)==1,:),1);
end


end