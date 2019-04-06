function [is_complete, ind_I_K] = check_complete(Q)

%%% Condition A: check if there is an item solely requiring attribute k
%%% for each attribute k=1,2,...,K
[J, K] = size(Q);
is_k_single = zeros(1, K);
ind_I_K = zeros(1, K);
for k = 1:K
    ee = zeros(1, K);
    ee(k) = 1;
    
    % absolute difference between rows of Q and ee_k
    absdiff = abs(Q-repmat(ee, J, 1));
    
    % find the indices of the single-attr-items requiring k
    find_single = find(sum(absdiff, 2) == 0);
    is_k_single(k) = ~isempty(find_single);
    if is_k_single(k)
        ind_I_K(k) = find_single(1);
    end
end
is_complete = prod(is_k_single);

end