function [is_double_comp] = check_double_complete(Q)

%%% This function checks if Q contains two disjoint complete submatrices
%%% i.e., check if each attr. is required by >=2 single attribute items

[J, K] = size(Q);
is_k_2single = zeros(1, K);
% single_index = zeros(1, K);
for k = 1:K
    ee = zeros(1, K);
    ee(k) = 1;
    
    % absolute difference between rows of Q and ee_k
    absdiff = abs(Q-repmat(ee, J, 1));
    
    % find the indices of the single-attr-items requiring k
    find_single = find(sum(absdiff, 2) == 0);
    is_k_2single(k) = (length(find_single)>=2);
%     if is_k_single(k)
%         single_index(k) = find_single(1);
%     end
end

is_double_comp = prod(is_k_2single);

end