function [is_mono, delta] = check_monotone(Q, theta)

% This function checks if (sparsely written) item parameter matrix theta
% satisfies the monotone constraint specifed by the Q-matrix under GIDNA

% First get delta matrix from theta matrix, then check if ALL hte nonzero
% elements of delta are positive; if so, monotone constraint is satisfied

[J, K] = size(Q);
attr_combo = get_I(K, K); % (2^K-1) * K
expand_Q = [ones(J,1), prod(bsxfun(@power, reshape(Q, [J 1 K]), ...
    reshape(attr_combo, [1 2^K-1 K])), 3)];

expand_Q_logic = (expand_Q==1);

S = sum(Q, 2);
expand_S = cell(K, 1);
for j = 1:J
    expand_S{j} = [zeros(1, S(j)); get_I(S(j), S(j))];
end


MM = cell(J, 1);
inv_MM = cell(J, 1);
for j = 1:J
    MM{j} = prod(bsxfun(@power, reshape(expand_S{j}, ...
        [2^S(j) 1 S(j)]), reshape(expand_S{j}, [1 2^S(j) S(j)])), 3);
    inv_MM{j} = inv(MM{j});
end


delta = zeros(J, 2^K);
for j = 1:J
    delta(j, expand_Q_logic(j, :)) = ( inv_MM{j} * theta(j, expand_Q_logic(j, :))' )';
end

is_mono = all(delta(delta ~= 0) > 0);


end