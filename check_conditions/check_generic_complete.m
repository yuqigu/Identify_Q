function [is_gen_com] = check_generic_complete(Q)

% This function checks if the Q-matrix is generic complete
% This funciton checks generic completeness using Hall's marriage theorem

% Q = [eye(3); 1 1 0; 1 0 1; 0 1 1; 1 1 1];

% Q = ones(7, 4);
% time0 = cputime;

[~, K] = size(Q);


% check if Hall's marriage condition is satisfied
% check if for any subset S of [K], if |N(S)| >= |S|
Halls = zeros(2^K-1, 1);

% get all the 2^K subsets of [K]
all_subset = binary(1:(2^K-1), K);

% compute the number of neighbors of each subset of [K]
num_neighbor = zeros(2^K-1, 1);
for kk = 1:(2^K-1)
    % extract columns of Q corresponding to the kk-th subset of attributes
    Q_cols = Q(:, all_subset(kk,:) ~= 0);
    % take the union of items that require the attributes in kk-th subset
    num_neighbor(kk) = sum(max(Q_cols, [],  2));
    
    Halls(kk) = ( num_neighbor(kk) >= sum(all_subset(kk,:)) );
end

is_gen_com = prod(Halls);

% mytime = cputime - time0;

end
