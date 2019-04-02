function [is_gid, is_local_gid] = check_Theorem2(Q)

% This function checks if the sufficient conditions in Theorme 2
% for generic identifiability of Q-matrix under DINA model are satisfied

% Theorem 2 treats the case where some attribute is required by only 2
% items

% Q = [eye(3); 1 1 0; 1 0 1; 0 1 1; 1 1 1];

% Q = [eye(8);...
%      0 0 1 1 1 0 1 1;...
%      0 1 0 1 0 1 1 1;...
%      1 0 0 0 1 1 1 1;...
%      1 1 1 1 1 1 0 1];
% time0 = cputime;


% Q =[0 0 0 0 0 1 0 0;...
%     0 0 0 1 0 0 1 0;...
%     0 0 0 1 0 0 1 0;...
%     0 1 1 0 1 0 1 0;...
%     0 1 0 1 0 0 1 1;...
%     0 0 0 0 0 0 1 0;...
%     1 1 0 0 0 0 1 0;...
%     0 0 0 0 0 0 1 0;...
%     0 1 0 0 0 0 0 0;...
%     0 1 0 0 1 0 1 1;...
%     0 1 0 0 1 0 1 0;...
%     0 0 0 0 0 0 1 1;...
%     0 1 0 1 1 0 1 0;...
%     0 1 0 0 0 0 1 0;...
%     1 0 0 0 0 0 1 0;...
%     0 1 0 0 0 0 1 0;...
%     0 1 0 0 1 0 1 0;...
%     0 1 0 0 1 1 1 0;...
%     1 1 1 0 1 0 1 0;...
%     0 1 1 0 1 0 1 0];


attr_count = sum(Q, 1);


% If all attributes are required by >= 3 item, then should check Theorem 1
if all(attr_count >= 3)
    is_gid = 0;
    is_local_gid = 0;
    fprintf('\n All attributes are required by >=3 items; check Theorem 1 instead.\n')
    return
end


% If some attribute is required by <= 1 item, then not generic ID.
if any(attr_count <= 1)
    is_gid = 0;
    is_local_gid = 0;
    fprintf('\n Under DINA, this Q-matrix IS NOT generically identifiable.\n')
    return
end



attr_only2 = find(attr_count == 2);

if length(attr_only2) > 1
    
    % Check if scenario (b.2) in Theorem 2 applies: 
    % Q contains two I_K (double-complete)
    % % if the function has not yet returned here, then Condition C is
    % % violated, and Q is at most generically identifiable
    is_double_comp = check_double_complete(Q);
    if is_double_comp
        is_gid = 1;
        is_local_gid = 1;
        fprintf('\n Under DINA, this Q-matrix IS generically identifiable, because Q is double-complete and conditions in Theorem 2(b.2) are satisfied.\n')
        return
    end
    
    fprintf('\n Theorem 2 is not applicable.\n')    
    is_gid = nan;
    is_local_gid = nan;
    return
end

[J, K] = size(Q);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get indices of the two items that require attribute attr_only2 (i.e., k)
ind_2items = find( Q(:, attr_only2) == 1 );
% Two subvectors of length K-1, corresponding to the only 2 items requiring
% attribute k
two_subvec = Q(ind_2items, setdiff(1:K, attr_only2));

% find standard basis vector from the two vectors
ind_item_1attr = find( all(two_subvec == 0, 2) );
% if the two vectors are not single-attribute items,
% then none of the two subvectors are all-zero, then not generically ID.
if isempty(ind_item_1attr)
    is_gid = 0;
    is_local_gid = 0;
    fprintf('\n Under DINA, this Q-matrix is NOT generically identifiable, because Q is not complete.\n')
    return
end



% Obtain the (J-2)*(K-1) submatrix Q_sub of Q
Q_sub = Q(setdiff(1:J, ind_2items), setdiff(1:K, attr_only2));

% Check if Q_sub satisfies the necessary and sufficient conditions in
% Theorem 1
[id_Q_sub, ~, ~, ~] = check_Theorem1(Q_sub, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If both q-vectors are standard basis vectors, (i.e., v=0_{K-1} in Theorem 2)
% then apply Theorem 2(b)
if length(ind_item_1attr)==2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if scenario (b.1) applies
    if id_Q_sub
        is_gid = 1;
        is_local_gid = 1;
        fprintf('\n Under DINA, this Q-matrix IS generically identifiable, because conditions in Theorem 2(b.1) are satisfied.\n')
        return
    end
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see Theorem 2 in the paper for the definition of subvector v
% If v = 1_{K-1}, then apply Theorem 2(a)
% If the program has not yet returned here, then Q must fall in the
% scenario of part (c) of Theorem 2

% Find the item with q-vector (1, v)
ind_item_v = setdiff(ind_2items, ind_item_1attr);
v = Q(ind_item_v, setdiff(1:K, attr_only2));


if id_Q_sub && ~all(v==1)
    is_gid = nan;
    is_local_gid = 1;
    fprintf('\n Under DINA, this Q-matrix IS LOCALLY generically identifiable, because conditions in Theorem 2(c) are satisfied.\n')
    return
end


% Check if Theorem 2(a) applies

if all(v==1)
    is_gid = 0;
    is_local_gid = 0;
    fprintf('\n Under DINA, this Q-matrix IS NOT LOCALLY generically identifiable, because conditions in Theorem 2(a) are satisfied.\n')
    return
end






% mytime = cputime - time0;

end