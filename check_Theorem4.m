function [cond_hold, item_set_1, item_set_2, item_set_3] = check_Theorem4(Q)

% This function checks if the sufficient conditions D and E in Theorem 4 for
% generic identifiability of Q under the GDINA model are satisfied
%%%%
% This funciton also checks if the necessary conditions in Theorem 3 are
% satisfied (each attribute is required by >=2 items)
%%%%
% This funciton also checks if the necessary conditions in Theorem 5 are
% satisfied (Q is generically complete)

% In particular, the function searches for two disjoint K*K submat Q1, Q2 
% that are generically complete; and the remaining Qprime has no zero column


% % The following are several examples

% Q = [0 0 1; 1 0 1; 1 1 0; 0 1 1; 1 0 1; 1 1 0; 1 1 1]; % size &*3

% % The following Q has size 10 * 4
% Q = [1 0 1 0; 1 1 0 0; 0 0 1 1; 0 1 0 1;...
%     1 0 0 1; 0 1 1 0; 0 1 1 1; 1 0 0 1;...
%     1 1 0 0; 0 0 1 1];
% Q = Q([1 5 9 2 6 10 4 8 7 3], :);
% % The search takes <0.3 second.

% % Fraction subtraction data Q-matrix
% % after deleting the 6th column and the 1st and 18th rows of it
% % the remaining submatrix of Q satisfies the Generic ID. conditions
% Q =[0 0 0 1 0 1 1 0;...
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
% Q = Q([2:17 19:20], :);
% Q = Q(:, [1:5 7:8]); 
% % this is a 18*7 submatrix of Q that is generically id.
% % the search takes 96 seconds.


% time0 = cputime;

[J, K] = size(Q);

attr_count = sum(Q, 1);

% If some attributes is required by <=2 item, then not generic ID. by
% Theorem 3
if any(attr_count <= 2)
    cond_hold = 0;
    item_set_1=[]; item_set_2=[]; item_set_3=[];
    fprintf('\n Under a general RLCM, this Q-matrix IS NOT generically identifiable by Theorem 3, because some attribute is required by <=2 items.\n')
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if Q is generically complete; if not, then not generically identifiable.
[is_gen_com] = check_generic_complete(Q);
if ~is_gen_com
    cond_hold = 0;
    item_set_1=[]; item_set_2=[]; item_set_3=[]; 
    fprintf('\n Under a general RLCM, this Q-matrix IS NOT generically identifiable by Theorem 5, since Q is not generically complete.\n')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if J>=2K+1; if not, conditions not satisfied
if J < (2*K+1)
    cond_hold = 0;
    item_set_1=[]; item_set_2=[]; item_set_3=[]; 
    fprintf('\n Conditions in Theorem 5 are not satisfied.\n')
    return
end






% compute the number of times each attr. is required in Q
attr_times = sum(Q, 1);

% if some attribute is required by <3 items, not generically identifiable.
if any(attr_times<3)
    cond_hold = 0;
    item_set_1=[]; item_set_2=[]; item_set_3=[]; 
    % fprintf('\n Under a general RLCM, this Q-matrix is not generically identifiable by Theorem 3, since some attribute is required by only two items.\n')
    return
end

% an indicator of whether the generic identifiability conditions hold
cond_hold = 0;

% get all K*K submatrices of Q
% each row of all_submat_K is a K-combination from 1:J
all_submat_K1 = nchoosek(1:J, K);
for ii = 1:size(all_submat_K1, 1)
    % ii is the index for a K*K submat
    Q1 = Q(all_submat_K1(ii,:), :);
    Q1_GC = check_generic_complete(Q1);
    
    % fprintf('Number %d search for Q1 started\n', ii);
    % if Q1 is not Gen. Complete, then go to next Q1
    if Q1_GC
        item_set_1 = all_submat_K1(ii,:);
    else
        continue
    end
    
    % if Q1 is Gen. Complete, search for Q2
    ind_remain = setdiff(1:J, all_submat_K1(ii,:));

    all_submat_K2 = nchoosek(ind_remain, K);

    for jj = 1:size(all_submat_K2, 1)
        Q2 = Q(all_submat_K2(jj,:), :);
        Q2_GC = check_generic_complete(Q2);
        
        fprintf('In numer %d search for Q1, number %d search for Q2 started\n',...
            ii, jj);

        if Q2_GC
            item_set_2 = all_submat_K2(jj,:);
        else
            continue
        end
        
        % if two Gen. Complete Q1 and Q2 are found, then find indices
        % of the remaining Q' part, and check Condition E
        ind_last = setdiff(ind_remain, all_submat_K2(jj,:));
        
        if all(sum(Q(ind_last,:), 1) >= 1)
            item_set_3 = ind_last;
        else
            continue
        end
        
        cond_hold = 1;
        fprintf('\n Generic identifiability conditions D and E satisfied! Q is generically identifiable!\n');
        break
    end
    
    if cond_hold == 1
        break
    end
    
end

fprintf('The first set of item indices (Q1):\n')
item_set_1

fprintf('The first set of item indices (Q2):\n')
item_set_2

fprintf('The first set of item indices (Qprime):\n')
item_set_3

fprintf('The Q1 is:\n')
Q(item_set_1,:)

fprintf('The Q2 is:\n')
Q(item_set_2,:)

fprintf('The Qprime is:\n')
Q(item_set_3,:)


% mytime = cputime - time0;

end