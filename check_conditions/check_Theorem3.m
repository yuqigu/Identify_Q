function [is_atleast_2] = check_Theorem3(Q)
% This function checks if Condition C is satisfied by matrix Q (Theorem 3)

attr_count = sum(Q, 1);


if any(attr_count <= 2)
    is_atleast_2 = 0;
    fprintf('\n Some attribute is required by <2 items; a general RLCM under this Q is not generically identifiable.\n')
else
    is_atleast_2 = 1;
    fprintf('\n All attributes are required by >=2 items.\n')
end

end